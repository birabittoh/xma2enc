// XMA2 encoder. The bitstream is WMA Pro with the XMA2 packet layer, always
// coded as one 512 sample block per frame. See docs/format.md.
#include "xma2enc/xma2enc.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <unordered_map>

#include "tables.h"

namespace xma2enc {
namespace {

constexpr int N = kFrameSamples;
constexpr int kMaxBands = 29;
constexpr int kMaxFrameBits = 12000;
constexpr double kLrScale = 181.0 / 128.0;  // decoder gain on non M/S bands
constexpr double kGamma = 0.25;              // noise follows band energy^gamma
constexpr double kPi = 3.14159265358979323846;

struct BitWriter {
  std::vector<uint8_t> bytes;
  size_t bits = 0;

  void Put(uint32_t value, int count) {
    for (int i = count - 1; i >= 0; --i) {
      if ((bits & 7) == 0) bytes.push_back(0);
      if ((value >> i) & 1) bytes.back() |= uint8_t(0x80 >> (bits & 7));
      ++bits;
    }
  }
  void Append(const BitWriter& other) {
    for (size_t i = 0; i < other.bits; ++i)
      Put((other.bytes[i >> 3] >> (7 - (i & 7))) & 1, 1);
  }
};

struct Code {
  uint32_t code = 0;
  uint8_t len = 0;
};

// Canonical assignment used by the decoder: entries in table order, each code
// is the previous one plus one at its own length.
template <typename LenFn>
std::vector<Code> BuildCodes(int count, LenFn len_of) {
  std::vector<Code> codes(count);
  uint64_t acc = 0;
  for (int i = 0; i < count; ++i) {
    int len = len_of(i);
    if (len <= 0) continue;
    codes[i].len = uint8_t(len);
    codes[i].code = uint32_t(acc >> (32 - len));
    acc += uint64_t(1) << (32 - len);
  }
  return codes;
}

struct RunLevelTable {
  Code esc, eob;
  std::unordered_map<uint32_t, Code> map;  // (run << 16) | level
};

struct Tables {
  std::vector<Code> sf_delta;  // index delta + 60
  RunLevelTable rl[2];
  std::vector<Code> vec4;      // by symbol
  std::vector<Code> vec2;
  std::vector<Code> vec1;

  Tables() {
    auto sf = BuildCodes(HUFF_SCALE_SIZE, [](int i) { return int(scale_table[i][1]); });
    sf_delta.assign(HUFF_SCALE_SIZE, Code());
    for (int i = 0; i < HUFF_SCALE_SIZE; ++i) sf_delta[scale_table[i][0]] = sf[i];

    auto c0 = BuildCodes(HUFF_COEF0_SIZE, [](int i) { return int(coef0_lens[i]); });
    for (int i = 0; i < HUFF_COEF0_SIZE; ++i)
      AddRl(rl[0], c0[i], coef0_syms[i], coef0_run, coef0_level);
    auto c1 = BuildCodes(HUFF_COEF1_SIZE, [](int i) { return int(coef1_table[i][1]); });
    for (int i = 0; i < HUFF_COEF1_SIZE; ++i)
      AddRl(rl[1], c1[i], coef1_table[i][0], coef1_run, coef1_level);

    auto v4 = BuildCodes(HUFF_VEC4_SIZE, [](int i) { return int(vec4_lens[i]); });
    vec4.assign(65538, Code());
    for (int i = 0; i < HUFF_VEC4_SIZE; ++i) vec4[vec4_syms[i]] = v4[i];
    auto v2 = BuildCodes(HUFF_VEC2_SIZE, [](int i) { return int(vec2_table[i][1]); });
    vec2.assign(512, Code());
    for (int i = 0; i < HUFF_VEC2_SIZE; ++i) vec2[vec2_table[i][0]] = v2[i];
    auto v1 = BuildCodes(HUFF_VEC1_SIZE, [](int i) { return int(vec1_table[i][1]); });
    vec1.assign(512, Code());
    for (int i = 0; i < HUFF_VEC1_SIZE; ++i) vec1[vec1_table[i][0]] = v1[i];
  }

  static void AddRl(RunLevelTable& t, Code c, int sym, const uint16_t* run,
                    const float* level) {
    if (sym == 0) {
      t.esc = c;
    } else if (sym == 1) {
      t.eob = c;
    } else {
      uint32_t key = (uint32_t(run[sym]) << 16) | uint32_t(level[sym]);
      auto it = t.map.find(key);
      if (it == t.map.end() || it->second.len > c.len) t.map[key] = c;
    }
  }
};

const Tables& T() {
  static const Tables tables;
  return tables;
}

void PutCode(BitWriter& w, const Code& c) { w.Put(c.code, c.len); }

void PutLargeVal(BitWriter& w, uint32_t v) {
  if (v < (1u << 8)) {
    w.Put(0, 1);
    w.Put(v, 8);
  } else if (v < (1u << 16)) {
    w.Put(2, 2);
    w.Put(v, 16);
  } else if (v < (1u << 24)) {
    w.Put(6, 3);
    w.Put(v, 24);
  } else {
    w.Put(7, 3);
    w.Put(v, 31);
  }
}

int BitLen(uint32_t v) {
  int n = 0;
  while (v) {
    ++n;
    v >>= 1;
  }
  return n;
}

struct BandLayout {
  int count = 0;
  int offsets[kMaxBands + 1] = {0};
};

BandLayout MakeBands(uint32_t rate) {
  int r = rate > 44100 ? 48000 : rate > 32000 ? 44100 : rate > 24000 ? 32000 : 24000;
  BandLayout b;
  int band = 1;
  b.offsets[0] = 0;
  for (int x = 0; x < kMaxBands - 1 && b.offsets[band - 1] < N; ++x) {
    int offset = (N * 2 * critical_freq[x]) / r + 2;
    offset &= ~3;
    if (offset > b.offsets[band - 1]) b.offsets[band++] = offset;
    if (offset >= N) break;
  }
  b.offsets[band - 1] = N;
  b.count = band - 1;
  return b;
}

// Writes quantized coefficients of one channel: table bit, vector part, run
// level part.
void PutCoefs(BitWriter& w, const int* q, int nv, int table) {
  const Tables& t = T();
  w.Put(table, 1);
  auto mag = [&](int i) { return uint32_t(std::abs(q[i])); };
  auto sign = [&](int i) { w.Put(q[i] < 0 ? 1 : 0, 1); };
  int i = 0;
  for (; i + 3 < nv; i += 4) {
    uint32_t m[4] = {mag(i), mag(i + 1), mag(i + 2), mag(i + 3)};
    bool small = m[0] < 16 && m[1] < 16 && m[2] < 16 && m[3] < 16;
    uint32_t sym = small ? ((m[0] << 12) | (m[1] << 8) | (m[2] << 4) | m[3]) + 1 : 0;
    if (small && T().vec4[sym].len) {
      PutCode(w, t.vec4[sym]);
    } else {
      PutCode(w, t.vec4[0]);
      for (int p = 0; p < 4; p += 2) {
        uint32_t a = m[p], b = m[p + 1];
        uint32_t s2 = (a < 16 && b < 16) ? ((a << 4) | b) + 1 : 0;
        if (s2 && t.vec2[s2].len) {
          PutCode(w, t.vec2[s2]);
        } else {
          PutCode(w, t.vec2[0]);
          for (uint32_t v : {a, b}) {
            if (v < HUFF_VEC1_SIZE - 1) {
              PutCode(w, t.vec1[v]);
            } else {
              PutCode(w, t.vec1[HUFF_VEC1_SIZE - 1]);
              PutLargeVal(w, v - (HUFF_VEC1_SIZE - 1));
            }
          }
        }
      }
    }
    for (int k = 0; k < 4; ++k)
      if (m[k]) sign(i + k);
  }

  const RunLevelTable& rl = t.rl[table];
  const int esc_len = BitLen(N - 1);
  int pos = nv;
  while (pos < N) {
    int p = pos;
    while (p < N && q[p] == 0) ++p;
    if (p >= N) {
      PutCode(w, rl.eob);
      break;
    }
    uint32_t run = uint32_t(p - pos), level = mag(p);
    auto it = rl.map.find((run << 16) | level);
    if (it != rl.map.end()) {
      PutCode(w, it->second);
    } else {
      PutCode(w, rl.esc);
      PutLargeVal(w, level);
      if (run == 0) {
        w.Put(0, 1);
      } else if (run <= 4) {
        w.Put(2, 2);
        w.Put(run - 1, 2);
      } else {
        w.Put(6, 3);
        w.Put(run - 4, esc_len);
      }
    }
    sign(p);
    pos = p + 1;
    if (pos >= N) break;
  }
}

struct ChannelPlan {
  int q[N];
  int sf[kMaxBands];      // transmitted scale factors
  int sf_step = 2;
  int quant_step = 90;
  bool transmit = false;
  int nv = 0;
  int table = 0;
  BitWriter coefs;
};

struct Frame {
  BitWriter body;  // without length prefix, pad bit and trailer
  uint32_t bits() const { return uint32_t(body.bits) + 17; }
};

int Quantize(double x, double step) {
  double a = std::fabs(x) / step;
  int q = int(a + 0.40);
  return x < 0 ? -q : q;
}

// Encodes one frame of MDCT coefficients (decoder scale) for 1 or 2 channels.
Frame EncodeFrame(const float* X[2], int channels, const BandLayout& bl, double g,
                  bool fast) {
  const int nb = bl.count;
  double target[2][kMaxBands];
  for (int c = 0; c < channels; ++c) {
    for (int b = 0; b < nb; ++b) {
      double e = 0;
      for (int k = bl.offsets[b]; k < bl.offsets[b + 1]; ++k) e += double(X[c][k]) * X[c][k];
      e /= (bl.offsets[b + 1] - bl.offsets[b]);
      target[c][b] = g * std::pow(e + 1.0, kGamma * 0.5);
    }
  }

  // Channel coding (L/R or per band M/S).
  static thread_local std::vector<float> coef_store;
  coef_store.assign(size_t(2) * N, 0.f);
  float* C[2] = {coef_store.data(), coef_store.data() + N};
  double step_target[2][kMaxBands];
  bool ms[kMaxBands] = {false};
  int ms_count = 0;
  if (channels == 2) {
    for (int b = 0; b < nb; ++b) {
      double tl = target[0][b], tr = target[1][b], tm = std::min(tl, tr) / kLrScale;
      double cost_lr = 0, cost_ms = 0;
      for (int k = bl.offsets[b]; k < bl.offsets[b + 1]; ++k) {
        double l = X[0][k], r = X[1][k];
        cost_lr += std::log2(1 + std::fabs(l) / tl) + std::log2(1 + std::fabs(r) / tr);
        cost_ms += std::log2(1 + std::fabs((l + r) * 0.5) / tm) +
                   std::log2(1 + std::fabs((r - l) * 0.5) / tm);
      }
      ms[b] = cost_ms < cost_lr;
      ms_count += ms[b];
    }
    for (int b = 0; b < nb; ++b) {
      for (int k = bl.offsets[b]; k < bl.offsets[b + 1]; ++k) {
        double l = X[0][k], r = X[1][k];
        if (ms_count == 0) {
          C[0][k] = float(l);
          C[1][k] = float(r);
        } else if (ms[b]) {
          C[0][k] = float((l + r) * 0.5);
          C[1][k] = float((r - l) * 0.5);
        } else {
          C[0][k] = float(l / kLrScale);
          C[1][k] = float(r / kLrScale);
        }
      }
      double scale = ms_count == 0 ? 1.0 : 1.0 / kLrScale;
      if (ms_count && ms[b]) {
        step_target[0][b] = step_target[1][b] = std::min(target[0][b], target[1][b]) * scale;
      } else {
        step_target[0][b] = target[0][b] * scale;
        step_target[1][b] = target[1][b] * scale;
      }
    }
  } else {
    for (int k = 0; k < N; ++k) C[0][k] = X[0][k];
    for (int b = 0; b < nb; ++b) step_target[0][b] = target[0][b];
  }

  ChannelPlan plan[2];
  bool any = false;
  for (int c = 0; c < channels; ++c) {
    ChannelPlan& p = plan[c];
    p.sf_step = 2;
    int s[kMaxBands];
    for (int b = 0; b < nb; ++b) {
      double l = 20.0 * std::log10(std::max(step_target[c][b], 1.0));
      s[b] = int(std::lround(l / p.sf_step));
      if (b) s[b] = std::clamp(s[b], s[b - 1] - 60, s[b - 1] + 60);
    }
    bool nonzero = false;
    bool band_nz[kMaxBands];
    for (int b = 0; b < nb; ++b) {
      double step = std::pow(10.0, s[b] * p.sf_step / 20.0);
      band_nz[b] = false;
      for (int k = bl.offsets[b]; k < bl.offsets[b + 1]; ++k) {
        p.q[k] = Quantize(C[c][k], step);
        if (p.q[k]) band_nz[b] = true;
      }
      nonzero |= band_nz[b];
    }
    p.transmit = nonzero;
    any |= nonzero;
    if (!nonzero) {
      for (int b = 0; b < nb; ++b) s[b] = s[0];
    } else {
      for (int b = 1; b < nb; ++b)
        if (!band_nz[b]) s[b] = s[b - 1];
      if (!band_nz[0]) {
        int f = 1;
        while (f < nb && !band_nz[f]) ++f;
        if (f < nb) s[0] = s[f];
      }
    }
    int maxs = *std::max_element(s, s + nb);
    p.quant_step = maxs * p.sf_step;
    for (int b = 0; b < nb; ++b) p.sf[b] = s[b] - s[0] + 45 / p.sf_step;
    // Re-clamp the first delta (starts from 45 / step).
  }

  BitWriter w;
  w.Put(1, 1);  // fixed channel layout
  w.Put(0, 1);  // one subframe of 512
  if (channels == 2) w.Put(0, 1);  // post processing
  w.Put(0, 8);  // drc gain
  w.Put(0, 1);  // trim
  w.Put(0, 1);  // extended header
  w.Put(0, 1);  // reserved
  if (channels == 2) {
    w.Put(0, 1);  // channel transform bit
    if (ms_count == 0) {
      w.Put(2, 2);  // no transform
    } else {
      w.Put(0, 1);
      if (ms_count == nb) {
        w.Put(1, 1);
      } else {
        w.Put(0, 1);
        for (int b = 0; b < nb; ++b) w.Put(ms[b], 1);
      }
    }
  }
  for (int c = 0; c < channels; ++c) w.Put(plan[c].transmit, 1);
  if (!any) {
    w.Put(0, 0);
    Frame f;
    f.body = w;
    return f;
  }

  // Choose the vector length and run level table per channel.
  for (int c = 0; c < channels; ++c) {
    ChannelPlan& p = plan[c];
    if (!p.transmit) continue;
    int last = 0;
    for (int k = 0; k < N; ++k)
      if (p.q[k]) last = k + 1;
    int max_nv = (last + 3) & ~3;
    int cand[48], nc = 0;
    if (fast) {
      cand[nc++] = max_nv;
      cand[nc++] = (max_nv / 2) & ~3;
    } else {
      for (int nv = 0; nv < max_nv; nv += 16) cand[nc++] = nv;
      cand[nc++] = max_nv;
    }
    size_t best = SIZE_MAX;
    BitWriter tmp;
    for (int ci = 0; ci < nc; ++ci) {
      for (int table = 0; table < 2; ++table) {
        tmp = BitWriter();
        PutCoefs(tmp, p.q, cand[ci], table);
        if (tmp.bits < best) {
          best = tmp.bits;
          p.nv = cand[ci];
          p.table = table;
          p.coefs = tmp;
        }
      }
    }
  }

  w.Put(1, 1);  // vector length transmitted
  for (int c = 0; c < channels; ++c) w.Put(plan[c].transmit ? plan[c].nv / 4 : 0, 8);

  int base = plan[0].quant_step;
  if (channels == 2) {
    // The non transmitting channel follows the other one.
    if (!plan[0].transmit) plan[0].quant_step = plan[1].quant_step;
    if (!plan[1].transmit) plan[1].quant_step = plan[0].quant_step;
    base = std::min(plan[0].quant_step, plan[1].quant_step);
  }
  int qd = base - 90;
  if (qd > -32 && qd < 31) {
    w.Put(uint32_t(qd) & 63, 6);
  } else {
    int e;
    if (qd >= 31) {
      w.Put(31, 6);
      e = qd - 31;
    } else {
      w.Put(32, 6);
      e = -32 - qd;
    }
    while (e >= 31) {
      w.Put(31, 5);
      e -= 31;
    }
    w.Put(e, 5);
  }
  if (channels == 2) {
    int d0 = plan[0].quant_step - base, d1 = plan[1].quant_step - base;
    int md = std::max(d0, d1);
    int len = md > 1 ? BitLen(uint32_t(md - 1)) : 0;
    if (len > 7) len = 7;
    w.Put(len, 3);
    for (int d : {d0, d1}) {
      if (d > 0) {
        w.Put(1, 1);
        if (len) w.Put(std::min(d - 1, (1 << len) - 1), len);
      } else {
        w.Put(0, 1);
      }
    }
  }
  for (int c = 0; c < channels; ++c) {
    const ChannelPlan& p = plan[c];
    w.Put(p.sf_step - 1, 2);
    int prev = 45 / p.sf_step;
    for (int b = 0; b < nb; ++b) {
      PutCode(w, T().sf_delta[p.sf[b] - prev + 60]);
      prev = p.sf[b];
    }
  }
  for (int c = 0; c < channels; ++c)
    if (plan[c].transmit) w.Append(plan[c].coefs);

  Frame f;
  f.body = std::move(w);
  return f;
}

// Forward MDCT (window 2N, hop N) in the decoder's coefficient scale:
// X[k] = sum x[n] w[n] cos(pi / N (n + 1/2 + N/2) (k + 1/2)). Folded into a
// DCT-IV, which runs as a complex FFT of N / 2 points.
class Mdct {
 public:
  Mdct() : win_(2 * N), tw_(N / 2), post_(N / 2), fft_tw_(N / 4) {
    for (int n = 0; n < 2 * N; ++n) win_[n] = float(std::sin(kPi * (n + 0.5) / (2 * N)));
    for (int n = 0; n < N / 2; ++n) {
      double a = -kPi * (4 * n + 1) / (4.0 * N);
      tw_[n] = {std::cos(a), std::sin(a)};
      post_[n] = {std::cos(-kPi * n / N), std::sin(-kPi * n / N)};
    }
    for (int n = 0; n < N / 4; ++n) {
      double a = -2 * kPi * n / (N / 2);
      fft_tw_[n] = {std::cos(a), std::sin(a)};
    }
  }

  void Run(const float* in, float* out) const {
    constexpr int H = N / 2;
    double u[N];
    // Quarters a, b, c, d of the windowed input.
    auto x = [&](int n) { return double(in[n]) * win_[n]; };
    for (int n = 0; n < H; ++n) {
      u[n] = -x(3 * H - 1 - n) - x(3 * H + n);  // -c_r - d
      u[H + n] = x(n) - x(N - 1 - n);           // a - b_r
    }
    std::complex<double> t[H];
    for (int n = 0; n < H; ++n) t[n] = std::complex<double>(u[2 * n], u[N - 1 - 2 * n]) * tw_[n];
    Fft(t);
    for (int k = 0; k < H; ++k) {
      std::complex<double> y = t[k] * post_[k];
      out[2 * k] = float(y.real());
      out[N - 1 - 2 * k] = float(-y.imag());
    }
  }

 private:
  void Fft(std::complex<double>* a) const {
    constexpr int H = N / 2;
    for (int i = 1, j = 0; i < H; ++i) {
      int bit = H >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= H; len <<= 1) {
      int step = H / len;
      for (int i = 0; i < H; i += len)
        for (int k = 0; k < len / 2; ++k) {
          std::complex<double> w = fft_tw_[k * step];
          std::complex<double> p = a[i + k], q = a[i + k + len / 2] * w;
          a[i + k] = p + q;
          a[i + k + len / 2] = p - q;
        }
    }
  }

  std::vector<float> win_;
  std::vector<std::complex<double>> tw_, post_, fft_tw_;
};

}  // namespace

bool Encode(const int16_t* pcm, size_t samples, uint32_t channels,
            uint32_t sample_rate, const Options& options, Result* out) {
  if (!out || !pcm || channels < 1 || channels > 2) return false;
  if (sample_rate == 0 || sample_rate > 48000) return false;
  if (options.block_size == 0 || options.block_size % kPacketSize) return false;

  const BandLayout bl = MakeBands(sample_rate);
  const size_t frames = (samples + kDecoderDelay + N - 1) / N;

  // Analysis.
  std::vector<float> coefs(frames * channels * N);
  {
    Mdct mdct;
    std::vector<float> sig(channels * (frames + 1) * N, 0.f);
    for (size_t i = 0; i < samples; ++i)
      for (uint32_t c = 0; c < channels; ++c)
        sig[c * (frames + 1) * N + kDecoderDelay + i] = pcm[i * channels + c];
    for (size_t f = 0; f < frames; ++f)
      for (uint32_t c = 0; c < channels; ++c)
        mdct.Run(&sig[c * (frames + 1) * N + f * N], &coefs[(f * channels + c) * N]);
  }

  // The rate search sizes every `stride`th frame and scales up, which is
  // accurate on long input and keeps the passes cheap.
  const size_t stride = frames >= 400 ? 8 : frames >= 100 ? 4 : 1;
  auto encode_all = [&](double g, std::vector<Frame>* store, bool fast, size_t step = 1) {
    uint64_t total = 0;
    size_t count = 0;
    if (store) store->resize(frames);
    for (size_t f = 0; f < frames; f += step) {
      ++count;
      const float* X[2] = {&coefs[(f * channels) * N],
                           channels > 1 ? &coefs[(f * channels + 1) * N] : nullptr};
      double gg = g;
      Frame fr = EncodeFrame(X, int(channels), bl, gg, fast);
      for (int tries = 0; fr.bits() > kMaxFrameBits && tries < 40; ++tries) {
        gg *= 1.2;
        fr = EncodeFrame(X, int(channels), bl, gg, fast);
      }
      total += fr.bits();
      if (store) (*store)[f] = std::move(fr);
    }
    return step == 1 ? total : uint64_t(double(total) * double(frames) / double(count));
  };

  // Bisection of the noise level for the target size.
  const double seconds = double(samples) / sample_rate;
  const double target_bits = options.bitrate * seconds * 0.97;
  double lo = 0.05, hi = 200.0;
  for (int i = 0; i < 9; ++i) {
    double mid = std::sqrt(lo * hi);
    if (double(encode_all(mid, nullptr, true, stride)) > target_bits)
      lo = mid;
    else
      hi = mid;
  }
  // The fast search overestimates sizes; settle on the full coder with a
  // secant step on log size over log noise level.
  std::vector<Frame> store;
  double g = hi;
  double g_prev = 0, bits_prev = 0;
  for (int i = 0; i < 4; ++i) {
    double bits = double(encode_all(g, nullptr, false, stride));
    double err = bits / target_bits;
    if (err <= 1.0 && err > 0.985) break;
    double slope = -1.2;
    if (g_prev > 0 && std::fabs(std::log(g / g_prev)) > 1e-6) {
      double sl = std::log(bits / bits_prev) / std::log(g / g_prev);
      if (sl < -0.3 && sl > -4) slope = sl;
    }
    g_prev = g;
    bits_prev = bits;
    g *= std::pow(target_bits * 0.99 / bits, 1.0 / slope);
  }
  // One exact pass, since the sampled sizes carry some error.
  bool verified = false;
  for (int i = 0; i < 3 && stride > 1; ++i) {
    const double bits = double(encode_all(g, &store, false));
    if (bits <= target_bits * 1.003 && bits > target_bits * 0.97) {
      verified = true;
      break;
    }
    g *= std::pow(target_bits * 0.985 / bits, 1.0 / -1.2);
  }
  if (!verified) {
    if (double(encode_all(g, nullptr, false)) > target_bits * 1.003) g *= 1.03;
    encode_all(g, &store, false);
  }


  // Packetize. Frames never cross a block boundary.
  const uint32_t packets_per_block = options.block_size / kPacketSize;
  const size_t payload_bits = size_t(kPacketSize) * 8 - 32;
  struct Packet {
    BitWriter payload;
    uint32_t count = 0;   // frames starting here
    uint32_t first = 0;   // bit offset of the first one
  };
  std::vector<Packet> pk(1);
  size_t cur = 0;
  for (size_t f = 0; f < frames; ++f) {
    const Frame& fr = store[f];
    BitWriter full;
    full.Put(fr.bits(), 15);
    full.Append(fr.body);
    full.Put(1, 1);
    full.Put(0, 1);  // trailer, patched below

    if (pk[cur].payload.bits == payload_bits) {
      ++cur;
      pk.emplace_back();
    }
    size_t block_end = (cur / packets_per_block + 1) * packets_per_block;
    size_t room = (block_end - cur) * payload_bits - pk[cur].payload.bits;
    if (full.bits > room) {
      cur = block_end;
      pk.resize(cur + 1);
    }
    if (pk[cur].count++ == 0) pk[cur].first = uint32_t(pk[cur].payload.bits);
    for (size_t i = 0; i < full.bits; ++i) {
      if (pk[cur].payload.bits == payload_bits) {
        ++cur;
        pk.emplace_back();
        pk[cur].first = uint32_t(full.bits - i);  // carried over bits
      }
      pk[cur].payload.Put((full.bytes[i >> 3] >> (7 - (i & 7))) & 1, 1);
    }
    // The decoder keeps reading frames from a packet while the trailer is
    // set, so it is set only when the next frame starts in this packet.
    if (f + 1 < frames) {
      size_t end_room = (cur / packets_per_block + 1) * packets_per_block;
      size_t room_left = (end_room - cur) * payload_bits - pk[cur].payload.bits;
      bool same_packet = pk[cur].payload.bits < payload_bits && store[f + 1].bits() <= room_left;
      if (same_packet) {
        size_t last_bit = pk[cur].payload.bits - 1;
        pk[cur].payload.bytes[last_bit >> 3] |= uint8_t(0x80 >> (last_bit & 7));
      }
    }
  }

  std::vector<uint8_t>& dst = out->packets;
  dst.assign(pk.size() * kPacketSize, 0);
  for (size_t i = 0; i < pk.size(); ++i) {
    uint8_t* d = &dst[i * kPacketSize];
    uint32_t h = (pk[i].count << 26) | ((pk[i].first & 0x7fff) << 11) | (1u << 8);
    d[0] = uint8_t(h >> 24);
    d[1] = uint8_t(h >> 16);
    d[2] = uint8_t(h >> 8);
    d[3] = uint8_t(h);
    std::copy(pk[i].payload.bytes.begin(), pk[i].payload.bytes.end(), d + 4);
  }

  out->channels = channels;
  out->sample_rate = sample_rate;
  out->frames = uint32_t(frames);
  out->samples = uint32_t(samples);
  out->block_size = options.block_size;
  out->block_count = uint32_t((pk.size() + packets_per_block - 1) / packets_per_block);
  out->seek_table.clear();
  uint64_t cumulative = 0;
  for (uint32_t blk = 0; blk < out->block_count; ++blk) {
    for (size_t i = blk * packets_per_block; i < std::min(pk.size(), size_t(blk + 1) * packets_per_block); ++i)
      cumulative += uint64_t(pk[i].count) * N;
    out->seek_table.push_back(uint32_t(cumulative));
  }
  return true;
}

}  // namespace xma2enc
