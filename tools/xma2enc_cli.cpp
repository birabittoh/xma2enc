// wav in, XMA2 out: xma2enc_cli in.wav out.xma2.wav [bitrate] [block_size]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "xma2enc/xma2enc.h"

static void Le16(std::vector<uint8_t>& v, uint32_t x) { v.push_back(x); v.push_back(x >> 8); }
static void Le32(std::vector<uint8_t>& v, uint32_t x) { Le16(v, x); Le16(v, x >> 16); }

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s in.wav out.xma2.wav [bitrate] [block_size]\n", argv[0]);
    return 2;
  }
  FILE* f = std::fopen(argv[1], "rb");
  if (!f) return 1;
  std::vector<uint8_t> in;
  uint8_t buf[65536];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) in.insert(in.end(), buf, buf + n);
  std::fclose(f);

  uint32_t channels = 0, rate = 0, bits = 0;
  const uint8_t* data = nullptr;
  size_t data_size = 0;
  for (size_t p = 12; p + 8 <= in.size();) {
    uint32_t size = in[p + 4] | in[p + 5] << 8 | in[p + 6] << 16 | uint32_t(in[p + 7]) << 24;
    if (!std::memcmp(&in[p], "fmt ", 4)) {
      channels = in[p + 10] | in[p + 11] << 8;
      rate = in[p + 12] | in[p + 13] << 8 | in[p + 14] << 16 | uint32_t(in[p + 15]) << 24;
      bits = in[p + 22] | in[p + 23] << 8;
    } else if (!std::memcmp(&in[p], "data", 4)) {
      data = &in[p + 8];
      data_size = std::min<size_t>(size, in.size() - p - 8);
      break;
    }
    p += 8 + size + (size & 1);
  }
  if (!data || bits != 16) {
    std::fprintf(stderr, "need 16 bit PCM wav\n");
    return 1;
  }

  xma2enc::Options opt;
  if (argc > 3) opt.bitrate = uint32_t(std::atoi(argv[3]));
  if (argc > 4) opt.block_size = uint32_t(std::strtoul(argv[4], nullptr, 0));
  size_t samples = data_size / (2 * channels);
  std::vector<int16_t> pcm(samples * channels);
  std::memcpy(pcm.data(), data, pcm.size() * 2);

  xma2enc::Result r;
  if (!xma2enc::Encode(pcm.data(), samples, channels, rate, opt, &r)) {
    std::fprintf(stderr, "encode failed\n");
    return 1;
  }

  std::vector<uint8_t> out;
  auto tag = [&](const char* t) { out.insert(out.end(), t, t + 4); };
  tag("RIFF");
  Le32(out, 0);
  tag("WAVE");
  tag("fmt ");
  Le32(out, 18 + 34);
  Le16(out, 0x166);
  Le16(out, channels);
  Le32(out, rate);
  Le32(out, rate * 4);
  Le16(out, xma2enc::kPacketSize);
  Le16(out, 16);
  Le16(out, 34);
  Le16(out, 1);                             // streams
  Le32(out, channels == 1 ? 4 : 3);         // channel mask
  Le32(out, r.frames * xma2enc::kFrameSamples);
  Le32(out, r.block_size);
  Le32(out, 0);
  Le32(out, uint32_t(samples));
  Le32(out, 0);
  Le32(out, 0);
  out.push_back(0);
  out.push_back(4);
  Le16(out, r.block_count);
  tag("data");
  Le32(out, uint32_t(r.packets.size()));
  out.insert(out.end(), r.packets.begin(), r.packets.end());
  uint32_t riff = uint32_t(out.size() - 8);
  std::memcpy(&out[4], &riff, 4);

  f = std::fopen(argv[2], "wb");
  std::fwrite(out.data(), 1, out.size(), f);
  std::fclose(f);
  std::printf("frames %u packets %zu blocks %u kbps %.1f\n", r.frames,
              r.packets.size() / xma2enc::kPacketSize, r.block_count,
              r.packets.size() * 8.0 / (double(samples) / rate) / 1000);
  return 0;
}
