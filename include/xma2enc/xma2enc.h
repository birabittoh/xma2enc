// General purpose XMA2 encoder: PCM in, XMA2 packets out.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace xma2enc {

struct Options {
  // Target average bitrate for the whole stream, bits per second.
  uint32_t bitrate = 128000;
  // Bytes per block. Frames never cross a block boundary, so every block can
  // be decoded on its own. Must be a multiple of 2048.
  uint32_t block_size = 0x10000;
};

struct Result {
  std::vector<uint8_t> packets;       // whole 2048 byte packets
  std::vector<uint32_t> seek_table;   // cumulative samples at the end of each block
  uint32_t channels = 0;
  uint32_t sample_rate = 0;
  uint32_t frames = 0;                // XMA frames, 512 samples each
  uint32_t samples = 0;               // input samples per channel
  uint32_t block_size = 0;
  uint32_t block_count = 0;
};

constexpr uint32_t kPacketSize = 2048;
constexpr uint32_t kFrameSamples = 512;
// Output samples before the first input sample when decoding from the start
// (one skipped frame plus the 64 sample decoder delay).
constexpr uint32_t kDecoderDelay = 576;

// pcm is interleaved, 1 or 2 channels, sample_rate up to 48000. The rate only
// picks the scale factor band layout, as in the decoder (above 44100 is
// treated as 48000, above 32000 as 44100, and so on). Returns false on bad
// parameters.
bool Encode(const int16_t* pcm, size_t samples, uint32_t channels,
            uint32_t sample_rate, const Options& options, Result* out);

}  // namespace xma2enc
