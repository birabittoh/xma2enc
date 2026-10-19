# xma2enc

A small XMA2 encoder: 16 bit PCM in, XMA2 packets out. The bitstream is WMA Pro
with the XMA2 packet layer, one 512 sample block per frame. Details are in
[docs/format.md](docs/format.md).

## Library

```cpp
#include "xma2enc/xma2enc.h"

xma2enc::Options options;          // bitrate, block_size
xma2enc::Result result;
bool ok = xma2enc::Encode(pcm, samples, channels, sample_rate, options, &result);
// result.packets, result.seek_table, result.frames, ...
```

`pcm` is interleaved, 1 or 2 channels, up to 48000 Hz. Add it to a CMake
project with `add_subdirectory(xma2enc)` and link the `xma2enc` target.

## Command line

```bash
cmake -B build && cmake --build build
build/xma2enc_cli in.wav out.xma2.wav [bitrate] [block_size]
```

The CLI target is only built when this directory is the top level project.

## License

LGPL 2.1 or later, see [COPYING.LGPLv2.1](COPYING.LGPLv2.1).
