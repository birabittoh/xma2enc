# XMA2 bitstream notes

What the encoder writes, as established against FFmpeg's `wmaprodec.c`.

## Stream

* WMA Pro derived. Decode flags 0x10d6, 512 sample frames, one stream
  (stereo or mono), 15 bit frame length, a DRC byte per frame, one pad bit
  and a trailer bit.
* The encoder codes one long 512 sample block per frame. Stereo bands are
  coded mid/side, coefficients as vector plus run/level symbols, with power
  law noise shaping (`Options::gamma`).
* Frame count is `ceil((n + 576) / 512)`: the first frame is skipped by the
  decoder and FFmpeg drops 64 more samples.

## Packets

2048 bytes. Header, big endian u32: frames starting in the packet (6 bits),
bit offset of the first of them after the header (15), metadata (3, always 1),
packet skip count (8, always 0).

* Frames never cross a block boundary; the rest of the block is zero padded.
* The trailer bit must be 0 on the last frame that starts in a packet,
  otherwise FFmpeg parses the next packet's header as a frame.

## Seek table

One u32 per block: 512 times the frames starting up to the block's end.

## Container

`tools/xma2enc_cli.cpp` wraps the packets in a RIFF 0x0166 XMA2WAVEFORMATEX
(34 byte extradata: streams, channel mask, samples encoded, bytes per block,
play begin and length, loop begin and length, loop count, encoder version 4,
block count), which FFmpeg decodes as is.
