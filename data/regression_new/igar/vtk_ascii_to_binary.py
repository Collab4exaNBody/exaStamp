#!/usr/bin/python3
"""
Converts an ASCII legacy VTK file (DATASET STRUCTURED_POINTS, one scalar per point) into a BINARY one, as read by
igar_read_field. The values are streamed in chunks: a 1024^3 field never has to fit in memory.

Usage:
    vtk_ascii_to_binary.py input.vtk output.vtk [--dtype auto|uint8|float32|float64] [--scale S]

--dtype
    auto (default)  uint8 if the file declares an integer type (unsigned_char, int, vtktypeint64, ...) or COLOR_SCALARS
                    (VTK's 8-bit arrays, written as byte/255 in ASCII: converted back to bytes), float64 otherwise
    uint8           1 byte per value (8-bit images). Every value times --scale must be an integer 0-255, otherwise the
                    conversion stops: values are never silently rounded.
    float32/float64 4/8 bytes per value
--scale S
    uint8 only: stored byte = value * S (e.g. 255 for an 8-bit image normalized to [0, 1]). Read it back with
    igar_read_field { storage: uint8, value_scale: 1/S }.

The point order (F or C) is kept as is: use the same grid_ordering in igar_read_field as for the ASCII file.
"""

import argparse
import sys
import numpy as np

INTEGER_TYPES = {"unsigned_char", "char", "unsigned_short", "short", "unsigned_int", "int", "unsigned_long", "long",
                 "vtktypeint8", "vtktypeuint8", "vtktypeint16", "vtktypeuint16", "vtktypeint32", "vtktypeuint32",
                 "vtktypeint64", "vtktypeuint64"}
# legacy VTK binary data is big-endian
OUT_TYPES = {"uint8": (">u1", "unsigned_char"), "float32": (">f4", "float"), "float64": (">f8", "double")}
CHUNK_BYTES = 1 << 26


def read_header(f):
    """Returns the header lines (bytes) up to the data, the point count, the declared data type and whether the data are
    COLOR_SCALARS (VTK's 8-bit arrays: byte/255 in ASCII files, no LOOKUP_TABLE)."""
    lines = []
    count = None
    dtype = None
    color = False
    while True:
        line = f.readline()
        if not line:
            sys.exit("error: LOOKUP_TABLE not found")
        lines.append(line)
        words = line.split()
        if not words:
            continue
        key = words[0].decode()
        if key == "BINARY":
            sys.exit("error: the input file is already BINARY")
        if key == "DATASET" and words[1].decode() != "STRUCTURED_POINTS":
            sys.exit("error: DATASET %s not supported (STRUCTURED_POINTS only)" % words[1].decode())
        if key == "POINT_DATA":
            count = int(words[1])
        if key == "SCALARS":
            dtype = words[2].decode()
            if len(words) > 3 and int(words[3]) != 1:
                sys.exit("error: scalar field expected, %s components found" % words[3].decode())
        if key == "LOOKUP_TABLE":
            break
        if key == "COLOR_SCALARS":
            if int(words[2]) != 1:
                sys.exit("error: scalar field expected, COLOR_SCALARS with %s values found" % words[2].decode())
            # rewritten as SCALARS unsigned_char + LOOKUP_TABLE
            lines[-1] = b"SCALARS %s unsigned_char\n" % words[1]
            lines.append(b"LOOKUP_TABLE default\n")
            dtype = "unsigned_char"
            color = True
            break
    if count is None or dtype is None:
        sys.exit("error: POINT_DATA or SCALARS line missing")
    return lines, count, dtype, color


def ascii_values(f):
    """Yields float64 arrays of the values, chunk by chunk (a chunk never cuts a number)."""
    tail = b""
    while True:
        block = f.read(CHUNK_BYTES)
        if not block:
            break
        block = tail + block
        cut = max(block.rfind(b" "), block.rfind(b"\n"), block.rfind(b"\t"), block.rfind(b"\r"))
        if cut < 0:
            tail = block
            continue
        tail = block[cut + 1:]
        text = block[:cut + 1].split()
        if text:
            yield np.array(text, dtype=np.float64)
    text = tail.split()
    if text:
        yield np.array(text, dtype=np.float64)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--dtype", default="auto", choices=["auto", "uint8", "float32", "float64"])
    parser.add_argument("--scale", type=float, default=1.0)
    args = parser.parse_args()

    with open(args.input, "rb") as fin:
        header, count, in_type, color = read_header(fin)
        out = args.dtype
        if out == "auto":
            out = "uint8" if in_type in INTEGER_TYPES else "float64"
        if args.scale != 1.0 and out != "uint8":
            sys.exit("error: --scale only applies to --dtype uint8")
        np_type, vtk_type = OUT_TYPES[out]

        with open(args.output, "wb") as fout:
            for line in header:
                words = line.split()
                if words and words[0] == b"ASCII":
                    line = b"BINARY\n"
                elif words and words[0] == b"SCALARS":
                    line = b"SCALARS %s %s\n" % (words[1], vtk_type.encode())
                fout.write(line)

            n = 0
            vmin, vmax = np.inf, -np.inf
            for values in ascii_values(fin):
                if n + values.size > count:
                    sys.exit("error: more than POINT_DATA=%d values" % count)
                if color:
                    values = np.rint(values * 255.0)
                if out == "uint8":
                    scaled = values * args.scale
                    levels = np.rint(scaled)
                    bad = (levels < 0) | (levels > 255) | (np.abs(levels - scaled) > 1.0e-6 * max(1.0, args.scale))
                    if bad.any():
                        i = int(np.argmax(bad))
                        sys.exit("error: value %.17g at index %d is not an integer 0-255 after --scale %g (use --dtype "
                                 "float64, or the right --scale)" % (values[i], n + i, args.scale))
                    values = levels
                vmin, vmax = min(vmin, values.min()), max(vmax, values.max())
                fout.write(values.astype(np_type).tobytes())
                n += values.size
            fout.write(b"\n")
            if n != count:
                sys.exit("error: %d values read, POINT_DATA=%d" % (n, count))

    print("%s -> %s : %d values, %s (%s), stored range [%g, %g]" % (args.input, args.output, count, out, vtk_type, vmin, vmax))
    if out == "uint8" and args.scale != 1.0:
        print("read with igar_read_field { storage: uint8, value_scale: %.17g }" % (1.0 / args.scale))


if __name__ == "__main__":
    main()
