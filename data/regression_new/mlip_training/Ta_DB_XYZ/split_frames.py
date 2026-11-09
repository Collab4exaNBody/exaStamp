#!/usr/bin/env python3
"""Split each multi-frame .xyz file in this dir into one file per frame,
under a folder named after the input file's prefix (stem)."""
import glob
import os


def split_xyz(path):
    prefix = os.path.splitext(os.path.basename(path))[0]
    os.makedirs(prefix, exist_ok=True)

    with open(path) as f:
        lines = f.readlines()

    i = 0
    frame = 0
    while i < len(lines):
        natoms = int(lines[i].strip())
        frame_lines = lines[i:i + 2 + natoms]
        out_path = os.path.join(prefix, f"{prefix}_{frame:04d}.xyz")
        with open(out_path, "w") as out:
            out.writelines(frame_lines)
        i += 2 + natoms
        frame += 1

    return frame


if __name__ == "__main__":
    for path in sorted(glob.glob("*.xyz")):
        n = split_xyz(path)
        print(f"{path}: {n} frames")
