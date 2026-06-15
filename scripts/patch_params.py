#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def replace_define(text: str, name: str, replacement: str) -> str:
    pat = re.compile(rf"^\s*#define\s+{re.escape(name)}\b.*$", re.MULTILINE)
    new_line = f"#define {name} {replacement}"
    if not pat.search(text):
        raise RuntimeError(f"Could not find #define {name}")
    return pat.sub(new_line, text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default="src/params.h")
    ap.add_argument("--nx", type=int, required=True, help="Total grid points, e.g. 401 or 801")
    args = ap.parse_args()

    path = Path(args.params)
    text = path.read_text()

    if args.nx < 3:
        raise ValueError("Nx must be at least 3")

    # Preserve the original style Nx=(factor*N_BS+1) for the single-droplet code.
    factor = args.nx - 1
    text = replace_define(text, "Nx", f"({factor} * N_BS + 1)")

    path.write_text(text)
    print(f"patched {path}: Nx -> ({factor} * N_BS + 1) = {args.nx} for N_BS=1")


if __name__ == "__main__":
    main()
