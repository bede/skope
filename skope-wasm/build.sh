#!/bin/sh
# Build for browsers into pkg/, or for Node into pkg-node/ with `./build.sh nodejs`
set -eu
cd "$(dirname "$0")"
target="${1:-web}"
out=pkg
[ "$target" = nodejs ] && out=pkg-node
wasm-pack build --release --target "$target" --out-dir "$out" --no-pack
