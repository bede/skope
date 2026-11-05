#!/usr/bin/env bash
# Deploy the demo to the Jekyll blog (bede.github.io) as a standalone static app at /skope/,
# once the compiled wasm matches the CLI. `--no-push` commits without pushing, and PAGES_ROOT
# overrides the blog checkout beside this repo
set -euo pipefail

SRC="$(cd "$(dirname "$0")" && pwd)"
REPO="$(dirname "$SRC")"
PAGES_ROOT="${PAGES_ROOT:-$(dirname "$REPO")/bede.github.io}"
DEST="$PAGES_ROOT/skope"

# The deployed app. skope/ holds only these, so a rename can't leave an orphan behind
FILES=(index.html fonts.css app.js pairing.js protocol.js worker.js pkg/skope_wasm.js pkg/skope_wasm_bg.wasm
  fonts/{fraunces,geist,geist-mono}-latin{,-ext}.woff2)

git -C "$PAGES_ROOT" rev-parse --git-dir > /dev/null

SOURCE_HASH="$(git -C "$REPO" rev-parse --short HEAD)"
# We deploy the working tree, so label it dirty when what it's built from is uncommitted
SOURCES=(src Cargo.toml Cargo.lock skope-wasm/src skope-wasm/Cargo.toml "${FILES[@]/#/skope-wasm/}")
if [[ -n "$(git -C "$REPO" status --porcelain -- "${SOURCES[@]}")" ]]; then
  SOURCE_HASH="$SOURCE_HASH-dirty-$(date -u +%Y%m%d%H%M%S)"
fi

# Deploy only wasm that matches the CLI, then build the browser package from the same source
cargo build --release --manifest-path "$REPO/Cargo.toml"
"$SRC/build.sh" nodejs
node "$SRC/parity.mjs"
"$SRC/build.sh"

rm -rf "$DEST"
mkdir -p "$DEST"
(cd "$SRC" && rsync -R "${FILES[@]}" "$DEST/")

# One deployment is one coherent module graph, so stamp every cache token with the revision
sed -E -i.bak "s/[?]v=dev\"/?v=$SOURCE_HASH\"/g" "$DEST"/*.html "$DEST"/*.css "$DEST"/*.js
rm "$DEST"/*.bak
echo "Deployed skope-wasm $SOURCE_HASH to $DEST"

cd "$PAGES_ROOT"
git add -A -- skope
if git diff --cached --quiet -- skope; then
  echo "No GitHub Pages changes to commit"
  exit 0
fi
# Commit skope/ alone, leaving anything else staged in the blog for its owner
git commit -q -m "skope-wasm $SOURCE_HASH" -- skope
echo "Committed skope-wasm $SOURCE_HASH to $PAGES_ROOT"

if [[ "${1:-}" != "--no-push" ]]; then
  git push origin "$(git rev-parse --abbrev-ref HEAD)"
  echo "Pushed to GitHub Pages: https://$(cat CNAME 2>/dev/null || echo bede.github.io)/skope/"
fi
