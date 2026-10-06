#!/usr/bin/env bash
# Deploy the demo to the Jekyll blog (bede.github.io) as a standalone static app at /skope/,
# once the compiled wasm matches the CLI. `--no-push` commits without pushing, and PAGES_ROOT
# overrides the blog checkout beside this repo
set -euo pipefail

SRC="$(cd "$(dirname "$0")" && pwd)"
REPO="$(dirname "$SRC")"
PAGES_ROOT="${PAGES_ROOT:-$(dirname "$REPO")/bede.github.io}"
DEST="$PAGES_ROOT/skope"
PUSH=1
[[ "${1:-}" == "--no-push" ]] && PUSH=0

# Files deployed into skope/ (paths relative to both SRC and DEST). This is the single source
# of truth: files copied, stamped, staged and pruned all derive from it, so a rename in SRC
# can't leave an orphaned copy behind in the deployed app
FILES=(
  index.html
  app.js
  pairing.js
  protocol.js
  worker.js
  pkg/skope_wasm.js
  pkg/skope_wasm_bg.wasm
)

# Everything the deployed bytes are built from, relative to REPO
SOURCES=(src Cargo.toml Cargo.lock skope-wasm/src skope-wasm/Cargo.toml "${FILES[@]/#/skope-wasm/}")

git -C "$PAGES_ROOT" rev-parse --git-dir > /dev/null

SOURCE_HASH="$(git -C "$REPO" rev-parse --short HEAD)"
# We deploy the working tree, not HEAD, so don't label uncommitted code with a clean hash.
# Repeated dirty deployments may contain different bytes at the same HEAD
if [[ -n "$(git -C "$REPO" status --porcelain -- "${SOURCES[@]}")" ]]; then
  SOURCE_HASH="$SOURCE_HASH-dirty-$(date -u +%Y%m%d%H%M%S)"
fi

# Deploy only wasm that matches the CLI, then build the browser package from the same source
cargo build --release --manifest-path "$REPO/Cargo.toml"
"$SRC/build.sh" nodejs
node "$SRC/parity.mjs"
"$SRC/build.sh"

mkdir -p "$DEST/pkg"
for f in "${FILES[@]}"; do
  cp "$SRC/$f" "$DEST/$f"
done

# One deployment is one coherent module graph. Stamp the build label and every cache token
# with the revision, so a new entry point can never import an older dependency
for f in index.html app.js worker.js; do
  sed -E -i.bak "s/(BUILD = \"|[?]v=)dev\"/\\1$SOURCE_HASH\"/g" "$DEST/$f"
  rm "$DEST/$f.bak"
done
if grep -nE '(BUILD = "|[?]v=)dev"' "$DEST"/*.html "$DEST"/*.js; then
  echo "Unstamped cache tokens remain in $DEST" >&2
  exit 1
fi
echo "Deployed skope-wasm $SOURCE_HASH to $DEST"

cd "$PAGES_ROOT"

# Prune any tracked file under skope/ that is no longer in FILES (e.g. left behind by a
# rename), then stage the current set
while IFS= read -r tracked; do
  rel="${tracked#skope/}"
  keep=0
  for f in "${FILES[@]}"; do [[ "$f" == "$rel" ]] && { keep=1; break; }; done
  [[ "$keep" == 0 ]] && git rm -q "$tracked"
done < <(git ls-files skope/)
git add "${FILES[@]/#/skope/}"

if git diff --cached --quiet -- skope; then
  echo "No GitHub Pages changes to commit"
  exit 0
fi
# Commit skope/ alone, leaving anything else staged in the blog for its owner
git commit -q -m "skope-wasm $SOURCE_HASH" -- skope
echo "Committed skope-wasm $SOURCE_HASH to $PAGES_ROOT"

if (( PUSH )); then
  git push origin "$(git rev-parse --abbrev-ref HEAD)"
  echo "Pushed to GitHub Pages: https://$(cat CNAME 2>/dev/null || echo bede.github.io)/skope/"
fi
