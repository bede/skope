#!/bin/sh
# Serve the demo at http://localhost:8000, as module workers cannot load from file://
set -eu
cd "$(dirname "$0")"
python3 -m http.server "${1:-8000}"
