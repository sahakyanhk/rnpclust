#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN_DIR="$SCRIPT_DIR/bin"

echo "Checking prerequisites..."

ok=true
for cmd in USalign parallel python3; do
  if command -v "$cmd" &>/dev/null; then
    echo "  $cmd: $(command -v "$cmd")"
  else
    echo "  $cmd: NOT FOUND" >&2
    ok=false
  fi
done

if [ "$ok" = false ]; then
  echo "Missing prerequisites. Install them and re-run." >&2
  exit 1
fi

chmod +x "$BIN_DIR"/*

echo ""
echo "Add this to your ~/.bashrc or ~/.bash_profile:"
echo ""
echo "  export PATH=\"$BIN_DIR:\$PATH\""
echo ""
echo "Then run: source ~/.bashrc"
echo ""
echo "Usage: uniclust <pdb_dir> [output_dir] [cutoff]"
