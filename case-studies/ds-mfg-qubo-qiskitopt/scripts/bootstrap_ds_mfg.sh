#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -z "${JULIA_BIN:-}" ]]; then
    if command -v julia >/dev/null 2>&1; then
        JULIA_BIN="$(command -v julia)"
    else
        echo "Could not find Julia. Install Julia 1.10 or newer, or set JULIA_BIN=/path/to/julia." >&2
        exit 1
    fi
fi

"$JULIA_BIN" --project="$ROOT" -e '
import Pkg
VERSION >= v"1.10" || error("Julia 1.10 or newer is required; found $(VERSION)")
Pkg.instantiate()
Pkg.precompile()
'

"$JULIA_BIN" --project="$ROOT" "$ROOT/scripts/smoke_test.jl"
