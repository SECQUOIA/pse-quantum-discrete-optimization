#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NOTEBOOK="${1:-"$ROOT/ds_mfg_qubo_qiskitopt.ipynb"}"
KERNEL_NAME="${KERNEL_NAME:-qiskitopt-julia}"
KERNEL_ROOT="${KERNEL_ROOT:-"$ROOT/.jupyter"}"
IJULIA_ENV="${IJULIA_ENV:-"$KERNEL_ROOT/ijulia-env"}"
KERNEL_DIR="$KERNEL_ROOT/kernels/$KERNEL_NAME"
TIMEOUT="${EXECUTE_TIMEOUT:-7200}"

if [[ -z "${JULIA_BIN:-}" ]]; then
    if command -v julia >/dev/null 2>&1; then
        JULIA_BIN="$(command -v julia)"
    else
        echo "Could not find Julia. Set JULIA_BIN=/path/to/julia. Julia 1.10 or newer is recommended." >&2
        exit 1
    fi
fi

if [[ -z "${JUPYTER_BIN:-}" ]]; then
    if command -v jupyter >/dev/null 2>&1; then
        JUPYTER_BIN="$(command -v jupyter)"
    else
        echo "Could not find jupyter. Set JUPYTER_BIN=/path/to/jupyter." >&2
        exit 1
    fi
fi

JULIA_DEPOT_PATH_VALUE="${JULIA_DEPOT_PATH_VALUE:-"$ROOT/.julia-depot-nbconvert:${HOME:-$ROOT/.julia-home}/.julia"}"

mkdir -p "$KERNEL_DIR" "$IJULIA_ENV"

cat > "$IJULIA_ENV/Project.toml" <<'EOF'
[deps]
IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
EOF

JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH_VALUE" \
JULIA_PKG_PRECOMPILE_AUTO=0 \
"$JULIA_BIN" --compiled-modules=no --project="$IJULIA_ENV" -e 'import Pkg; Pkg.instantiate(); import IJulia; println("IJulia nbconvert kernel ready")'

cat > "$KERNEL_DIR/kernel.json" <<EOF
{
  "display_name": "QiskitOpt Julia",
  "language": "julia",
  "interrupt_mode": "signal",
  "argv": [
    "$JULIA_BIN",
    "-i",
    "--color=yes",
    "--compiled-modules=no",
    "--project=$IJULIA_ENV",
    "-e",
    "import IJulia; IJulia.run_kernel()",
    "{connection_file}"
  ],
  "env": {
    "JULIA_DEPOT_PATH": "$JULIA_DEPOT_PATH_VALUE",
    "JULIA_PKG_PRECOMPILE_AUTO": "0"
  }
}
EOF

JUPYTER_PATH="$KERNEL_ROOT${JUPYTER_PATH:+:$JUPYTER_PATH}" \
"$JUPYTER_BIN" nbconvert \
    --to notebook \
    --execute "$NOTEBOOK" \
    --output "$(basename "${NOTEBOOK%.ipynb}")_executed.ipynb" \
    --output-dir "$ROOT" \
    --ExecutePreprocessor.kernel_name="$KERNEL_NAME" \
    --ExecutePreprocessor.timeout="$TIMEOUT"
