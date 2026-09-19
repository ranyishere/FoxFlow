#!/usr/bin/env bash
#
# run_example.sh — thin wrapper around the FoxFlow compiler driver (ir.jl).
#
# Usage:
#   ./run_example.sh [command] [example_name] [options] [-- exec args...]
#
# Arguments:
#   command       gen | build | run | gdb   (default: gen)
#   example_name  Folder under src/tests/ (default: microtubules)
#   exec args     Anything after "--" is forwarded to the built model (run only)
#
# Extra options are forwarded straight to ir.jl, e.g. --sundials, --project, --target.
#
# Examples:
#   ./run_example.sh                          # generate C++ for "microtubules"
#   ./run_example.sh gen fracture_network     # generate for another example
#   ./run_example.sh build microtubules       # generate + compile
#   ./run_example.sh run microtubules --project "$PWD/src/tests/generated_tests"
#   ./run_example.sh gdb microtubules         # build (Debug) + launch under gdb
#
#   # In a second terminal, keep a ParaView .pvd in sync with the running sim:
#   ./run_example.sh watch <build_dir>/my_results
#
# Environment:
#   FOXFLOW_SUNDIALS_DIR   passed through to CMake as -DSUNDIALS_DIR
#
set -euo pipefail

# Resolve the repository root (directory containing this script).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IR_SCRIPT="$SCRIPT_DIR/src/ir_generation/ir.jl"

# First positional may be a command; otherwise default to "gen".
COMMAND="gen"
# Extra flags injected by convenience commands (e.g. "gdb" -> run --gdb).
EXTRA_FLAGS=()
case "${1:-}" in
    gen|build|run)
        COMMAND="$1"; shift ;;
    gdb)
        # Convenience: build a Debug model and drop into gdb.
        COMMAND="run"; EXTRA_FLAGS+=("--gdb"); shift ;;
    watch)
        # ./run_example.sh watch <results_dir> [--once ...]
        shift
        julia "$IR_SCRIPT" watch "$@"
        exit $? ;;
    help|-h|--help)
        julia "$IR_SCRIPT" help; exit 0 ;;
esac

EXAMPLE="${1:-microtubules}"
if [[ $# -gt 0 ]]; then shift; fi   # drop example name; keep remaining args/options

INPUT_DIR="$SCRIPT_DIR/src/tests/$EXAMPLE"
OUTPUT_DIR="$SCRIPT_DIR/src/tests/generated_tests/$EXAMPLE"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: example '$EXAMPLE' not found at: $INPUT_DIR"
    echo
    echo "Available examples:"
    for d in "$SCRIPT_DIR"/src/tests/*/; do
        # Only list folders that look like FoxFlow examples.
        [[ -f "$d/simulation.fflow" ]] && echo "  - $(basename "$d")"
    done
    exit 1
fi

echo "Command:    $COMMAND"
echo "Example:    $EXAMPLE"
echo "Input dir:  $INPUT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo

# Remaining "$@" are extra options/flags (e.g. --sundials, --project, -- exec args).

# Print out the julia command for debugging purposes.

echo "Running:"
echo "  julia \"$IR_SCRIPT\" \"$COMMAND\" \"$INPUT_DIR\" \"$OUTPUT_DIR\" ${EXTRA_FLAGS[*]:-} \"$@\""
echo
julia "$IR_SCRIPT" "$COMMAND" "$INPUT_DIR" "$OUTPUT_DIR" ${EXTRA_FLAGS[@]+"${EXTRA_FLAGS[@]}"} "$@"
