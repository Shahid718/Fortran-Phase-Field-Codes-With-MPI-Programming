#!/usr/bin/env bash
#-------------------------------------------------------------------------------
# run_grid_sweep.sh
#
# Compiles and runs main.f90 for a set of grid sizes (Nx = Ny = N),
# saving each run's output data and timing log separately.
#
# Usage:
#   ./run_grid_sweep.sh                     # uses default sizes below
#   ./run_grid_sweep.sh 64 128 256 512 1024  # or pass sizes explicitly
#
# Requires: main.f90 in the same directory, mpif90/mpirun on PATH.
#-------------------------------------------------------------------------------
set -euo pipefail

SRC="main.f90"
NP=1                    # number of MPI ranks
RESULTS_DIR="grid_sweep_results"

# Grid sizes to test. Edit this list or pass sizes as command-line args.
SIZES=(64 128 256 512 1024)
if [ "$#" -gt 0 ]; then
    SIZES=("$@")
fi

if [ ! -f "$SRC" ]; then
    echo "Error: $SRC not found in current directory." >&2
    exit 1
fi

mkdir -p "$RESULTS_DIR"

for N in "${SIZES[@]}"; do
    echo "==================================================================="
    echo " Grid size: ${N} x ${N}"
    echo "==================================================================="

    WORKSRC="main_N${N}.f90"
    EXE="main_N${N}"

    # Substitute Nx and Ny parameter values
    sed -e "s/integer, parameter :: Nx = [0-9]\+/integer, parameter :: Nx = ${N}/" \
        -e "s/integer, parameter :: Ny = [0-9]\+/integer, parameter :: Ny = ${N}/" \
        "$SRC" > "$WORKSRC"

    # Compile
    if ! mpif90 -std=gnu "$WORKSRC" -O2 -o "$EXE" ; then
        echo "Compilation failed for N=${N}, skipping." >&2
        rm -f "$WORKSRC"
        continue
    fi

    # Run and capture output + wall time
    LOGFILE="${RESULTS_DIR}/log_N${N}.txt"
    START=$(date +%s.%N)
    mpirun -np "$NP" "./${EXE}" | tee "$LOGFILE"
    END=$(date +%s.%N)
    WALL=$(echo "$END - $START" | bc)
    echo "Wall clock (script-measured) = ${WALL} s" | tee -a "$LOGFILE"

    # Save the output data under a size-specific name
    if [ -f ch.dat ]; then
        mv ch.dat "${RESULTS_DIR}/ch_N${N}.dat"
    else
        echo "Warning: ch.dat not produced for N=${N}" | tee -a "$LOGFILE"
    fi

    # Clean up the per-size source and executable
    rm -f "$WORKSRC" "$EXE"

    echo
done

echo "All runs complete. Results in ${RESULTS_DIR}/"