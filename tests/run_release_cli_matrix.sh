#!/usr/bin/env bash
set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MBS="${MBS:-$ROOT_DIR/cpu/bin/mbs_po_mpi}"
TIMEOUT_SECONDS="${MBS_RELEASE_TIMEOUT_SECONDS:-30}"
KEEP_WORK="${MBS_RELEASE_KEEP_WORK:-0}"
MPI_RUNNER="${MBS_RELEASE_MPIRUN:-}"

if [[ "$MBS" != /* ]]; then
    MBS="$(cd "$(dirname "$MBS")" 2>/dev/null && pwd)/$(basename "$MBS")"
fi

if [[ ! -x "$MBS" ]]; then
    printf 'ERROR: release binary is not executable: %s\n' "$MBS" >&2
    printf '  Fix: build cpu/bin/mbs_po_mpi or set MBS to that executable.\n' >&2
    exit 2
fi
if ! command -v timeout >/dev/null 2>&1; then
    printf 'ERROR: GNU timeout is required by the release CLI matrix.\n' >&2
    printf '  Fix: install coreutils and rerun this script.\n' >&2
    exit 2
fi

if [[ -z "$MPI_RUNNER" ]]; then
    for candidate in /usr/bin/mpirun "$(command -v mpirun 2>/dev/null || true)"; do
        [[ -n "$candidate" && -x "$candidate" ]] || continue
        if "$candidate" --version 2>/dev/null | grep -qi 'Open MPI'; then
            MPI_RUNNER="$candidate"
            break
        fi
    done
fi
if [[ -n "$MPI_RUNNER" ]]; then
    if [[ ! -x "$MPI_RUNNER" ]] \
        || ! "$MPI_RUNNER" --version 2>/dev/null | grep -qi 'Open MPI'; then
        printf 'ERROR: MBS_RELEASE_MPIRUN is not an executable Open MPI launcher: %s\n' \
            "$MPI_RUNNER" >&2
        printf '  Fix: point MBS_RELEASE_MPIRUN to mpirun from the same Open MPI installation used by mpicxx.\n' >&2
        exit 2
    fi
fi
if [[ ! "$TIMEOUT_SECONDS" =~ ^[1-9][0-9]*([.][0-9]+)?$ ]]; then
    printf 'ERROR: MBS_RELEASE_TIMEOUT_SECONDS must be a positive number; got %q.\n' \
        "$TIMEOUT_SECONDS" >&2
    printf '  Fix: use a value such as 30.\n' >&2
    exit 2
fi

HELP_DEBUG_OUTPUT="$("$MBS" --help-debug 2>&1)"
if (($? != 0)); then
    printf 'ERROR: the release binary failed to generate --help-debug output.\n' >&2
    printf '  Fix: repair help generation before running the release matrix.\n' >&2
    exit 2
fi
if ! "$MBS" --help >/dev/null 2>&1; then
    printf 'ERROR: the release binary failed to generate --help output.\n' >&2
    printf '  Fix: repair production help generation before running the release matrix.\n' >&2
    exit 2
fi
VERSION_OUTPUT="$("$MBS" --version 2>&1)"
if (($? != 0)) || ! grep -q '^MBS-fast ' <<<"$VERSION_OUTPUT" \
    || ! grep -q '^backend-build: ' <<<"$VERSION_OUTPUT" \
    || ! grep -q '^compiler: ' <<<"$VERSION_OUTPUT" \
    || ! grep -q '^openmp: ' <<<"$VERSION_OUTPUT" \
    || ! grep -q '^mpi: ' <<<"$VERSION_OUTPUT"; then
    printf 'ERROR: --version output is incomplete or failed.\n' >&2
    printf '  Fix: report revision, backend, compiler, OpenMP, and MPI state.\n' >&2
    exit 2
fi
missing_coverage=()
while IFS= read -r option; do
    if ! grep -Eq -- "(^|[[:space:]\"'\\(])${option}([[:space:]\"'\\)]|$)" "$0"; then
        missing_coverage+=("$option")
    fi
done < <(printf '%s\n' "$HELP_DEBUG_OUTPUT" \
    | sed -n 's/^  \(--[a-z0-9-]*\).*/\1/p')
if ((${#missing_coverage[@]} != 0)); then
    printf 'ERROR: release matrix does not mention registered option(s):' >&2
    printf ' %s' "${missing_coverage[@]}" >&2
    printf '\n  Fix: add a success or expected-error case for every listed option.\n' >&2
    exit 2
fi

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/mbs-release-cli.XXXXXX")" || exit 2
RUNTIME_DIR="$WORK_DIR/runtime"
LOG_DIR="$WORK_DIR/logs"
RESULT_DIR="$WORK_DIR/results"
mkdir -p "$RUNTIME_DIR" "$LOG_DIR" "$RESULT_DIR" || exit 2

cleanup()
{
    local status=$?
    if [[ "$KEEP_WORK" == "1" ]]; then
        printf 'Release CLI matrix artifacts: %s\n' "$WORK_DIR"
    else
        rm -rf -- "$WORK_DIR"
    fi
    return "$status"
}
trap cleanup EXIT

ORIENTATION_FILE="$WORK_DIR/orientations.dat"
THETA_FILE="$WORK_DIR/theta.dat"
K_EQ_FILE="$WORK_DIR/k_eq.dat"
ADAPTIVE_CONFIG="$WORK_DIR/adaptive.conf"
PARTICLE_FILE="$WORK_DIR/cube.dat"
EXAMPLE_PARTICLE="$ROOT_DIR/examples/cube.particle"
EMPTY_PARTICLE_FILE="$WORK_DIR/empty-particle.dat"
BAD_PARTICLE_FILE="$WORK_DIR/bad-particle.dat"
ZERO_VOLUME_PARTICLE_FILE="$WORK_DIR/zero-volume-particle.dat"
DUPLICATE_THETA_FILE="$WORK_DIR/theta-duplicate.dat"
TRAJECTORY_FILE="$WORK_DIR/trajectories.dat"
BAD_TRAJECTORY_FILE="$WORK_DIR/trajectories-bad.dat"
MISSING_TRAJECTORY_GROUP_FILE="$WORK_DIR/trajectories-missing-group.dat"
OVERSIZED_TRAJECTORY_GROUP_FILE="$WORK_DIR/trajectories-oversized-group.dat"

printf '0 0\n45 30\n' >"$ORIENTATION_FILE" || exit 2
printf '0\n90\n180\n' >"$THETA_FILE" || exit 2
printf '0.10\n0.20\n' >"$K_EQ_FILE" || exit 2
printf '%s\n' \
    'reflections.min = 1' \
    'reflections.max = 2' \
    'reflections.tolerance = 0.9' \
    'alpha.min_points = 12' \
    'alpha.max_points = 24' \
    'alpha.tolerance = 0.9' \
    'theta.min_points = 17' \
    'theta.max_points = 33' \
    'theta.tolerance = 0.9' \
    'orientations.min = 16' \
    'orientations.max = 64' \
    'orientations.tolerance = 0.9' \
    'beta.min_points = 2' \
    'beta.max_points = 4' \
    'beta.tolerance = 0.9' \
    'gamma.min_points = 6' \
    'gamma.max_points = 12' \
    'gamma.tolerance = 0.9' \
    'controller.stable_passes = 1' \
    'controller.min_pilot_orientations = 16' \
    'controller.max_pilot_orientations = 16' \
    'controller.max_joint_sweeps = 1' \
    'controller.max_euler_sweeps = 1' \
    >"$ADAPTIVE_CONFIG" || exit 2
printf '%s\n' \
    '0' \
    '0' \
    '90 90' \
    '' \
    '-0.5 -0.5 -0.5' '0.5 -0.5 -0.5' '0.5 0.5 -0.5' '-0.5 0.5 -0.5' \
    '' \
    '-0.5 -0.5 0.5' '-0.5 0.5 0.5' '0.5 0.5 0.5' '0.5 -0.5 0.5' \
    '' \
    '-0.5 -0.5 -0.5' '-0.5 -0.5 0.5' '0.5 -0.5 0.5' '0.5 -0.5 -0.5' \
    '' \
    '0.5 -0.5 -0.5' '0.5 -0.5 0.5' '0.5 0.5 0.5' '0.5 0.5 -0.5' \
    '' \
    '0.5 0.5 -0.5' '0.5 0.5 0.5' '-0.5 0.5 0.5' '-0.5 0.5 -0.5' \
    '' \
    '-0.5 0.5 -0.5' '-0.5 0.5 0.5' '-0.5 -0.5 0.5' '-0.5 -0.5 -0.5' \
    >"$PARTICLE_FILE" || exit 2
printf '' >"$EMPTY_PARTICLE_FILE" || exit 2
printf '%s\n' '0' '0' '90 90' '' '0 0' >"$BAD_PARTICLE_FILE" || exit 2
printf '%s\n' \
    '0' '0' '90 90' '' \
    '0 0 0' '1 0 0' '0 1 0' \
    >"$ZERO_VOLUME_PARTICLE_FILE" || exit 2
printf '%s\n' '0' '90' '90' '180' >"$DUPLICATE_THETA_FILE" || exit 2
printf '%s\n' '0' '1 2 : 0' >"$TRAJECTORY_FILE" || exit 2
printf '%s\n' '0 : 1024' >"$BAD_TRAJECTORY_FILE" || exit 2
printf '%s\n' '0 :' >"$MISSING_TRAJECTORY_GROUP_FILE" || exit 2
for ((i = 0; i <= 1024; ++i)); do
    printf '%s\n' '0 : 0'
done >"$OVERSIZED_TRAJECTORY_GROUP_FILE" || exit 2

TOTAL=0
PASSED=0
FAILED=0
RUN_TOTAL=0
RUN_PASSED=0
ERROR_TOTAL=0
ERROR_PASSED=0

CRASH_PATTERN='AddressSanitizer|UndefinedBehaviorSanitizer|Segmentation fault|segfault|^Aborted( \(core dumped\))?$|core dumped|Process received signal|MPI_ABORT|Fatal signal|terminate called'

run_command()
{
    local log_file=$1
    shift
    (
        cd "$RUNTIME_DIR" || exit 125
        OMP_NUM_THREADS=1 timeout --signal=TERM --kill-after=2s \
            "${TIMEOUT_SECONDS}s" "$@"
    ) >"$log_file" 2>&1
}

print_failure()
{
    local case_number=$1
    local name=$2
    local reason=$3
    local log_file=$4
    shift 4

    printf 'FAIL [%03d] %s: %s\n' "$case_number" "$name" "$reason" >&2
    printf '  command:' >&2
    printf ' %q' "$@" >&2
    printf '\n' >&2

    if [[ -s "$log_file" ]]; then
        local excerpt
        excerpt="$(grep -Ei '^(ERROR|WARNING)|Fix:|Adaptive .*selected|Unified convergence selected|Segmentation|segfault|Aborted|signal|done$' "$log_file" | tail -n 10)"
        if [[ -z "$excerpt" ]]; then
            excerpt="$(tail -n 8 "$log_file")"
        fi
        printf '  log excerpt:\n' >&2
        while IFS= read -r line; do
            printf '    %s\n' "$line" >&2
        done <<<"$excerpt"
    fi
}

record_pass()
{
    local case_number=$1
    local name=$2
    ((PASSED += 1))
    if [[ "${MBS_RELEASE_VERBOSE:-0}" == "1" ]]; then
        printf 'PASS [%03d] %s\n' "$case_number" "$name"
    fi
}

expect_success()
{
    local name=$1
    shift
    ((TOTAL += 1))
    ((RUN_TOTAL += 1))
    local case_number=$TOTAL
    local log_file
    local output_dir
    printf -v log_file '%s/%03d.log' "$LOG_DIR" "$case_number"
    printf -v output_dir '%s/%03d' "$RESULT_DIR" "$case_number"
    local -a command=("$@" --output "$output_dir")

    local status
    run_command "$log_file" "${command[@]}"
    status=$?

    local reason=''
    if ((status == 124)); then
        reason="timed out after ${TIMEOUT_SECONDS}s"
    elif ((status >= 128)); then
        reason="terminated by signal (status $status)"
    elif ((status != 0)); then
        reason="expected exit 0, got $status"
    elif grep -Eqi "$CRASH_PATTERN" "$log_file"; then
        reason='log contains an abort/crash signature'
    elif grep -Eq '(^|[[:space:]])ERROR:' "$log_file"; then
        reason='exit 0 log contains ERROR:'
    elif [[ -z "$(find "$output_dir" -type f -size +0c -print -quit 2>/dev/null)" ]]; then
        reason='calculation produced no nonempty output file'
    fi

    if [[ -n "$reason" ]]; then
        ((FAILED += 1))
        print_failure "$case_number" "$name" "$reason" "$log_file" \
            "${command[@]}"
        return
    fi

    ((RUN_PASSED += 1))
    record_pass "$case_number" "$name"
}

expect_error()
{
    local name=$1
    shift
    ((TOTAL += 1))
    ((ERROR_TOTAL += 1))
    local case_number=$TOTAL
    local log_file
    printf -v log_file '%s/%03d.log' "$LOG_DIR" "$case_number"
    local -a command=("$@")
    local have_output=0
    local argument
    for argument in "${command[@]}"; do
        if [[ "$argument" == '--output' || "$argument" == '-o' ]]; then
            have_output=1
            break
        fi
    done
    if ((have_output == 0)); then
        command+=(--output "$RESULT_DIR/error-$case_number")
    fi

    local status
    run_command "$log_file" "${command[@]}"
    status=$?

    local reason=''
    if ((status == 0)); then
        reason='invalid command unexpectedly succeeded'
    elif ((status == 124)); then
        reason="invalid command timed out after ${TIMEOUT_SECONDS}s"
    elif ((status >= 128)); then
        reason="invalid command terminated by signal (status $status)"
    elif grep -Eqi "$CRASH_PATTERN" "$log_file"; then
        reason='invalid command produced an abort/crash signature'
    elif ! grep -q 'ERROR' "$log_file"; then
        reason='diagnostic does not contain ERROR'
    elif ! grep -q 'Fix:' "$log_file"; then
        reason='diagnostic does not contain Fix:'
    fi

    if [[ -n "$reason" ]]; then
        ((FAILED += 1))
        print_failure "$case_number" "$name" "$reason" "$log_file" \
            "${command[@]}"
        return
    fi

    ((ERROR_PASSED += 1))
    record_pass "$case_number" "$name"
}

expect_scan_results()
{
    local name=$1
    local file_glob=$2
    local expected_count=$3
    local required_log_pattern=$4
    shift 4
    ((TOTAL += 1))
    ((RUN_TOTAL += 1))
    local case_number=$TOTAL
    local log_file output_dir
    printf -v log_file '%s/%03d.log' "$LOG_DIR" "$case_number"
    printf -v output_dir '%s/%03d' "$RESULT_DIR" "$case_number"
    local -a command=("$@" --output "$output_dir")

    local status
    run_command "$log_file" "${command[@]}"
    status=$?
    local reason=''
    local file_count=0
    if ((status == 0)); then
        file_count=$(find "$output_dir" -type f -name "$file_glob" -printf . 2>/dev/null \
            | wc -c)
    fi
    if ((status == 124)); then
        reason="timed out after ${TIMEOUT_SECONDS}s"
    elif ((status >= 128)); then
        reason="terminated by signal (status $status)"
    elif ((status != 0)); then
        reason="expected exit 0, got $status"
    elif grep -Eqi "$CRASH_PATTERN" "$log_file"; then
        reason='log contains an abort/crash signature'
    elif grep -Eq '(^|[[:space:]])ERROR:' "$log_file"; then
        reason='exit 0 log contains ERROR:'
    elif ((file_count != expected_count)); then
        reason="expected $expected_count files matching $file_glob, found $file_count"
    elif [[ -n "$required_log_pattern" ]] \
        && ! grep -Eq "$required_log_pattern" "$log_file"; then
        reason="log does not confirm applied values: $required_log_pattern"
    fi

    if [[ -n "$reason" ]]; then
        ((FAILED += 1))
        print_failure "$case_number" "$name" "$reason" "$log_file" \
            "${command[@]}"
        return
    fi
    ((RUN_PASSED += 1))
    record_pass "$case_number" "$name"
}

PHYSICS=(
    --particle 1 1 1
    --refractive-index 1.31 0
    --wavelength-um 10
    --max-reflections 2
)
PHYSICS_WITHOUT_DEPTH=(
    --particle 1 1 1
    --refractive-index 1.31 0
    --wavelength-um 10
)
GRID=(--scattering-grid 0 180 6 2)
PO_CORE=("$MBS" --method po "${PHYSICS[@]}" --backend cpu --threads 1 --close)
GO_CORE=("$MBS" --method go "${PHYSICS[@]}" --backend cpu --threads 1 --close)

ORIENTATION_NAMES=(
    fixed-orientation
    euler-grid
    monte-carlo
    orientation-file
    diffraction-grid
    sobol
    so3-quaternion
    sobol-seed
    sobol-ring
    hammersley
    lattice
    lattice-generator
    euler-quadrature
    euler-adaptive
    adaptive-euler-grid
    adaptive-orientations
    auto
    autofull
    diffraction-autofull
)

append_orientation()
{
    local index=$1
    local target_name=$2
    local -n target=$target_name
    case "$index" in
        0) target+=(--fixed-orientation 0 0) ;;
        1) target+=(--euler-grid 2 3) ;;
        2) target+=(--monte-carlo 4) ;;
        3) target+=(--orientation-file "$ORIENTATION_FILE") ;;
        4) target+=(--diffraction-grid 8) ;;
        5) target+=(--sobol 4) ;;
        6) target+=(--so3-quaternion 4) ;;
        7) target+=(--sobol-seed 4 7) ;;
        8) target+=(--sobol-ring 2 3) ;;
        9) target+=(--hammersley 4) ;;
        10) target+=(--lattice 5) ;;
        11) target+=(--lattice-generator 5 2) ;;
        12) target+=(--euler-quadrature 2 3) ;;
        13) target+=(--euler-adaptive 2 6) ;;
        14) target+=(--adaptive-euler-grid 0.9) ;;
        15) target+=(--adaptive-orientations 0.9) ;;
        16) target+=(--auto 0.9) ;;
        17) target+=(--autofull 0.9) ;;
        18) target+=(--diffraction-autofull 0.9) ;;
        *) return 2 ;;
    esac
}

printf 'Running real production calculations...\n'
for ((i = 0; i < ${#ORIENTATION_NAMES[@]}; ++i)); do
    command=("${PO_CORE[@]}")
    append_orientation "$i" command
    case "$i" in
        14|15)
            command+=("${GRID[@]}" --adaptive-config "$ADAPTIVE_CONFIG")
            ;;
        16|17|18)
            command+=(--adaptive-config "$ADAPTIVE_CONFIG")
            ;;
        *)
            command+=("${GRID[@]}")
            ;;
    esac
    expect_success "PO orientation: ${ORIENTATION_NAMES[i]}" "${command[@]}"
done

expect_success 'PO adaptive grid: auto-theta-grid' \
    "${PO_CORE[@]}" --sobol 16 --auto-theta-grid 0.9 \
    --adaptive-config "$ADAPTIVE_CONFIG"
expect_success 'PO adaptive grid: auto-phi' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --auto-phi \
    --adaptive-config "$ADAPTIVE_CONFIG"
expect_success 'PO adaptive grid: adaptive-phi' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --adaptive-phi 0.9 \
    --adaptive-config "$ADAPTIVE_CONFIG"
expect_success 'PO adaptive grid: adaptive-reflections' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --adaptive-reflections 0.9 \
    --adaptive-config "$ADAPTIVE_CONFIG"

GO_SUPPORTED=(0 1 2 4 5 7)
for i in "${GO_SUPPORTED[@]}"; do
    command=("${GO_CORE[@]}")
    append_orientation "$i" command
    command+=("${GRID[@]}")
    expect_success "GO orientation: ${ORIENTATION_NAMES[i]}" "${command[@]}"
done

expect_success 'theta grid: file' \
    "${PO_CORE[@]}" --sobol 4 --theta-grid-file "$THETA_FILE"
expect_success 'theta grid: backscatter cone' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 --scattering-grid 5 6 2
expect_success 'theta grid: production default' \
    "${PO_CORE[@]}" --fixed-orientation 0 0
expect_success 'backend: auto in CPU build' \
    "$MBS" --method po "${PHYSICS[@]}" --backend auto --threads 1 --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'multi-size: dmax grid with Sobol' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --dmax-grid 0.5 1 2
expect_success 'multi-size: k_eq grid with diffraction grid' \
    "${PO_CORE[@]}" --diffraction-grid 8 "${GRID[@]}" --k-eq-grid 0.1 0.2 2
expect_success 'multi-size: k_eq list with diffraction grid' \
    "${PO_CORE[@]}" --diffraction-grid 8 "${GRID[@]}" --k-eq-list "$K_EQ_FILE"
expect_scan_results 'multi-size semantics: Sobol Dmax values and files' \
    '*_x*.dat' 2 'x range: 0[.]314159.*0[.]628319' \
    "$MBS" --method po --particle 1 10 10 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --threads 1 --close \
    --sobol 4 "${GRID[@]}" --dmax-grid 1 2 2
expect_scan_results 'multi-size semantics: Euler-grid consumes Dmax scan' \
    '*_D*.dat' 2 'Euler-grid serial scan: Dmax 1[.]000000[.][.]2[.]000000' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" --dmax-grid 1 2 2
expect_scan_results 'multi-size semantics: Euler-grid consumes k_eq scan' \
    '*_keq*.dat' 2 'Euler-grid serial scan: k_eq 0p1[.][.]0p2' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" --k-eq-grid 0.1 0.2 2
expect_scan_results 'adaptive semantics: explicit theta and phi stay fixed' \
    '*_convergence.tsv' 1 'Unified convergence selected: n=2, N_phi=18, N_theta=3' \
    "${PO_CORE[@]}" --auto 0.9 "${GRID[@]}" --phi-points 18 \
    --adaptive-config "$ADAPTIVE_CONFIG"
expect_scan_results 'geometry semantics: beam-area ratio reaches Sobol tracing' \
    '*_log.txt' 1 'Beam area restriction ratio: 7[.]000000' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --beam-area-ratio 7
expect_success 'output: Jones matrices in fixed PO' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --jones-output
expect_success 'output: full and no-shadow PO results' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --no-shadow-output
expect_success 'particle file: strict valid cube' \
    "$MBS" --method po --particle-file "$PARTICLE_FILE" \
    --refractive-index 1.31 0 --wavelength-um 10 --max-reflections 2 \
    --backend cpu --threads 1 --close --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'particle file: packaged cube example' \
    "$MBS" --method po --particle-file "$EXAMPLE_PARTICLE" \
    --refractive-index 1.31 0 --wavelength-um 10 --max-reflections 2 \
    --backend cpu --threads 1 --close --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'particle: bullet' \
    "$MBS" --method po --particle 2 1 1 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --threads 1 --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'particle: bullet rosette' \
    "$MBS" --method po --particle 3 1 1 0.2 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --threads 1 --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'particle: canonical droxtal' \
    "$MBS" --method po --particle 4 1 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --threads 1 --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'particle: concave hexagonal' \
    "$MBS" --method po --particle 10 1 1 20 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --threads 1 --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'particle: two-column aggregate' \
    "$MBS" --method po --particle 12 1 1 2 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --threads 1 --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'particle: fixed built-in aggregate' \
    "$MBS" --method po --particle 999 0.01 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --threads 1 --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_success 'method modifier: Karczewski Monte Carlo' \
    "${PO_CORE[@]}" --monte-carlo 4 "${GRID[@]}" --karczewski
expect_success 'accuracy/performance and diagnostic controls' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --absorption --abs-points all --beam-area-ratio 100 \
    --beam-cutoff 0 --beam-cutoff-jones 0 --beam-cutoff-area 0 \
    --beam-cutoff-importance 0 --trace-cutoff 0 --trace-cutoff-jones 0 \
    --trace-cutoff-area 0 --trace-cutoff-importance 0 --trace-max-beams 0 \
    --trace-prefilter --trace-prefilter-margin 0 --trace-prefilter-stats \
    --progress-interval 0 --legacy-sign --full-only --shadow \
    --ot-phase-average --ot-phase-shift 0 --ot-ping-distance 0
expect_success 'tracing prefilter explicit disable' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --no-trace-prefilter
expect_success 'trajectory selection and grouping' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --trajectories "$TRAJECTORY_FILE" --all-trajectories --trajectory-groups
expect_success 'multi-size: process-parallel dmax file scan' \
    "$MBS" --method po --particle-file "$PARTICLE_FILE" \
    --refractive-index 1.31 0 --wavelength-um 10 --max-reflections 1 \
    --backend cpu --close --fixed-orientation 0 0 "${GRID[@]}" \
    --dmax-grid 0.5 1 2 --scan-jobs 2 --scan-threads 1

printf 'Running pairwise selector conflicts...\n'

METHOD_NAMES=(canonical-method legacy-po legacy-go)
append_method()
{
    local index=$1
    local target_name=$2
    local -n target=$target_name
    case "$index" in
        0) target+=(--method po) ;;
        1) target+=(--po) ;;
        2) target+=(--go) ;;
    esac
}
for ((i = 0; i < ${#METHOD_NAMES[@]}; ++i)); do
    for ((j = i + 1; j < ${#METHOD_NAMES[@]}; ++j)); do
        command=("$MBS" "${PHYSICS[@]}" --backend cpu --threads 1 --close \
            --fixed-orientation 0 0 "${GRID[@]}")
        append_method "$i" command
        append_method "$j" command
        expect_error "method selectors: ${METHOD_NAMES[i]} + ${METHOD_NAMES[j]}" \
            "${command[@]}"
    done
done
expect_error 'method selectors: duplicate canonical values' \
    "$MBS" --method po --method go "${PHYSICS[@]}" --backend cpu --close \
    --fixed-orientation 0 0 "${GRID[@]}"

BACKEND_NAMES=(canonical-backend legacy-cpu legacy-gpu)
append_backend()
{
    local index=$1
    local target_name=$2
    local -n target=$target_name
    case "$index" in
        0) target+=(--backend cpu) ;;
        1) target+=(--cpu) ;;
        2) target+=(--gpu) ;;
    esac
}
for ((i = 0; i < ${#BACKEND_NAMES[@]}; ++i)); do
    for ((j = i + 1; j < ${#BACKEND_NAMES[@]}; ++j)); do
        command=("$MBS" --method po "${PHYSICS[@]}" --threads 1 --close \
            --fixed-orientation 0 0 "${GRID[@]}")
        append_backend "$i" command
        append_backend "$j" command
        expect_error "backend selectors: ${BACKEND_NAMES[i]} + ${BACKEND_NAMES[j]}" \
            "${command[@]}"
    done
done
expect_error 'backend selectors: duplicate canonical values' \
    "$MBS" --method po "${PHYSICS[@]}" --backend cpu --backend auto --close \
    --fixed-orientation 0 0 "${GRID[@]}"

GEOMETRY_NAMES=(canonical-geometry force-convex force-nonconvex)
append_geometry()
{
    local index=$1
    local target_name=$2
    local -n target=$target_name
    case "$index" in
        0) target+=(--geometry auto) ;;
        1) target+=(--force-convex) ;;
        2) target+=(--force-nonconvex) ;;
    esac
}
for ((i = 0; i < ${#GEOMETRY_NAMES[@]}; ++i)); do
    for ((j = i + 1; j < ${#GEOMETRY_NAMES[@]}; ++j)); do
        command=("${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}")
        append_geometry "$i" command
        append_geometry "$j" command
        expect_error "geometry selectors: ${GEOMETRY_NAMES[i]} + ${GEOMETRY_NAMES[j]}" \
            "${command[@]}"
    done
done

expect_error 'particle selectors: built-in + file' \
    "$MBS" --method po --particle 1 1 1 --particle-file /dev/null \
    --refractive-index 1.31 0 --wavelength-um 10 --backend cpu --close \
    --fixed-orientation 0 0 "${GRID[@]}"

for ((i = 0; i < ${#ORIENTATION_NAMES[@]}; ++i)); do
    for ((j = i + 1; j < ${#ORIENTATION_NAMES[@]}; ++j)); do
        command=("${PO_CORE[@]}" "${GRID[@]}")
        append_orientation "$i" command
        append_orientation "$j" command
        expect_error "orientation selectors: ${ORIENTATION_NAMES[i]} + ${ORIENTATION_NAMES[j]}" \
            "${command[@]}"
    done
done

THETA_NAMES=(scattering-grid theta-grid-file auto-theta-grid)
append_theta_grid()
{
    local index=$1
    local target_name=$2
    local -n target=$target_name
    case "$index" in
        0) target+=("${GRID[@]}") ;;
        1) target+=(--theta-grid-file "$THETA_FILE") ;;
        2) target+=(--auto-theta-grid 0.9) ;;
    esac
}
for ((i = 0; i < ${#THETA_NAMES[@]}; ++i)); do
    for ((j = i + 1; j < ${#THETA_NAMES[@]}; ++j)); do
        command=("${PO_CORE[@]}" --sobol 16)
        append_theta_grid "$i" command
        append_theta_grid "$j" command
        expect_error "theta selectors: ${THETA_NAMES[i]} + ${THETA_NAMES[j]}" \
            "${command[@]}"
    done
done

SCAN_NAMES=(dmax-grid k-eq-grid k-eq-list)
append_scan()
{
    local index=$1
    local target_name=$2
    local -n target=$target_name
    case "$index" in
        0) target+=(--dmax-grid 0.5 1 2) ;;
        1) target+=(--k-eq-grid 0.1 0.2 2) ;;
        2) target+=(--k-eq-list "$K_EQ_FILE") ;;
    esac
}
for ((i = 0; i < ${#SCAN_NAMES[@]}; ++i)); do
    for ((j = i + 1; j < ${#SCAN_NAMES[@]}; ++j)); do
        command=("${PO_CORE[@]}" --sobol 4 "${GRID[@]}")
        append_scan "$i" command
        append_scan "$j" command
        expect_error "size selectors: ${SCAN_NAMES[i]} + ${SCAN_NAMES[j]}" \
            "${command[@]}"
    done
done

expect_error 'phi selectors: auto-phi + adaptive-phi' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --auto-phi --adaptive-phi 0.9
expect_error 'particle scaling selectors: resize + k_eq' \
    "$MBS" --method po --particle-file /dev/null --resize-dmax-um 1 --k-eq 1 \
    --refractive-index 1.31 0 --wavelength-um 10 --backend cpu --close \
    --fixed-orientation 0 0 "${GRID[@]}"

printf 'Running invalid cross-domain combinations...\n'

expect_error 'backend/method: CUDA requested from CPU binary' \
    "$MBS" --method po "${PHYSICS[@]}" --backend cuda --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_error 'backend/method: FFT on CPU backend' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --fft
expect_error 'backend/method: FFT with GO' \
    "${GO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --fft
expect_error 'backend/method: GO rejects ignored incoherent mode' \
    "${GO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --incoherent
expect_error 'backend/method: GO rejects ignored output-beam cutoff' \
    "${GO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --beam-cutoff 0.01

GO_UNSUPPORTED=(3 6 8 9 10 11 12 13 14 15 16 17 18)
for i in "${GO_UNSUPPORTED[@]}"; do
    command=("${GO_CORE[@]}")
    append_orientation "$i" command
    command+=("${GRID[@]}")
    expect_error "method/orientation: GO + ${ORIENTATION_NAMES[i]}" \
        "${command[@]}"
done

expect_error 'grid/adaptive: fixed + auto-theta-grid' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 --auto-theta-grid 0.9
expect_error 'grid/adaptive: fixed + auto-phi' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --auto-phi
expect_error 'grid/adaptive: fixed + adaptive-phi' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --adaptive-phi 0.9
expect_error 'grid/adaptive: fixed + adaptive-reflections' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --adaptive-reflections 0.9
expect_error 'grid/adaptive: adaptive-phi + explicit phi-points' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --adaptive-phi 0.9 --phi-points 12
expect_error 'grid/adaptive: auto-phi + explicit phi-points' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --auto-phi --phi-points 12
expect_error 'adaptive: auto + adaptive-phi' \
    "${PO_CORE[@]}" --auto 0.9 --adaptive-phi 0.9
expect_error 'adaptive: auto + adaptive-reflections' \
    "${PO_CORE[@]}" --auto 0.9 --adaptive-reflections 0.9
expect_error 'adaptive: autofull + adaptive-phi' \
    "${PO_CORE[@]}" --autofull 0.9 --adaptive-phi 0.9
expect_error 'adaptive: autofull + adaptive-reflections' \
    "${PO_CORE[@]}" --autofull 0.9 --adaptive-reflections 0.9
expect_error 'adaptive: config without adaptive search' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --adaptive-config "$ADAPTIVE_CONFIG"
expect_error 'adaptive: max-orientations with fixed Sobol count' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --max-orientations 32
expect_error 'adaptive: Owen averaging outside autofull' \
    "${PO_CORE[@]}" --auto 0.9 --owen-average 2
expect_error 'adaptive: Owen seeds outside autofull' \
    "${PO_CORE[@]}" --adaptive-orientations 0.9 "${GRID[@]}" --owen-seeds 1 2
expect_error 'adaptive: ring-points with fixed orientation' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --ring-points 3
expect_error 'adaptive: pole modifier with Sobol' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --pole
expect_error 'adaptive: mirror-gamma with fixed orientation' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --mirror-gamma
expect_error 'adaptive: symmetry with fixed orientation' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --symmetry 2 6
expect_error 'adaptive: max-theta-points without theta search' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --max-theta-points 33
expect_error 'adaptive: max-phi-points without phi search' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --max-phi-points 24
expect_error 'adaptive: max-beta-points without Euler search' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --max-beta-points 4
expect_error 'adaptive: max-gamma-points without Euler search' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --max-gamma-points 12
expect_error 'adaptive: convergence-passes without adaptive search' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --convergence-passes 1
expect_error 'adaptive runtime limit remains actionable' \
    "$MBS" --method po "${PHYSICS_WITHOUT_DEPTH[@]}" --backend cpu --threads 1 --close \
    --adaptive-euler-grid 1e-12 "${GRID[@]}" --max-beta-points 4 \
    --max-gamma-points 12 --max-orientations 128 --convergence-passes 1

expect_error 'multi-size: scan-jobs without a scan' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --scan-jobs 1
expect_error 'multi-size: scan-threads without a scan' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --scan-threads 1
expect_error 'multi-size/backend: gpu-devices without parallel scan' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --gpu-devices 0
expect_error 'multi-size: shared batches with dmax scan' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --dmax-grid 0.5 1 2 \
    --scan-jobs 1 --shared-k-eq-batches
expect_error 'multi-size: shared batches without scan-jobs' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --k-eq-grid 0.1 0.2 2 \
    --shared-k-eq-batches
expect_error 'multi-size/backend: shared batches on CPU' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --k-eq-grid 0.1 0.2 2 \
    --scan-jobs 1 --shared-k-eq-batches
expect_error 'multi-size/backend: gpu-devices on CPU' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --k-eq-grid 0.1 0.2 2 \
    --scan-jobs 1 --gpu-devices 0
expect_error 'multi-size: batch ratio without shared batches' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --k-eq-grid 0.1 0.2 2 \
    --k-eq-batch-ratio 1.05
expect_error 'multi-size: single k_eq plus k_eq scan' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --k-eq 0.2 --k-eq-grid 0.1 0.2 2
expect_error 'multi-size: parallel dmax scan with built-in particle' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --dmax-grid 0.5 1 2 --scan-jobs 1
expect_error 'multi-size: dmax scan with fixed orientation' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --dmax-grid 0.5 1 2
expect_error 'multi-size: dmax scan with so3-quaternion' \
    "${PO_CORE[@]}" --so3-quaternion 4 "${GRID[@]}" --dmax-grid 0.5 1 2
expect_error 'multi-size: dmax scan with sobol-seed' \
    "${PO_CORE[@]}" --sobol-seed 4 7 "${GRID[@]}" --dmax-grid 0.5 1 2
expect_error 'multi-size: dmax scan with sobol-ring' \
    "${PO_CORE[@]}" --sobol-ring 2 3 "${GRID[@]}" --dmax-grid 0.5 1 2
expect_error 'multi-size: dmax scan with Hammersley' \
    "${PO_CORE[@]}" --hammersley 4 "${GRID[@]}" --dmax-grid 0.5 1 2
expect_error 'multi-size: dmax scan with adaptive-gamma Euler' \
    "${PO_CORE[@]}" --euler-adaptive 2 6 "${GRID[@]}" --dmax-grid 0.5 1 2
expect_error 'multi-size: k_eq scan with fixed Sobol rule' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --k-eq-grid 0.1 0.2 2

expect_error 'output/method: Jones output with GO' \
    "${GO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --jones-output
expect_error 'output/method: no-shadow output with GO' \
    "${GO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --no-shadow-output
expect_error 'output/orientation: save-betas with fixed orientation' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --save-betas
expect_error 'output/orientation: checkpoint with Sobol' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --checkpoint
expect_error 'output: empty path' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --output ''
expect_error 'output: no-shadow result after disabling shadow beam' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --no-shadow-beam --no-shadow-output
expect_error 'legacy: unavailable single-backscatter-point mode' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" --backscatter-point
expect_error 'backend: GPU trace prefilter from CPU binary' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --gpu-trace-prefilter
expect_error 'tracing prefilter: enable and disable together' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --trace-prefilter --no-trace-prefilter
expect_error 'absorption: samples while absorption is disabled' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --abs-points all
expect_error 'geometry: nonpositive beam-area ratio' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --beam-area-ratio 0
expect_error 'cutoff: trace value outside relative range' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --trace-cutoff 1.1
expect_error 'tracing: negative beam limit' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --trace-max-beams -1
expect_error 'tracing: negative prefilter margin' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --trace-prefilter-margin -1
expect_error 'output: negative progress interval' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --progress-interval -1
expect_error 'trajectory: all requires a trajectory file' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --all-trajectories
expect_error 'trajectory: groups require a trajectory file' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --trajectory-groups
expect_error 'trajectory: malformed group index is actionable' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --trajectories "$BAD_TRAJECTORY_FILE"
expect_error 'trajectory: missing group index is actionable' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --trajectories "$MISSING_TRAJECTORY_GROUP_FILE"
expect_error 'trajectory: oversized group is bounded' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --trajectories "$OVERSIZED_TRAJECTORY_GROUP_FILE"

printf 'Running release hardening cases...\n'

expect_error 'particle file: empty input is actionable and does not crash' \
    "$MBS" --method po --particle-file "$EMPTY_PARTICLE_FILE" \
    --refractive-index 1.31 0 --wavelength-um 10 --max-reflections 1 \
    --backend cpu --threads 1 --close --fixed-orientation 0 0 "${GRID[@]}"
expect_error 'particle file: malformed vertex is actionable and does not crash' \
    "$MBS" --method po --particle-file "$BAD_PARTICLE_FILE" \
    --refractive-index 1.31 0 --wavelength-um 10 --max-reflections 1 \
    --backend cpu --threads 1 --close --fixed-orientation 0 0 "${GRID[@]}"
expect_error 'particle: droxtal rejects ambiguous three-value form' \
    "$MBS" --method po --particle 4 1 1 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_error 'particle: rosette rejects zero explicit cap' \
    "$MBS" --method po --particle 3 1 1 0 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_error 'particle: concave hexagonal rejects flat cavity' \
    "$MBS" --method po --particle 10 1 1 0 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_error 'particle: concave hexagonal rejects intersecting cavities' \
    "$MBS" --method po --particle 10 1 1 60 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_error 'particle: type 12 rejects unsupported component count' \
    "$MBS" --method po --particle 12 1 1 3 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 1 --backend cpu --close \
    --fixed-orientation 0 0 "${GRID[@]}"

expect_error 'theta file: duplicate rows' \
    "${PO_CORE[@]}" --sobol 4 --theta-grid-file "$DUPLICATE_THETA_FILE"
expect_error 'theta grid: equal four-value endpoints' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 \
    --scattering-grid 180 180 6 2
expect_error 'theta grid: detector cell product overflow' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 \
    --scattering-grid 0 180 50000 50000
expect_error 'theta grid: phi override endpoint overflow' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --phi-points 2147483647
expect_error 'adaptive phi: candidate growth overflow' \
    "${PO_CORE[@]}" --sobol 16 "${GRID[@]}" --adaptive-phi 0.9 \
    --max-phi-points 2147483647
expect_error 'orientation grid: Euler product overflow' \
    "${PO_CORE[@]}" --euler-grid 50000 50000 "${GRID[@]}"
expect_error 'orientation grid: adaptive Euler limit product overflow' \
    "${PO_CORE[@]}" --adaptive-euler-grid 0.1 "${GRID[@]}" \
    --max-beta-points 50000 --max-gamma-points 50000
expect_error 'orientation grid: lattice generator must be coprime' \
    "${PO_CORE[@]}" --lattice-generator 64 32 "${GRID[@]}"
expect_error 'orientation range: beta below zero' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" --beta-range-deg -1 90
expect_error 'orientation range: beta above 180' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" --beta-range-deg 0 181
expect_error 'orientation range: gamma equal endpoints' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" --gamma-range-deg 20 20
expect_error 'orientation range: gamma spans over 360 degrees' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" --gamma-range-deg -1 361
expect_error 'orientation modifier: symmetry with GO Sobol' \
    "${GO_CORE[@]}" --sobol 4 "${GRID[@]}" --symmetry 2 6
expect_error 'orientation modifier: symmetry with PO Euler grid' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" --symmetry 2 6
expect_error 'orientation modifier: pole with Monte Carlo' \
    "${PO_CORE[@]}" --monte-carlo 4 "${GRID[@]}" --pole
expect_error 'orientation modifier: chunk with fixed orientation' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --orientation-chunk 4
expect_error 'orientation modifier: chunk with Monte Carlo' \
    "${PO_CORE[@]}" --monte-carlo 4 "${GRID[@]}" --orientation-chunk 4
expect_error 'orientation modifier: chunk with orientation file' \
    "${PO_CORE[@]}" --orientation-file "$ORIENTATION_FILE" "${GRID[@]}" \
    --orientation-chunk 4
expect_error 'orientation modifier: coherent orientations with Sobol' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --coherent-orientations
expect_error 'orientation modifier: coherent and incoherent conflict' \
    "${PO_CORE[@]}" --euler-grid 2 3 "${GRID[@]}" \
    --coherent-orientations --incoherent
expect_error 'orientation modifier: Jones output with Sobol' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --jones-output
expect_error 'orientation modifier: Jones output with incoherent fixed PO' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --jones-output --incoherent
expect_error 'orientation modifier: backscatter filter with Sobol' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --backscatter-filter-deg 5
expect_error 'orientation modifier: no-shadow beam with GO' \
    "${GO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" --no-shadow-beam
expect_error 'adaptive: unified auto plus redundant auto-phi' \
    "${PO_CORE[@]}" --auto 0.9 --auto-phi
expect_error 'adaptive: unified auto plus redundant auto-theta-grid' \
    "${PO_CORE[@]}" --auto 0.9 --auto-theta-grid 0.9
expect_error 'adaptive: Owen generated and explicit seeds together' \
    "${PO_CORE[@]}" --autofull 0.9 --owen-average 2 --owen-seeds 1 2
expect_error 'adaptive: duplicate explicit Owen seeds' \
    "${PO_CORE[@]}" --autofull 0.9 --owen-seeds 1 1

expect_error 'multi-size: single k_eq plus dmax scan' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --k-eq 0.2 \
    --dmax-grid 0.5 1 2
expect_error 'multi-size: single dmax plus k_eq scan' \
    "$MBS" --method po --particle-file "$PARTICLE_FILE" --resize-dmax-um 1 \
    --refractive-index 1.31 0 --wavelength-um 10 --max-reflections 1 \
    --backend cpu --close --sobol 4 "${GRID[@]}" --k-eq-grid 0.1 0.2 2
expect_error 'multi-size: threads and scan-threads conflict' \
    "$MBS" --method po --particle-file "$PARTICLE_FILE" \
    --refractive-index 1.31 0 --wavelength-um 10 --max-reflections 1 \
    --backend cpu --close --sobol 4 "${GRID[@]}" --dmax-grid 0.5 1 2 \
    --scan-jobs 1 --threads 1 --scan-threads 1
expect_error 'multi-size: excessive dmax scan count' \
    "${PO_CORE[@]}" --sobol 4 "${GRID[@]}" --dmax-grid 0.5 1 100001
expect_error 'backend: excessive OpenMP thread count' \
    "$MBS" --method po "${PHYSICS[@]}" --backend cpu --threads 65537 --close \
    --fixed-orientation 0 0 "${GRID[@]}"
expect_error 'tracing: reflection depth exceeds location history' \
    "$MBS" --method po --particle 1 1 1 --refractive-index 1.31 0 \
    --wavelength-um 10 --max-reflections 31 --backend cpu --threads 1 --close \
    --fixed-orientation 0 0 "${GRID[@]}"

D_UNSUPPORTED=(0 2 3 6 7 8 9 13 14 15 16 17)
for i in "${D_UNSUPPORTED[@]}"; do
    command=("${PO_CORE[@]}" "${GRID[@]}")
    append_orientation "$i" command
    command+=(--dmax-grid 0.5 1 2)
    expect_error "all serial dmax paths: ${ORIENTATION_NAMES[i]} rejected" \
        "${command[@]}"
done
K_UNSUPPORTED=(0 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17)
for i in "${K_UNSUPPORTED[@]}"; do
    command=("${PO_CORE[@]}" "${GRID[@]}")
    append_orientation "$i" command
    command+=(--k-eq-grid 0.1 0.2 2)
    expect_error "all serial k_eq paths: ${ORIENTATION_NAMES[i]} rejected" \
        "${command[@]}"
done

expect_error 'output: unterminated placeholder' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --output "$RESULT_DIR/bad%"
expect_error 'output: unknown placeholder' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --output "$RESULT_DIR/%wat_run"
expect_error 'output: placeholder references a missing option' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --output "$RESULT_DIR/%0k-eq_run"
expect_error 'output: self-referencing placeholder' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --output '%o_run'
expect_error 'output: trailing slash lacks a result prefix' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --output "$RESULT_DIR/"
expect_error 'output: unwritable parent is actionable' \
    "${PO_CORE[@]}" --fixed-orientation 0 0 "${GRID[@]}" \
    --output /proc/mbs-release-matrix/run

if [[ -n "$MPI_RUNNER" ]]; then
    expect_success 'MPI: two-rank Sobol calculation' \
        env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
        "$MPI_RUNNER" --oversubscribe -np 2 "$MBS" --method po \
        "${PHYSICS[@]}" --backend cpu --threads 1 --close --sobol 4 "${GRID[@]}"
    expect_error 'MPI: fixed orientation is rejected before output races' \
        env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
        "$MPI_RUNNER" --oversubscribe -np 2 "$MBS" --method po \
        "${PHYSICS[@]}" --backend cpu --threads 1 --close \
        --fixed-orientation 0 0 "${GRID[@]}"
    expect_error 'MPI: rank-zero output failure is collective' \
        env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
        "$MPI_RUNNER" --oversubscribe -np 2 "$MBS" --method po \
        "${PHYSICS[@]}" --backend cpu --threads 1 --close --sobol 4 \
        "${GRID[@]}" --output /proc/mbs-release-mpi/run
    expect_error 'MPI: late zero-volume error finalizes every rank' \
        env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
        "$MPI_RUNNER" --oversubscribe -np 2 "$MBS" --method po \
        --particle-file "$ZERO_VOLUME_PARTICLE_FILE" --k-eq 1 \
        --refractive-index 1.31 0 --wavelength-um 10 --max-reflections 1 \
        --backend cpu --threads 1 --close --sobol 4 "${GRID[@]}"
    expect_error 'multi-size: MPI ranks cannot each fork a scan' \
        env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
        "$MPI_RUNNER" --oversubscribe -np 2 "$MBS" --method po \
        --particle-file "$PARTICLE_FILE" --refractive-index 1.31 0 \
        --wavelength-um 10 --max-reflections 1 --backend cpu --close \
        --sobol 4 "${GRID[@]}" --dmax-grid 0.5 1 2 --scan-jobs 1
fi

((TOTAL += 1))
runtime_leaks=$(find "$RUNTIME_DIR" -maxdepth 1 -type f \
    \( -name 'particle_for_check.dat' -o -name 'beam_log.dat' \
       -o -name 'log1.txt' -o -name 'energy.dat' \) -print)
if [[ -n "$runtime_leaks" ]]; then
    ((FAILED += 1))
    printf 'FAIL [%03d] output contract: diagnostics leaked into working directory:\n%s\n' \
        "$TOTAL" "$runtime_leaks" >&2
else
    record_pass "$TOTAL" 'output contract: no unsolicited working-directory files'
fi

printf '\nRelease CLI matrix summary\n'
printf '  checks:          %d\n' "$TOTAL"
printf '  passed:          %d\n' "$PASSED"
printf '  failed:          %d\n' "$FAILED"
printf '  real runs:       %d/%d passed\n' "$RUN_PASSED" "$RUN_TOTAL"
printf '  expected errors: %d/%d passed\n' "$ERROR_PASSED" "$ERROR_TOTAL"

if ((FAILED != 0)); then
    printf 'RELEASE CLI MATRIX FAILED\n' >&2
    exit 1
fi

printf 'RELEASE CLI MATRIX PASSED\n'
