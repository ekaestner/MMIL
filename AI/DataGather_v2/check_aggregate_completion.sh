#!/bin/bash
# check_aggregate_completion.sh
# Checks whether aggregate scripts have completed across multiple BIDS datasets.
# Scripts checked: fsqc.py, gather_cat12.py, gather_euler.py,
#                  gather_freesurfer.py, mriqc_aggregate.py

# ── Dataset list (edit as needed) ────────────────────────────
BIDS_ROOT="/space/mcdonald-syn01/1/BIDS"
DATASETS=(
    "ADNI-3"
    "BrACE"
    "CAPES"
    "DLBS_dataset"
    "ECP_dataset_bids"
    "ENIGMA_dataset"
    "HEP_dataset_bids"
    "IDEAS_dataset"
    "MTL_proj"
)

# ── Column labels ─────────────────────────────────────────────
COL_LABELS=("mriqc" "fsqc" "euler" "FS" "cat12")
N_COLS=${#COL_LABELS[@]}

# ── Colors & symbols (disabled when stdout is not a terminal) ──
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[1;33m'
    CYAN='\033[0;36m'
    BOLD='\033[1m'
    NC='\033[0m'
else
    GREEN=''; RED=''; YELLOW=''; CYAN=''; BOLD=''; NC=''
fi

# ── Counters for detailed section ────────────────────────────
PASS=0
FAIL=0

pass()  { echo -e "    ${GREEN}[PASS]${NC} $1"; ((PASS++)); }
fail()  { echo -e "    ${RED}[FAIL]${NC} $1"; ((FAIL++)); }
info()  { echo -e "    ${YELLOW}[INFO]${NC} $1"; }

# ── Summary table storage (flat array, row-major) ─────────────
# Values: 1=pass, 0=fail, -1=n/a (no relevant derivative dir)
declare -a SUMMARY   # SUMMARY[row*N_COLS + col]

# ── Helpers ───────────────────────────────────────────────────
latest_datedir() {
    ls -1 "${1}" 2>/dev/null | grep -E '^[0-9]{8}$' | sort -n | tail -1
}

# ── Per-script check functions ────────────────────────────────
# Each prints detail lines and echoes the result code: 1, 0, or -1

check_mriqc() {
    local derivdir="$1"
    local mriqc_dir
    mriqc_dir=$(ls -1d "${derivdir}"/mriqc* 2>/dev/null | head -1)

    if [[ -z "${mriqc_dir}" ]]; then
        fail "mriqc: No mriqc* directory found"
        echo -1; return
    fi
    info "mriqc: ${mriqc_dir}"

    local group_files
    group_files=$(find "${mriqc_dir}" -maxdepth 1 -name "group_*.tsv" -size +0c 2>/dev/null)
    if [[ -z "${group_files}" ]]; then
        fail "mriqc: No non-empty group_*.tsv files found"
        echo 0; return
    fi

    local ok=1
    while IFS= read -r f; do
        local rows
        rows=$(tail -n +2 "${f}" | wc -l)
        if [[ ${rows} -gt 0 ]]; then
            pass "mriqc: $(basename "${f}") — ${rows} rows"
        else
            fail "mriqc: $(basename "${f}") is empty"
            ok=0
        fi
    done <<< "${group_files}"
    echo ${ok}
}

check_fsqc() {
    local derivdir="$1"
    local fs_dir
    fs_dir=$(ls -1d "${derivdir}"/freesurfer* 2>/dev/null | head -1)

    if [[ -z "${fs_dir}" ]]; then
        fail "fsqc: No freesurfer* directory found"
        echo -1; return
    fi

    local fsqc_base="${fs_dir}/fsqc"
    local latest
    latest=$(latest_datedir "${fsqc_base}")
    if [[ -z "${latest}" ]]; then
        fail "fsqc: No date-stamped output in ${fsqc_base}"
        echo 0; return
    fi
    info "fsqc: most recent run ${latest}"

    local found=0
    for tpType in crossSectional longitudinal longTemplate; do
        local tpdir="${fsqc_base}/${latest}/${tpType}"
        if [[ -d "${tpdir}" ]]; then
            local n
            n=$(ls "${tpdir}" 2>/dev/null | wc -l)
            if [[ ${n} -gt 0 ]]; then
                pass "fsqc: ${tpType} — ${n} item(s)"
                found=1
            else
                fail "fsqc: ${tpType} directory is empty"
            fi
        fi
    done
    [[ ${found} -eq 0 ]] && { fail "fsqc: no tpType output found under ${fsqc_base}/${latest}"; echo 0; return; }
    echo 1
}

check_euler() {
    local derivdir="$1"
    local fs_dir
    fs_dir=$(ls -1d "${derivdir}"/freesurfer* 2>/dev/null | head -1)

    if [[ -z "${fs_dir}" ]]; then
        fail "euler: No freesurfer* directory found"
        echo -1; return
    fi

    local euler_base="${fs_dir}/euler"
    local latest
    latest=$(latest_datedir "${euler_base}")
    if [[ -z "${latest}" ]]; then
        fail "euler: No date-stamped output in ${euler_base}"
        echo 0; return
    fi
    info "euler: most recent run ${latest}"

    local found=0
    for tpType in crossSectional longTemplate; do
        local csv="${euler_base}/${latest}/${tpType}/euler_${tpType}.csv"
        if [[ -f "${csv}" && -s "${csv}" ]]; then
            local rows
            rows=$(tail -n +2 "${csv}" | wc -l)
            pass "euler: ${tpType}/euler_${tpType}.csv — ${rows} rows"
            found=1
        elif [[ -f "${csv}" ]]; then
            fail "euler: ${tpType}/euler_${tpType}.csv exists but is empty"
        fi
    done
    [[ ${found} -eq 0 ]] && { fail "euler: no euler_*.csv found under ${euler_base}/${latest}"; echo 0; return; }
    echo 1
}

check_gather_fs() {
    local derivdir="$1"
    local fs_dir
    fs_dir=$(ls -1d "${derivdir}"/freesurfer* 2>/dev/null | head -1)

    if [[ -z "${fs_dir}" ]]; then
        fail "gather_fs: No freesurfer* directory found"
        echo -1; return
    fi

    local measure_base="${fs_dir}/Measure"
    local latest
    latest=$(latest_datedir "${measure_base}")
    if [[ -z "${latest}" ]]; then
        fail "gather_fs: No date-stamped output in ${measure_base}"
        echo 0; return
    fi
    info "gather_fs: most recent run ${latest}"

    local found=0
    for tpType in crossSectional longitudinal longTemplate; do
        local tpdir="${measure_base}/${latest}/${tpType}"
        if [[ -d "${tpdir}" ]]; then
            local n_long n_wide
            n_long=$(find "${tpdir}" -name "*_longTable.csv" -size +0c 2>/dev/null | wc -l)
            n_wide=$(find "${tpdir}" -name "*_wideTable.csv" -size +0c 2>/dev/null | wc -l)
            if [[ ${n_long} -gt 0 && ${n_wide} -gt 0 ]]; then
                pass "gather_fs: ${tpType} — ${n_long} long + ${n_wide} wide CSVs"
                found=1
            else
                fail "gather_fs: ${tpType} — missing tables (long=${n_long}, wide=${n_wide})"
            fi
        fi
    done
    [[ ${found} -eq 0 ]] && { fail "gather_fs: no tpType output found under ${measure_base}/${latest}"; echo 0; return; }
    echo 1
}

check_cat12() {
    local derivdir="$1"
    local cat_dirs
    cat_dirs=$(ls -1d "${derivdir}"/cat* 2>/dev/null)

    if [[ -z "${cat_dirs}" ]]; then
        fail "cat12: No cat* directory found"
        echo -1; return
    fi

    # Use the last (most version-advanced) cat directory
    local cat_dir
    cat_dir=$(echo "${cat_dirs}" | tail -1)
    info "cat12: ${cat_dir}"

    local measures_base="${cat_dir}/Measures"
    local latest
    latest=$(latest_datedir "${measures_base}")
    if [[ -z "${latest}" ]]; then
        fail "cat12: No date-stamped output in ${measures_base}"
        echo 0; return
    fi
    info "cat12: most recent run ${latest}"

    local outdir="${measures_base}/${latest}"
    local result=1

    local qc_csv="${outdir}/aggregate_report_QC_qualitymeasures.csv"
    if [[ -f "${qc_csv}" && -s "${qc_csv}" ]]; then
        local rows
        rows=$(tail -n +2 "${qc_csv}" | wc -l)
        pass "cat12: aggregate_report_QC_qualitymeasures.csv — ${rows} rows"
    else
        fail "cat12: aggregate_report_QC_qualitymeasures.csv missing or empty"
        result=0
    fi

    local n_long n_wide
    n_long=$(find "${outdir}" -name "aggregate_biovariates_*_long.csv" -size +0c 2>/dev/null | wc -l)
    n_wide=$(find "${outdir}" -name "aggregate_biovariates_*_wide.csv" -size +0c 2>/dev/null | wc -l)
    if [[ ${n_long} -gt 0 && ${n_wide} -gt 0 ]]; then
        pass "cat12: ROI aggregates — ${n_long} long + ${n_wide} wide atlas CSVs"
    else
        fail "cat12: ROI aggregate files missing (long=${n_long}, wide=${n_wide})"
        result=0
    fi

    echo ${result}
}

# ── Helper: render a summary cell ────────────────────────────
render_cell() {
    case "$1" in
        1)  printf "${GREEN}  ✓  ${NC}" ;;
        0)  printf "${RED}  ✗  ${NC}" ;;
        -1) printf "${YELLOW}  -  ${NC}" ;;
        *)  printf "  ?  " ;;
    esac
}

# ═════════════════════════════════════════════════════════════
# MAIN LOOP
# ═════════════════════════════════════════════════════════════

echo ""
echo -e "${BOLD}══════════════════════════════════════════════════════════════${NC}"
echo -e "${BOLD} BrACE Pipeline — Aggregate Script Completion Check${NC}"
echo -e " $(date)"
echo -e "${BOLD}══════════════════════════════════════════════════════════════${NC}"

declare -a CHECK_FNS=("check_mriqc" "check_fsqc" "check_euler" "check_gather_fs" "check_cat12")

for row_idx in "${!DATASETS[@]}"; do
    ds="${DATASETS[${row_idx}]}"
    derivdir="${BIDS_ROOT}/${ds}/derivatives"

    echo ""
    echo -e "${CYAN}${BOLD}[ ${ds} ]${NC}  (${derivdir})"

    if [[ ! -d "${derivdir}" ]]; then
        echo -e "    ${RED}derivatives/ directory not found — skipping${NC}"
        for col_idx in $(seq 0 $((N_COLS - 1))); do
            SUMMARY[$((row_idx * N_COLS + col_idx))]=-1
        done
        continue
    fi

    for col_idx in "${!CHECK_FNS[@]}"; do
        fn="${CHECK_FNS[${col_idx}]}"
        result=$("${fn}" "${derivdir}" 2>&1)
        # The last line of output is the numeric code; everything else is detail
        code=$(echo "${result}" | tail -1)
        detail=$(echo "${result}" | head -n -1)
        [[ -n "${detail}" ]] && echo "${detail}"
        SUMMARY[$((row_idx * N_COLS + col_idx))]="${code}"
    done
done

# ═════════════════════════════════════════════════════════════
# SUMMARY TABLE
# ═════════════════════════════════════════════════════════════

echo ""
echo -e "${BOLD}══════════════════════════════════════════════════════════════${NC}"
echo -e "${BOLD} SUMMARY TABLE${NC}"
echo -e " ${GREEN}✓${NC} = complete   ${RED}✗${NC} = incomplete   ${YELLOW}–${NC} = no relevant derivative"
echo -e "${BOLD}══════════════════════════════════════════════════════════════${NC}"
echo ""

# Column widths
DS_W=18   # dataset name width
COL_W=9   # each check column width

# Header
printf "%-${DS_W}s" "Dataset"
for label in "${COL_LABELS[@]}"; do
    printf "%-${COL_W}s" "${label}"
done
echo ""

# Separator
printf '%*s' "${DS_W}" '' | tr ' ' '-'
for label in "${COL_LABELS[@]}"; do
    printf '%*s' "${COL_W}" '' | tr ' ' '-'
done
echo ""

# Data rows
for row_idx in "${!DATASETS[@]}"; do
    ds="${DATASETS[${row_idx}]}"
    printf "%-${DS_W}s" "${ds}"
    for col_idx in $(seq 0 $((N_COLS - 1))); do
        code="${SUMMARY[$((row_idx * N_COLS + col_idx))]}"
        render_cell "${code}"
        # Pad to col width (render_cell prints 5 visible chars including color escapes)
        printf "%-$((COL_W - 5))s" ""
    done
    echo ""
done

echo ""

# Totals
total_pass=0; total_fail=0; total_na=0
for val in "${SUMMARY[@]}"; do
    case "${val}" in
        1)  ((total_pass++)) ;;
        0)  ((total_fail++)) ;;
        -1) ((total_na++)) ;;
    esac
done
total_checks=$(( total_pass + total_fail ))
echo -e " Across all datasets: ${GREEN}${total_pass} passed${NC}  ${RED}${total_fail} failed${NC}  ${YELLOW}${total_na} n/a${NC}  (${total_checks} applicable checks)"
echo -e "${BOLD}══════════════════════════════════════════════════════════════${NC}"
