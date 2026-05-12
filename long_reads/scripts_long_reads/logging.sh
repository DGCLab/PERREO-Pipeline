# ==============================
# PERREO terminal UI helpers
# ==============================

# Color is enabled only for interactive terminals by default. Set FORCE_COLOR=1
# to force ANSI colors, or NO_COLOR=1 to disable them explicitly.
if [[ -n "${NO_COLOR:-}" ]]; then
  PERREO_COLOR=0
elif [[ -t 1 || -n "${FORCE_COLOR:-}" ]]; then
  PERREO_COLOR=1
else
  PERREO_COLOR=0
fi

if [[ "$PERREO_COLOR" -eq 1 ]]; then
  RED="\e[31m"
  GREEN="\e[32m"
  YELLOW="\e[33m"
  BLUE="\e[34m"
  MAGENTA="\e[35m"
  CYAN="\e[36m"
  DIM="\e[2m"
  BOLD="\e[1m"
  RESET="\e[0m"
else
  RED="" GREEN="" YELLOW="" BLUE="" MAGENTA="" CYAN="" DIM="" BOLD="" RESET=""
fi

OK="✔"
WARN="⚠"
ERR="✖"
INFO="ℹ"
STEP="▶"
SKIP="↷"
DONE="✓"

_perreo_timestamp() {
  date '+%H:%M:%S'
}

_perreo_strip_ansi() {
  sed -E 's/\x1B\[[0-9;]*[[:alpha:]]//g'
}

_perreo_visible_length() {
  printf '%b' "$1" | _perreo_strip_ansi | awk '{ print length }'
}

_perreo_repeat() {
  local char="$1"
  local width="$2"
  local out=""
  local i
  for ((i=0; i<width; i++)); do
    out+="$char"
  done
  printf '%s' "$out"
}

_perreo_hr() {
  local char="${1:-─}"
  local width="${PERREO_WIDTH:-82}"
  _perreo_repeat "$char" "$width"
  printf '\n'
}

msg_info() {
  printf '%b\n' "${BLUE}${BOLD}${INFO} [$(_perreo_timestamp)]${RESET} $*"
}

msg_ok() {
  printf '%b\n' "${GREEN}${BOLD}${OK} [$(_perreo_timestamp)]${RESET} $*"
}

msg_warn() {
  printf '%b\n' "${YELLOW}${BOLD}${WARN} [$(_perreo_timestamp)]${RESET} $*"
}

msg_error() {
  printf '%b\n' "${RED}${BOLD}${ERR} [$(_perreo_timestamp)]${RESET} $*" >&2
}

msg_skip() {
  printf '%b\n' "${CYAN}${BOLD}${SKIP} [$(_perreo_timestamp)]${RESET} $*"
}

_perreo_logo() {
  if [[ "${PERREO_SHOW_LOGO:-1}" == "0" ]]; then
    return
  fi

  printf '%b' "${CYAN}${BOLD}"
  cat <<'PERREO_LOGO'
    ____  ______ ____   ____   ______ ____
   / __ \/ ____// __ \ / __ \ / ____// __ \
  / /_/ / __/  / /_/ // /_/ // __/  / / / /
 / ____/ /___ / _, _// _, _// /___ / /_/ /
/_/   /_____//_/ |_|/_/ |_|/_____/ \____/
PERREO_LOGO
  printf '%b\n' "${RESET}"
}

perreo_banner() {
  local mode="$1"
  local description="$2"
  local width="${PERREO_WIDTH:-82}"
  local line
  line=$(_perreo_repeat '═' "$width")

  _perreo_logo
  printf '%b\n' "${CYAN}${BOLD}╔${line}╗${RESET}"
  printf '%b\n' "${CYAN}${BOLD}║${RESET} ${BOLD}PERREO Pipeline${RESET} ${DIM}|${RESET} ${MAGENTA}${mode}${RESET}"
  printf '%b\n' "${CYAN}${BOLD}║${RESET} ${description}"
  printf '%b\n' "${CYAN}${BOLD}║${RESET} ${DIM}Coffee status:${RESET} above threshold ☕"
  printf '%b\n' "${CYAN}${BOLD}╚${line}╝${RESET}"
}

perreo_stage() {
  local current="$1"
  local total="$2"
  local title="$3"
  printf '\n%b\n' "${BLUE}${BOLD}${STEP} Stage ${current}/${total}:${RESET} ${BOLD}${title}${RESET}"
  printf '%b\n' "${DIM}$(_perreo_hr '─')${RESET}"
}

perreo_summary() {
  local title="$1"
  shift
  local width="${PERREO_WIDTH:-82}"
  local key value pair key_len value_len pad

  printf '\n%b\n' "${MAGENTA}${BOLD}┌─ ${title} ${RESET}"
  for pair in "$@"; do
    key="${pair%%=*}"
    value="${pair#*=}"
    key_len=$(_perreo_visible_length "$key")
    value_len=$(_perreo_visible_length "$value")
    pad=$(( width - key_len - value_len - 7 ))
    (( pad < 1 )) && pad=1
    printf '%b' "${MAGENTA}${BOLD}│${RESET} ${key} "
    _perreo_repeat '.' "$pad"
    printf '%b\n' " ${GREEN}${value}${RESET}"
  done
  printf '%b\n\n' "${MAGENTA}${BOLD}└$(_perreo_hr '─')${RESET}"
}

perreo_finished() {
  local message="${1:-PERREO successfully completed}"
  printf '\n%b\n' "${GREEN}${BOLD}${DONE} ${message}${RESET}"
}
