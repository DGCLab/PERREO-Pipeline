# ==============================
# Logging functions
# ==============================

if [[ -t 1 ]]; then
  RED="\e[31m"
  GREEN="\e[32m"
  YELLOW="\e[33m"
  BLUE="\e[34m"
  BOLD="\e[1m"
  RESET="\e[0m"

  OK="✔"
  WARN="⚠"
  ERR="✖"
else
  RED="" GREEN="" YELLOW="" BLUE="" BOLD="" RESET=""
  OK="[OK]"
  WARN="[WARN]"
  ERR="[ERROR]"
fi

msg_info() {
  echo -e "${BLUE}${BOLD}$*${RESET}"
}

msg_ok() {
  echo -e "${GREEN}${BOLD}[✔] $*${RESET}"
}

msg_warn() {
  echo -e "${YELLOW}${BOLD}[⚠] $*${RESET}"
}

msg_error() {
  echo -e "${RED}${BOLD}[✖] $*${RESET}" >&2
}

