#!/usr/bin/env bash
# Service launcher only. All numerical decisions remain inside --fullauto.
set -euo pipefail
if [[ $# -lt 4 ]]; then
    printf '%s\n' 'Usage: run_fullauto_service.sh UNIT BINARY LOG --fullauto EPS [native options...]' >&2
    exit 2
fi
task_unit=$1
task_binary=$(realpath "$2")
task_log=$(realpath -m "$3")
shift 3
[[ $task_unit =~ ^[a-zA-Z0-9_.-]+$ ]] || { printf '%s\n' 'Invalid unit name' >&2; exit 2; }
task_unit=${task_unit%.service}
task_directory=${XDG_CONFIG_HOME:-$HOME/.config}/systemd/user
mkdir -p "$task_directory" "$(dirname "$task_log")"
task_quote() {
    local task_value=$1
    task_value=${task_value//\\/\\\\}
    task_value=${task_value//\"/\\\"}
    task_value=${task_value//%/%%}
    task_value=${task_value//\$/\$\$}
    task_value=${task_value//$'\n'/\\n}
    printf '"%s"' "$task_value"
}
task_temporary=$(mktemp)
trap 'rm -f "$task_temporary"' EXIT
{
    printf '%s\n' '[Unit]' 'Description=MBS-fast native fullauto with durable checkpoints' 'StartLimitIntervalSec=3600' 'StartLimitBurst=3' '[Service]' 'Type=simple'
    for task_name in LD_LIBRARY_PATH CUDA_VISIBLE_DEVICES OMP_NUM_THREADS OMP_PROC_BIND OMP_PLACES; do
        if [[ -n ${!task_name:-} ]]; then
            printf 'Environment='; task_quote "$task_name=${!task_name}"; printf '\n'
        fi
    done
    printf 'WorkingDirectory=%s\nExecStart=' "${PWD//%/%%}";task_quote "$task_binary"
    for task_argument in "$@"; do printf ' ';task_quote "$task_argument";done
    printf '\n';printf '%s\n' 'Restart=on-failure' 'RestartPreventExitStatus=2 3 4' 'RestartSec=30' 'TimeoutStopSec=60' 'KillMode=control-group'
    printf 'StandardOutput=append:%s\n' "${task_log//%/%%}"
    printf '%s\n' 'StandardError=inherit' '[Install]' 'WantedBy=default.target'
} > "$task_temporary"
task_path=$task_directory/$task_unit.service
if [[ -e $task_path ]] && ! cmp -s "$task_path" "$task_temporary"; then
    printf 'Existing unit differs: %s; use a new unit name.\n' "$task_path" >&2
    exit 2
fi
task_user=$(id -un)
if [[ $(loginctl show-user "$task_user" -p Linger --value) != yes ]]; then
    loginctl enable-linger "$task_user"
fi
cp "$task_temporary" "$task_path"
systemctl --user daemon-reload
systemctl --user enable --now "$task_unit.service"
systemctl --user show "$task_unit.service" -p ActiveState -p MainPID -p UnitFileState
