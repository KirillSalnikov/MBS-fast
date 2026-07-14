#!/usr/bin/env bash
set -euo pipefail

nvcc="${1:-nvcc}"
preferred="${2:-g++}"

resolve_command() {
    if [[ "$1" == */* ]]; then
        [[ -x "$1" ]] && readlink -f "$1"
    else
        command -v "$1" 2>/dev/null
    fi
}

nvcc_path="$(resolve_command "$nvcc" || true)"
preferred_path="$(resolve_command "$preferred" || true)"
if [[ -z "$preferred_path" ]]; then
    preferred_path="$preferred"
fi

if [[ -z "$nvcc_path" ]]; then
    printf '%s\n' "$preferred_path"
    exit 0
fi

host_config="$(dirname "$nvcc_path")/../include/crt/host_config.h"
max_gcc=""
if [[ -r "$host_config" ]]; then
    max_gcc="$(sed -nE \
        's/^[[:space:]]*#if[[:space:]]+__GNUC__[[:space:]]*>[[:space:]]*([0-9]+).*/\1/p' \
        "$host_config" | head -1)"
fi

if [[ -z "$max_gcc" ]]; then
    printf '%s\n' "$preferred_path"
    exit 0
fi

compiler_major() {
    "$1" -dumpfullversion -dumpversion 2>/dev/null | cut -d. -f1
}

preferred_major="$(compiler_major "$preferred_path" || true)"
if [[ "$preferred_major" =~ ^[0-9]+$ ]] && (( preferred_major <= max_gcc )); then
    printf '%s\n' "$preferred_path"
    exit 0
fi

for ((version = max_gcc; version >= 5; --version)); do
    candidate="$(resolve_command "g++-$version" || true)"
    if [[ -n "$candidate" ]]; then
        printf '%s\n' "$candidate"
        exit 0
    fi
done

printf 'MBS-fast: CUDA supports GCC <= %s, but no compatible g++ was found.\n' \
    "$max_gcc" >&2
printf '%s\n' "$preferred_path"
