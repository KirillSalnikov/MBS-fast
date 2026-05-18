#!/usr/bin/env bash
set -euo pipefail

CXX="${1:-${CXX:-g++}}"

supports_flags() {
    local flags="$1"
    local obj="/tmp/mbs_arch_test_$$.o"
    printf 'int main(){return 0;}\n' \
        | "$CXX" -x c++ - $flags -c -o "$obj" >/dev/null 2>&1
    local rc=$?
    rm -f "$obj"
    return "$rc"
}

cpu_model="$(lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)"

case "$cpu_model" in
    *"EPYC 7H12"*)
        echo "-march=znver2 -mtune=znver2"
        ;;
    *"EPYC 9"*|*"Ryzen AI"*|*"Ryzen 9 99"*|*"Ryzen 9 98"*|*"Ryzen 7 98"*|*"Ryzen 7 97"*|*"Ryzen 5 96"*)
        if supports_flags "-march=znver5 -mtune=znver5 -mavx512f -mavx512dq -mavx512vl"; then
            echo "-march=znver5 -mtune=znver5 -mavx512f -mavx512dq -mavx512vl"
        elif supports_flags "-march=znver4 -mtune=znver4 -mavx512f -mavx512dq -mavx512vl"; then
            echo "-march=znver4 -mtune=znver4 -mavx512f -mavx512dq -mavx512vl"
        else
            echo "-march=native -mtune=native"
        fi
        ;;
    *)
        echo "-march=native -mtune=native"
        ;;
esac
