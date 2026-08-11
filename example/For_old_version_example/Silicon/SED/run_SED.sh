#!/usr/bin/env bash
set -euo pipefail

input_file="input.in"
run_input="input.run.in"
trap 'rm -f "${run_input}"' EXIT

run_action() {
    local action="$1"
    sed -E "s/^action[[:space:]]*=.*/action = ${action}/" "${input_file}" > "${run_input}"
    pySED "${run_input}"
}

run_action compute
run_action plot
