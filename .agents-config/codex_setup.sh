#!/usr/bin/env bash

set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"

SKILLS="$ROOT/.agents/skills"
EM="$ROOT/.agent-config/emw"
NLA="$ROOT/.agent-config/nla"

mkdir -p "$SKILLS"

for skill in \
    electromagnetics-notation \
    surface-integral-equations \
    volume-integral-equations \
    operator-discretization \
    uniform-grid-vie \
    structured-em-operators
do
    ln -sfn \
        "../../.agent-config/electromagnetics/$skill" \
        "$SKILLS/$skill"
done

for skill in \
    linear-algebra-backends \
    structured-matrices \
    block-linear-algebra \
    performance-engineering \
    numerical-validation \
    toeplitz-preconditioning \
    sparse-direct-solvers
do
    ln -sfn \
        "../../.agent-config/nla/$skill" \
        "$SKILLS/$skill"
done
