#!/usr/bin/env bash

# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

# apply-rulesets.sh — create or update repository rulesets from
# .github/rulesets/*.json. Idempotent: matches existing rulesets by name
# and updates them in place, creating only what's missing.
#
# Requires: gh (authenticated as a repo admin) and jq.
# Usage: ./apply-rulesets.sh
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(gh repo view --json nameWithOwner --jq .nameWithOwner)"

for file in "$script_dir"/*.json; do
  name="$(jq --raw-output .name "$file")"
  id="$(gh api "repos/$repo/rulesets" \
    --jq ".[] | select(.name == \"$name\") | .id" \
    | head -n 1)"

  if [[ -n "$id" ]]; then
    gh api --method PUT "repos/$repo/rulesets/$id" --input "$file" > /dev/null
    echo "updated ruleset: $name (id $id)"
  else
    gh api --method POST "repos/$repo/rulesets" --input "$file" > /dev/null
    echo "created ruleset: $name"
  fi
done
