# eig3x3 — project recipes. `just --list` shows everything.

py := "uv run python"
pytest := "uv run pytest"
parity_dir := "parity"

default: build test lint typecheck parity

# ---- Lean ----

test:
    lake test

build-lean:
    lake build

build-bench:
    lake build eig3x3_bench

build-cli:
    lake build eig3x3_cli

build: build-lean build-bench build-cli

# ---- Python parity harness ----

parity suite="" n="1000" seed="42":
    {{pytest}} {{parity_dir}}/parity.py {{ if suite == "" { "" } else { "-k " + suite } }} --n {{n}} --seed {{seed}} -v

parity-parallel suite="" n="1000" seed="42" workers="auto":
    {{pytest}} {{parity_dir}}/parity.py -n {{workers}} {{ if suite == "" { "" } else { "-k " + suite } }} --n {{n}} --seed {{seed}} -v

parity-ci suite="" n="1000" seed="42" workers="auto" xml="generated/junit.xml":
    {{pytest}} {{parity_dir}}/parity.py -n {{workers}} --junitxml={{xml}} {{ if suite == "" { "" } else { "-k " + suite } }} --n {{n}} --seed {{seed}} -v

errata n="2000" seed="42":
    cd {{parity_dir}} && {{py}} errata.py {{n}} {{seed}}

golden:
    cd {{parity_dir}} && {{py}} golden.py

bench:
    cd {{parity_dir}} && {{py}} bench.py

# ---- documentation ----

docs-update:
    cd docbuild && lake update doc-gen4

docs-dev:
    cd docbuild && DOCGEN_SRC="vscode" lake build Eig3x3:docs

docs:
    cd docbuild && lake build Eig3x3:docs

# docs-serve port="8000":
#     cd docbuild/.lake/build/doc && {{py}} -m http.server {{port}}
docs-serve port="8000":
    #!/usr/bin/env bash
    if command -v lsof >/dev/null 2>&1; then
        pids=$(lsof -t -i :{{port}});
        if [ -n "$pids" ]; then
            echo "Port {{port}} is already in use by PID(s): $pids";
            echo "Run 'just docs-serve-stop {{port}}' first.";
            exit 1;
        fi;
    elif command -v fuser >/dev/null 2>&1; then
        if fuser {{port}}/tcp >/dev/null 2>&1; then
            echo "Port {{port}} is already in use.";
            echo "Run 'just docs-serve-stop {{port}}' first.";
            exit 1;
        fi;
    else
        if {{py}} -c "import socket; s = socket.socket(); s.bind(('', {{port}})); s.close()" 2>/dev/null; then
            true;
        else
            echo "Port {{port}} is already in use (detected via socket bind).";
            echo "Run 'just docs-serve-stop {{port}}' first.";
            exit 1;
        fi;
    fi
    cd docbuild/.lake/build/doc && {{py}} -m http.server {{port}}

docs-serve-stop port="8000":
    #!/usr/bin/env bash
    echo "Stopping http.server on port 8000..."
    if command -v lsof >/dev/null 2>&1; then
        pids=$(lsof -t -i :{{port}});
        if [ -n "$pids" ]; then
            echo "PIDs: $pids";
            kill $pids;
        else
            echo "Nothing running on port {{port}}";
        fi
    elif command -v fuser >/dev/null 2>&1; then
        fuser -k {{port}}/tcp || true;
    else
        echo "No lsof or fuser; please install one of them or kill manually.";
        exit 1;
    fi

# ---- hygiene ----

spell:
    typos

spell-diff:
    typos --diff

spell-fix:
    typos --write-changes

lint:
    uv run ruff check {{parity_dir}}
    uv run ruff format --check {{parity_dir}}

typecheck:
    uv run pyrefly check {{parity_dir}}

# ---- full CI gate ----

ci: test lint typecheck parity-ci

# --- release utilities ---

# All recipes are idempotent.

# Install jq and the GitHub CLI (Debian/Ubuntu), then print next steps.
release-tooling: release-jq release-gh
    @echo ""
    @echo "Tooling ready. If this was a fresh gh install, run: gh auth login"

# jq — JSON processor used by apply-rulesets.sh
release-jq:
    #!/usr/bin/env bash
    set -euo pipefail
    if command -v jq >/dev/null 2>&1; then
        echo "jq already installed: $(jq --version)"
    else
        apt-get update && apt-get install -y jq
    fi

# gh — GitHub CLI from Debian's repositories
# note: may need proper permissions with `chmod 1777 /tmp`
release-gh:
    #!/usr/bin/env bash
    set -euo pipefail
    if command -v gh >/dev/null 2>&1; then
        echo "gh already installed: $(gh --version | head -n 1)"
    else
        apt-get update && apt-get install -y gh
    fi

# Moves the CHANGELOG [Unreleased] section to the new version, bumps
# version in pyproject.toml + lakefile.toml, regenerates uv.lock,
# and runs pre-flight checks. The commit and tag stay manual and
# deliberate. Includes guards to prevent accidental double-release
# or empty release notes.
release bump="patch":
    #!/usr/bin/env bash
    set -euo pipefail

    # --- Guard 0: clean tree, so the bump commit stays atomic ---
    if [ -n "$(git status --porcelain)" ]; then
        echo "error: working tree is dirty; commit or stash first" >&2
        exit 1
    fi

    # --- Compute the candidate version WITHOUT writing pyproject.toml ---
    version="$(uv version --bump "{{ bump }}" --dry-run --short)"
    echo "candidate version: $version"

    # --- Guard 1: exactly one Unreleased heading (catches format drift) ---
    count="$(grep -cE '^## \[Unreleased\]$' CHANGELOG.md || true)"
    if [ "$count" -ne 1 ]; then
        echo "error: expected exactly one '## [Unreleased]' heading, found $count" >&2
        exit 1
    fi

    # --- Guard 2: idempotency — refuse if this version section exists ---
    if grep -qE "^## \[$version\] - " CHANGELOG.md; then
        echo "error: CHANGELOG.md already has a [$version] section" >&2
        exit 1
    fi

    # --- Guard 3: Unreleased must contain at least one bullet ---
    # (Doubles as a re-run guard: after a successful rotation, Unreleased
    # is empty, so running `just release` again aborts here.)
    bullets="$(awk '
        /^## \[Unreleased\]$/ { found=1; next }
        /^## \[/ { found=0 }
        found' CHANGELOG.md | grep -cE '^[[:space:]]*[-*] ' || true)"
    if [ "$bullets" -eq 0 ]; then
        echo "error: [Unreleased] has no bullet entries; write release notes first" >&2
        exit 1
    fi

    # --- All guards passed; now mutate, in dependency order ---

    # 1. Bump pyproject.toml (uv version writes the file itself).
    uv version --bump "{{ bump }}"

    # 2. Mirror into lakefile.toml (plain semver, no v! literal).
    sed -i -E "s/^(version[[:space:]]*=[[:space:]]*)\"[^\"]*\"/\1\"$version\"/" lakefile.toml
    grep -E "^version[[:space:]]*=" lakefile.toml

    sed -i -E "s/^(version[[:space:]]*=[[:space:]]*)\"[^\"]*\"/\1\"$version\"/" docbuild/lakefile.toml
    grep -E "^version[[:space:]]*=" docbuild/lakefile.toml

    # 3. Regenerate uv.lock against the bumped pyproject.
    uv sync

    # 4. CHANGELOG rotation via temp file; swap only after validation.
    today="$(date +%Y-%m-%d)"
    tmp="$(mktemp)"
    awk -v version="$version" -v date="$today" '
        /^## \[Unreleased\]$/ {
            print "## [Unreleased]"
            print ""
            print "## [" version "] - " date
            next
        }
        { print }
    ' CHANGELOG.md > "$tmp"

    if ! grep -qE "^## \[$version\] - $today\$" "$tmp"; then
        rm -f "$tmp"
        echo "error: CHANGELOG rewrite failed validation; file untouched" >&2
        exit 1
    fi
    mv "$tmp" CHANGELOG.md
    grep -A4 "^## \[$version\]" CHANGELOG.md | head -6

    echo ""
    echo "Bumped to $version. Review the diff, then:"
    echo "  git add -A && git commit -m \"Release v$version\""
    echo "  git push, open the dev -> main PR, merge, then:"
    echo "  git tag -a v$version -m \"v$version\" && git push origin v$version"


# Create or update GitHub rulesets from .github/rulesets/*.json.
# Requires gh authentication: run 'gh auth login' once first.
rulesets: release-tooling
    ./.github/rulesets/apply-rulesets.sh
