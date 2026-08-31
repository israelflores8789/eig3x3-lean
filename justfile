# eig3x3 — project recipes. `just --list` shows everything.

py := "uv run python"
pytest := "uv run pytest"
parity_dir := "parity"

default: build_all test lint typecheck parity

# ---- Lean ----

build:
    lake build

test:
    lake test

build_bench:
    lake build eig3x3_bench

build_cli:
    lake build eig3x3_cli

build_all: build build_bench build_cli

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

docs-build:
    cd docbuild && lake build Eig3x3:docs

docs: docs-build

docs-serve port="8000":
    cd docbuild/.lake/build/doc && {{py}} -m http.server {{port}}

# ---- hygiene ----

lint:
    uv run ruff check {{parity_dir}}
    uv run ruff format --check {{parity_dir}}

typecheck:
    uv run pyrefly check {{parity_dir}}

# ---- the full local / CI gate ----

ci: test lint typecheck parity-ci