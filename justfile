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
    cd {{parity_dir}} && {{pytest}} parity.py {{ if suite == "" { "" } else { "-k " + suite } }} --n {{n}} --seed {{seed}} -v

parity-parallel suite="" n="1000" seed="42" workers="auto":
    cd {{parity_dir}} && {{pytest}} parity.py -n {{workers}} {{ if suite == "" { "" } else { "-k " + suite } }} --n {{n}} --seed {{seed}} -v

parity-ci suite="" n="1000" seed="42" workers="auto" xml="generated/junit.xml":
    cd {{parity_dir}} && {{pytest}} parity.py -n {{workers}} --junitxml={{xml}} {{ if suite == "" { "" } else { "-k " + suite } }} --n {{n}} --seed {{seed}} -v

errata n="2000" seed="42":
    cd {{parity_dir}} && {{py}} errata.py {{n}} {{seed}}

golden:
    cd {{parity_dir}} && {{py}} golden.py

bench:
    cd {{parity_dir}} && {{py}} bench.py

# ---- hygiene ----

lint:
    uv run ruff check {{parity_dir}}
    uv run ruff format --check {{parity_dir}}

typecheck:
    uv run pyrefly check {{parity_dir}}

# ---- the full local / CI gate ----

ci: test lint typecheck parity-ci