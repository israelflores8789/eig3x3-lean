# eig3x3 — project recipes. `just --list` shows everything.

py := "python"
parity_dir := "parity"

default: test

# ---- Lean ----

build:
    lake build

test:
    lake test

build_bench:
    lake build eig3x3_bench

build_cli:
    lake build eig3x3_cli

bulid_all: build build_bench build_cli

# ---- Python parity harness ----
# (after the pytest migration, these become: pytest parity/tests -q)

parity n="1000" seed="42":
    cd {{parity_dir}} && {{py}} compare.py --n {{n}} --seed {{seed}}

errata n="2000" seed="42":
    cd {{parity_dir}} && {{py}} errata.py {{n}} {{seed}}

golden:
    cd {{parity_dir}} && {{py}} golden.py

bench:
    cd {{parity_dir}} && {{py}} bench.py

# ---- hygiene ----

lint:
    ruff check {{parity_dir}}
    ruff format --check {{parity_dir}}

typecheck:
    pyrefly check {{parity_dir}}

# ---- the full local gate ----

ci: test lint typecheck parity