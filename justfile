# eig3x3 — project recipes. `just --list` shows everything.

py := "python"
parity_dir := "parity"

default: test

# ---- Lean ----

build:
    lake build

test:
    lake test

cli:
    lake build eig3x3-cli

# ---- Python parity harness ----
# (after the pytest migration, these become: pytest parity/tests -q)

parity n="1000" seed="42":
    cd {{parity_dir}} && {{py}} compare.py --impl mirror --n {{n}} --seed {{seed}}

exact n="2000":
    cd {{parity_dir}} && {{py}} exact.py {{n}}

props n="2000":
    cd {{parity_dir}} && {{py}} properties.py {{n}}

golden:
    cd {{parity_dir}} && {{py}} golden.py

golden-check:
    cd {{parity_dir}} && {{py}} golden.py --check

bench:
    cd {{parity_dir}} && {{py}} bench.py

# ---- hygiene ----

lint:
    ruff check {{parity_dir}}
    ruff format --check {{parity_dir}}

typecheck:
    pyrefly check {{parity_dir}}

# ---- the full local gate ----

ci: test lint typecheck exact props parity golden-check