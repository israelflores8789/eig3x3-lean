# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""conftest.py — Pytest configuration, fixtures, and auto-scaffolding."""

import os
import subprocess

import numpy as np
import pytest

import golden
import lean_cli

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--n",
        action="store",
        type=int,
        default=1000,
        help="Batch size for random/logscale matrix suites",
    )
    parser.addoption(
        "--seed",
        action="store",
        type=int,
        default=42,
        help="Random seed for reproducible suite generation",
    )
    parser.addoption(
        "--lean-binary",
        action="store",
        default=None,
        help="Path to eig3x3_cli binary (defaults to auto-scaffolded Lake output)",
    )
    parser.addoption(
        "--no-auto-build",
        action="store_true",
        default=False,
        help="Disable automatic lake build of eig3x3_cli if missing",
    )


def pytest_sessionstart(session: pytest.Session) -> None:
    """Pre-scaffold CLI binary and golden artifacts on the controller process."""
    if not hasattr(session.config, "workerinput"):
        custom_bin = session.config.getoption("--lean-binary")
        binary = custom_bin if custom_bin else lean_cli.DEFAULT_BINARY
        no_auto_build = session.config.getoption("--no-auto-build")

        if not lean_cli.available(binary) and not no_auto_build:
            try:
                subprocess.run(
                    ["lake", "build", "eig3x3_cli"],
                    cwd=PROJECT_ROOT,
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except subprocess.CalledProcessError as e:
                pytest.fail(f"Failed to auto-build eig3x3_cli via Lake:\n{e.stderr}")

        if not os.path.isfile(golden.GOLDEN_JSON):
            golden.generate(golden.GOLDEN_JSON)


@pytest.fixture(scope="session")
def lean_cli_binary(request: pytest.FixtureRequest) -> str:
    """Auto-scaffolding fixture: ensures the Lean CLI binary is compiled."""
    custom_bin = request.config.getoption("--lean-binary")
    binary = custom_bin if custom_bin else lean_cli.DEFAULT_BINARY
    no_auto_build = request.config.getoption("--no-auto-build")

    if not lean_cli.available(binary):
        if no_auto_build:
            pytest.fail(
                f"eig3x3_cli not found at {binary!r} and --no-auto-build specified."
            )
        try:
            subprocess.run(
                ["lake", "build", "eig3x3_cli"],
                cwd=PROJECT_ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Failed to auto-build eig3x3_cli via Lake:\n{e.stderr}")

    if not lean_cli.available(binary):
        pytest.fail(f"eig3x3_cli still unavailable at {binary!r} after build.")

    return binary


@pytest.fixture(scope="session")
def golden_json() -> str:
    """Ensures generated/golden.json exists, generating it if absent."""
    path = golden.GOLDEN_JSON
    if not os.path.isfile(path):
        golden.generate(path)
    return path


@pytest.fixture
def batch_n(request: pytest.FixtureRequest) -> int:
    """Batch size configured via --n (default 1000)."""
    return int(request.config.getoption("--n"))


@pytest.fixture
def seed(request: pytest.FixtureRequest) -> int:
    """Random seed configured via --seed (default 42)."""
    return int(request.config.getoption("--seed"))


@pytest.fixture
def rng(seed: int) -> np.random.Generator:
    """Seeded numpy random Generator."""
    return np.random.default_rng(seed)
