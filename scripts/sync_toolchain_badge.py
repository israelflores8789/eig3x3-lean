"""Keep the README Lean badge synchronized with the root toolchain."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent.parent
TOOLCHAIN_FILE = ROOT / "lean-toolchain"
README_FILE = ROOT / "README.md"
BADGE_PATTERN = re.compile(r'https://img\.shields\.io/badge/Lean-[^" ]+\.svg')


def main() -> None:
    toolchain = TOOLCHAIN_FILE.read_text(encoding="utf-8").strip()
    prefix = "leanprover/lean4:"
    if not toolchain.startswith(prefix):
        raise SystemExit(f"unsupported lean-toolchain value: {toolchain!r}")

    version = toolchain.removeprefix(prefix)
    if not version.startswith("v") or any(char.isspace() for char in version):
        raise SystemExit(f"invalid Lean version in lean-toolchain: {version!r}")

    # Shields uses doubled hyphens to represent literal hyphens in badge text.
    badge_url = f"https://img.shields.io/badge/Lean-{version.replace('-', '--')}-blue.svg"
    readme = README_FILE.read_text(encoding="utf-8")
    updated, replacements = BADGE_PATTERN.subn(badge_url, readme)
    if replacements != 1:
        raise SystemExit(
            f"expected one Lean Shields badge in README.md, found {replacements}"
        )

    if updated == readme:
        print(f"README Lean badge is already synchronized with {version}")
        return

    README_FILE.write_text(updated, encoding="utf-8")
    print(f"updated README Lean badge to {version}")


if __name__ == "__main__":
    main()
