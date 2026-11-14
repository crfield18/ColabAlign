"""Keep the release dependency pins consistent across installation routes."""

import re
import tomllib
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]


def _read_project():
    with (ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def _conda_pins(path):
    """Read exact package/version entries, not general YAML or Jinja syntax."""
    pattern = re.compile(r"^\s*-\s+([\w-]+)\s*=\s*([\w.]+)\s*$")
    return dict(
        match.groups()
        for line in path.read_text(encoding="utf-8").splitlines()
        if (match := pattern.match(line))
    )


def _python_pins(specs):
    pins = {}
    for spec in specs:
        name, version = spec.split("==")  # reject unpinned or range dependencies
        pins["python-build" if name == "build" else name] = version
    return pins


def test_python_dependency_pins_match_conda():
    project = _read_project()
    env = _conda_pins(ROOT / "colabalign.yml")
    specs = project["project"]["dependencies"] + project["build-system"]["requires"]
    for extra in project["project"]["optional-dependencies"].values():
        specs += extra
    assert _python_pins(specs).items() <= env.items()
    assert project["project"]["requires-python"] == ">=3.12,<3.13"
    assert env["python"].startswith("3.12.")


def test_recipe_pins_match_reference_environment():
    env = _conda_pins(ROOT / "colabalign.yml")
    recipe = _conda_pins(ROOT / "conda-recipe/colabalign/meta.yaml")
    assert recipe
    assert recipe.items() <= env.items()
    assert set(env) - {"pytest", "python-build"} == set(recipe)


def test_python_lock_has_hashes_and_matches_project_pins():
    project = _read_project()
    specs = project["project"]["dependencies"] + project["build-system"]["requires"]
    for extra in project["project"]["optional-dependencies"].values():
        specs += extra
    lines = (ROOT / "locks/requirements-py312.txt").read_text(encoding="utf-8").splitlines()
    resolved_specs = []
    for index, line in enumerate(lines):
        match = re.match(r"^([\w-]+==[\w.]+)", line)
        if match:
            resolved_specs.append(match.group(1))
            assert re.fullmatch(
                r"\s+--hash=sha256:[0-9a-f]{64}(?: \\)?", lines[index + 1]
            )
    assert _python_pins(specs).items() <= _python_pins(resolved_specs).items()


def test_explicit_locks_have_hashes_and_match_all_direct_pins():
    env = _conda_pins(ROOT / "colabalign.yml")
    for platform in ("linux-64", "osx-64", "osx-arm64"):
        lock = ROOT / f"locks/conda-{platform}.lock"
        lines = lock.read_text(encoding="utf-8").splitlines()
        assert "@EXPLICIT" in lines
        resolved = {}
        for line in lines:
            if not line or line.startswith(("#", "@")):
                continue
            url = urlparse(line)
            assert url.scheme == "https"
            assert url.netloc == "conda.anaconda.org"
            assert re.fullmatch(r"[0-9a-f]{32}", url.fragment)
            filename = Path(url.path).name.removesuffix(".conda").removesuffix(".tar.bz2")
            name, version, _build = filename.rsplit("-", 2)
            resolved[name] = version
        assert env.items() <= resolved.items()
