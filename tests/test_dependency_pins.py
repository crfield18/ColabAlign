"""Keep the release dependency pins consistent across installation routes."""

import re
import tomllib
from pathlib import Path
from urllib.parse import urlparse

from packaging.specifiers import SpecifierSet
from packaging.version import Version

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


def _python_ranges(specs):
    """Split runtime requirements into {name: specifier}, rejecting exact pins."""
    ranges = {}
    for spec in specs:
        name, specifier = re.fullmatch(r"([\w-]+)(.*)", spec).groups()
        assert "==" not in specifier, f"{name} is exactly pinned"
        ranges[name] = specifier
    return ranges


def _tooling_specs(project):
    specs = list(project["build-system"]["requires"])
    for extra in project["project"]["optional-dependencies"].values():
        specs += extra
    return specs


def test_python_dependency_pins_match_conda():
    project = _read_project()
    env = _conda_pins(ROOT / "colabalign.yml")
    assert _python_pins(_tooling_specs(project)).items() <= env.items()
    for name, specifier in _python_ranges(project["project"]["dependencies"]).items():
        assert Version(env[name]) in SpecifierSet(specifier), (name, specifier, env[name])
    assert project["project"]["requires-python"] == ">=3.14"
    assert env["python"].startswith("3.14.")


def _recipe_requirements(path):
    """Read recipe host and run requirements as {section: {name: specifier}}."""
    sections, current = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if match := re.match(r"^  (host|run):\s*$", line):
            current = sections.setdefault(match.group(1), {})
        elif current is not None and (match := re.match(r"^    -\s+([\w-]+)\s*(\S*)\s*$", line)):
            current[match.group(1)] = match.group(2)
        elif not line.startswith(("    ", "  #")):
            current = None
    return sections


def test_recipe_ranges_contain_reference_environment():
    """The recipe uses ranges; every exact reference pin must satisfy them."""
    env = _conda_pins(ROOT / "colabalign.yml")
    recipe = _recipe_requirements(ROOT / "conda-recipe/colabalign/meta.yaml")
    assert set(recipe) == {"host", "run"}
    for section, requirements in recipe.items():
        for name, spec in requirements.items():
            assert not re.fullmatch(r"=?=?[\w.]+", spec), f"{section} {name} is exactly pinned"
            assert Version(env[name]) in SpecifierSet(spec), (section, name, spec, env[name])
    assert set(env) - {"pytest", "python-build", "pip", "setuptools", "wheel"} == set(recipe["run"])
    project = _read_project()
    for name, specifier in _python_ranges(project["project"]["dependencies"]).items():
        assert SpecifierSet(recipe["run"][name]) == SpecifierSet(specifier), name
    requires_python = SpecifierSet(project["project"]["requires-python"])
    assert SpecifierSet(recipe["run"]["python"]) == requires_python
    assert SpecifierSet(recipe["host"]["python"]) == requires_python


def test_python_lock_has_hashes_and_matches_project_pins():
    project = _read_project()
    env = _conda_pins(ROOT / "colabalign.yml")
    lines = (ROOT / "locks/requirements-py314.txt").read_text(encoding="utf-8").splitlines()
    resolved_specs = []
    for index, line in enumerate(lines):
        match = re.match(r"^([\w-]+==[\w.]+)", line)
        if match:
            resolved_specs.append(match.group(1))
            assert re.fullmatch(
                r"\s+--hash=sha256:[0-9a-f]{64}(?: \\)?", lines[index + 1]
            )
    resolved = _python_pins(resolved_specs)
    assert _python_pins(_tooling_specs(project)).items() <= resolved.items()
    # The lock holds the exact reference versions of the ranged runtime dependencies.
    for name in _python_ranges(project["project"]["dependencies"]):
        assert resolved[name] == env[name], name


def test_reference_constraints_match_conda_pins():
    project = _read_project()
    env = _conda_pins(ROOT / "colabalign.yml")
    lines = (ROOT / "locks/reference-constraints.txt").read_text(encoding="utf-8").splitlines()
    constraints = _python_pins(line for line in lines if line and not line.startswith("#"))
    assert set(constraints) == set(_python_ranges(project["project"]["dependencies"]))
    assert constraints.items() <= env.items()


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
