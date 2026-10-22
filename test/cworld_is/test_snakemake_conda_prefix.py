from hashlib import md5
from pathlib import Path
from types import SimpleNamespace

from snakemake.deployment.conda import Env


def test_direct_mamba_prefix_matches_snakemake_environment_address(tmp_path: Path) -> None:
    spec = tmp_path / "cworld.yaml"
    spec.write_text("name: cworld\ndependencies:\n  - perl=5.32\n")
    prefix_dir = tmp_path / "conda-prefix"
    prefix_dir.mkdir()

    workflow = SimpleNamespace(
        deployment_settings=SimpleNamespace(apptainer_args=None),
    )
    environment = Env(
        workflow,
        env_file=spec,
        envs_dir=str(prefix_dir),
    )
    environment._content = spec.read_bytes()
    environment._content_deploy = False
    environment._content_pin = False

    expected_hash = md5(
        str(prefix_dir.resolve()).encode() + spec.read_bytes()
    ).hexdigest()
    assert environment.hash == expected_hash
    assert Path(environment.address) == prefix_dir / f"{expected_hash}_"
