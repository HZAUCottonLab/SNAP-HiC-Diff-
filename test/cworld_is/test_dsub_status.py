from pathlib import Path
import os
import subprocess


REPO = Path(__file__).resolve().parents[2]
STATUS = REPO / "config" / "dsub_status.sh"


def run_status(tmp_path, state):
    fake_djob = tmp_path / "djob"
    fake_djob.write_text(
        "#!/usr/bin/env bash\n"
        "printf 'ID NAME STATE USER\\n'\n"
        f"printf '123 rule {state} user\\n'\n"
    )
    fake_djob.chmod(0o755)
    env = os.environ | {"PATH": f"{tmp_path}:{os.environ['PATH']}"}
    return subprocess.run(["bash", str(STATUS), "123"], text=True, capture_output=True, env=env)


def test_dsub_status_maps_scheduler_states_to_cluster_generic_protocol(tmp_path):
    assert run_status(tmp_path, "PENDING").stdout.strip() == "running"
    assert run_status(tmp_path, "RUNNING").stdout.strip() == "running"
    assert run_status(tmp_path, "SUCCEEDED").stdout.strip() == "success"
    assert run_status(tmp_path, "FAILED").stdout.strip() == "failed"
