"""Construct CARBON and verify its output contract without research data."""

from pathlib import Path
import sys

import torch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from CARBON.models.carbon_model import CARBON  # noqa: E402


def main() -> None:
    torch.manual_seed(0)
    model = CARBON(in_channels=5)
    model.eval()
    with torch.no_grad():
        output = model(torch.zeros((2, 101, 5), dtype=torch.float32))

    profile_shape = tuple(output["logits_profile"].shape)
    count_shape = tuple(output["log_total_count"].shape)
    if profile_shape != (2, 101) or count_shape != (2, 1):
        raise RuntimeError(
            f"Unexpected output shapes: logits_profile={profile_shape}, "
            f"log_total_count={count_shape}"
        )

    print(
        "CARBON smoke test passed: "
        f"logits_profile={profile_shape}, log_total_count={count_shape}"
    )


if __name__ == "__main__":
    main()

