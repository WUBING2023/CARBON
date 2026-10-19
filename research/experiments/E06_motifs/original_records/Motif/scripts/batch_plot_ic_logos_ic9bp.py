#!/usr/bin/env python3
"""
对每个 outputs/<TASK>/tfmodisco（已有 modisco_results.h5）在 **mode A** 下
绘制与 ``ic_scaled_logos_mode_A`` 相同原理的 IC logo，但在全长 PPM 上截取 **9 bp**：
选取 **per-position IC（相对均匀背景）之和最大** 的连续窗口。

默认：每个任务 **仅 1 张**（优先 ``pos_patterns/pattern_0``，否则首个 pos pattern，再否则首个 neg pattern），
输出到 ``tfmodisco/ic_scaled_logos_mode_A_ic9bp/``。

  cd 可解释性分析/Motif
  python3 scripts/batch_plot_ic_logos_ic9bp.py --jobs 8
  python3 scripts/batch_plot_ic_logos_ic9bp.py --mirror-mode-a --jobs 8   # 与已有 *.ic_fwd.png 一一对应（张数更多）

说明：222 个任务中若有 MoDISco 未检出任何 pos/neg pattern（pos 为空），则该任务跳过。
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Literal

import h5py

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from motif_analysis.modisco_h5_patterns import load_ppm_for_ic_logo
from motif_analysis.ppm_ic_trim import crop_ppm_max_ic_sum_window


def _load_ic_plot_module():
    path = ROOT / "scripts" / "plot_modisco_ic_logos.py"
    spec = importlib.util.spec_from_file_location("_plot_modisco_ic_logos", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_ic_plot = _load_ic_plot_module()
plot_ic_logo_png = _ic_plot.plot_ic_logo_png
group_path_to_filename_stem = _ic_plot.group_path_to_filename_stem

SKIP_DIR_NAMES = frozenset({"batch_motif_logs", "fixtures", "demo", ".git"})
IC_WIN = 9
MODE: Literal["A"] = "A"


def _discover(outputs_root: Path) -> list[Path]:
    out: list[Path] = []
    for child in sorted(outputs_root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIR_NAMES or child.name.startswith("."):
            continue
        out.append(child)
    return out


def stem_to_h5_group(stem: str) -> str:
    """pos_patterns.pattern_0.subpattern_2 -> pos_patterns/pattern_0/subpattern_2"""
    parts = stem.split(".")
    if len(parts) < 2:
        raise ValueError(f"bad stem: {stem!r}")
    return f"{parts[0]}/" + "/".join(parts[1:])


def _pick_primary_group(h5_path: Path) -> str | None:
    with h5py.File(h5_path, "r") as hf:
        if "pos_patterns/pattern_0" in hf:
            return "pos_patterns/pattern_0"
        if "pos_patterns" in hf:
            keys = sorted(k for k in hf["pos_patterns"].keys() if k.startswith("pattern_"))
            if keys:
                return f"pos_patterns/{keys[0]}"
        if "neg_patterns" in hf:
            keys = sorted(k for k in hf["neg_patterns"].keys() if k.startswith("pattern_"))
            if keys:
                return f"neg_patterns/{keys[0]}"
    return None


def _plot_one_ic9(
    h5_path: Path,
    group_path: str,
    out_png: Path,
) -> tuple[bool, str]:
    try:
        ppm = load_ppm_for_ic_logo(h5_path, group_path, MODE)
    except (KeyError, ValueError, TypeError) as e:
        return False, str(e)
    cropped, _ = crop_ppm_max_ic_sum_window(ppm, IC_WIN)
    plot_ic_logo_png(out_png, cropped, title=None)
    return True, ""


def _run_task_primary(payload: dict) -> dict:
    task_dir = Path(payload["task_dir"])
    h5 = task_dir / "tfmodisco" / "modisco_results.h5"
    out_dir = task_dir / "tfmodisco" / "ic_scaled_logos_mode_A_ic9bp"
    t0 = datetime.now().isoformat(timespec="seconds")
    if not h5.is_file():
        return {"task": task_dir.name, "exit": "no_h5", "stderr": "", "time": t0}
    gpath = _pick_primary_group(h5)
    if gpath is None:
        return {"task": task_dir.name, "exit": "no_pattern", "stderr": "empty pos/neg", "time": t0}
    stem = group_path_to_filename_stem(gpath)
    out_png = out_dir / f"{stem}.ic_fwd.png"
    if payload.get("skip_existing") and out_png.is_file():
        return {"task": task_dir.name, "exit": "skipped", "stderr": "", "time": t0}
    ok, err = _plot_one_ic9(h5, gpath, out_png)
    if not ok:
        return {"task": task_dir.name, "exit": "fail", "stderr": err, "time": t0}
    return {"task": task_dir.name, "exit": "0", "stderr": str(out_png), "time": t0}


def _run_task_mirror(payload: dict) -> dict:
    task_dir = Path(payload["task_dir"])
    h5 = task_dir / "tfmodisco" / "modisco_results.h5"
    mode_dir = task_dir / "tfmodisco" / "ic_scaled_logos_mode_A"
    out_dir = task_dir / "tfmodisco" / "ic_scaled_logos_mode_A_ic9bp"
    t0 = datetime.now().isoformat(timespec="seconds")
    if not h5.is_file() or not mode_dir.is_dir():
        return {"task": task_dir.name, "exit": "skip", "stderr": "no h5 or mode_A dir", "time": t0}
    stems: list[str] = []
    for p in sorted(mode_dir.glob("*.ic_fwd.png")):
        s = p.stem
        if s.endswith(".ic_fwd"):
            stems.append(s[: -len(".ic_fwd")])
    if not stems:
        return {"task": task_dir.name, "exit": "skip", "stderr": "no ic_fwd in mode_A", "time": t0}
    errs: list[str] = []
    n_ok = 0
    for stem in stems:
        # stem is like pos_patterns.pattern_0
        out_png = out_dir / f"{stem}.ic_fwd.png"
        if payload.get("skip_existing") and out_png.is_file():
            n_ok += 1
            continue
        try:
            gpath = stem_to_h5_group(stem)
        except ValueError as e:
            errs.append(f"{stem}: {e}")
            continue
        ok, err = _plot_one_ic9(h5, gpath, out_png)
        if ok:
            n_ok += 1
        else:
            errs.append(f"{stem}: {err}")
    if errs:
        return {
            "task": task_dir.name,
            "exit": "partial" if n_ok else "fail",
            "stderr": "; ".join(errs[:5]),
            "time": t0,
        }
    return {"task": task_dir.name, "exit": "0", "stderr": f"pngs>={n_ok}", "time": t0}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument(
        "--mirror-mode-a",
        action="store_true",
        help="与 ic_scaled_logos_mode_A 下每张 *.ic_fwd.png 一一生成 ic9bp 对照（而非每任务一张）",
    )
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        default=True,
        help="若目标 PNG 已存在则跳过（默认真）",
    )
    ap.add_argument("--no-skip-existing", action="store_true")
    ap.add_argument("--summary", action="store_true")
    args = ap.parse_args()
    out_root = args.outputs_root.resolve()
    skip_existing = args.skip_existing and not args.no_skip_existing

    queued: list[Path] = []
    for td in _discover(out_root):
        h5 = td / "tfmodisco" / "modisco_results.h5"
        if h5.is_file():
            queued.append(td)

    print(
        f"[ic9bp] outputs={out_root} tasks_with_h5={len(queued)} "
        f"mirror_mode_a={args.mirror_mode_a} ic_window={IC_WIN} jobs={args.jobs}",
        flush=True,
    )

    if args.summary:
        for td in queued[:40]:
            print(f"  {td.name}")
        if len(queued) > 40:
            print(f"  ... +{len(queued) - 40}")
        return

    worker = _run_task_mirror if args.mirror_mode_a else _run_task_primary
    payload_base: dict = {"skip_existing": skip_existing}
    jobs = min(max(1, args.jobs), len(queued) or 1)
    rows: list[dict[str, str]] = []

    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(worker, {**payload_base, "task_dir": str(td)}): td.name for td in queued}
        for fut in as_completed(futs):
            name = futs[fut]
            try:
                row = fut.result()
            except Exception as e:
                row = {"task": name, "exit": "exc", "stderr": repr(e), "time": ""}
            rows.append(row)
            exv = row["exit"]
            tag = "OK" if exv == "0" else ("SKIP" if exv in ("skipped", "skip") else f"FAIL({exv})")
            print(f"[ic9bp] {tag} {row['task']}", flush=True)
            if exv not in ("0", "skipped", "skip") and row.get("stderr"):
                for ln in str(row["stderr"]).splitlines()[-2:]:
                    if ln.strip():
                        print(f"    {ln[:500]}", flush=True)

    n_ok = sum(1 for r in rows if r["exit"] == "0")
    n_skip_pat = sum(1 for r in rows if r["exit"] == "no_pattern")
    n_skipped = sum(1 for r in rows if r["exit"] == "skipped")
    print(
        f"[ic9bp] done ok={n_ok} no_pattern={n_skip_pat} skipped_existing={n_skipped} "
        f"total_tasks={len(rows)}",
        flush=True,
    )

    if args.mirror_mode_a:
        n_out = sum(
            1
            for _ in out_root.glob("**/tfmodisco/ic_scaled_logos_mode_A_ic9bp/*.ic_fwd.png")
            if "fixtures" not in str(_.parent)
        )
        print(f"[ic9bp] png count under ic_scaled_logos_mode_A_ic9bp: {n_out}", flush=True)


if __name__ == "__main__":
    main()
