"""Prepare measured lossy wall-film n,k data from the user Excel workbook."""

from __future__ import annotations

import argparse
import os
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/private/tmp") / "matplotlib-codex-cache"))

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.io_utils import ensure_dir, project_path, save_figure


QUANT_LO = 8.1014
QUANT_HI = 12.962


def _paths() -> dict[str, Path]:
    return {
        "csv": project_path(
            "data", "processed", "measured_lossy_wall_film_nk_sheet2.csv",
        ),
        "report": project_path(
            "results", "diagnostics_v2", "reports",
            "V02_measured_wall_film_nk_report.md",
        ),
        "figure": project_path(
            "results", "diagnostics_v2", "figures",
            "V02_measured_wall_film_nk.png",
        ),
    }


def _write_failure_report(path: Path, exc: Exception, args: argparse.Namespace) -> None:
    ensure_dir(path.parent)
    path.write_text(
        "# V02 Measured Wall-Film n,k Preparation\n\n"
        "Overall result level: **FAIL**\n\n"
        f"Input workbook: `{args.input}`\n\n"
        f"Sheet: `{args.sheet}`\n\n"
        f"Failure: `{type(exc).__name__}: {exc}`\n\n"
        "The measured_lossy_wall_film material cannot be used until the Excel "
        "file is present and the first three Sheet2 columns contain numeric "
        "wavelength_um, n, k data starting on row 3.\n",
        encoding="utf-8",
    )


def load_and_clean_excel(input_path: Path, sheet_name: str) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    try:
        raw = pd.read_excel(input_path, sheet_name=sheet_name, header=None, skiprows=2)
    except ImportError:
        raw = _read_xlsx_sheet_without_openpyxl(input_path, sheet_name, skiprows=2)
    if raw.shape[1] < 3:
        raise ValueError("Sheet must contain at least three columns.")
    df = raw.iloc[:, :3].copy()
    df.columns = ["wavelength_um", "n", "k"]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(how="any").copy()
    df = df.sort_values("wavelength_um").reset_index(drop=True)
    if df.empty:
        raise ValueError("No numeric n,k rows found after cleaning.")
    if not (df["wavelength_um"] > 0).all():
        raise ValueError("All wavelength_um values must be positive.")
    if not (df["n"] > 0).all():
        raise ValueError("All n values must be positive.")
    if not (df["k"] >= 0).all():
        raise ValueError("All k values must be non-negative.")
    if df.duplicated("wavelength_um").any():
        dupes = df.loc[df.duplicated("wavelength_um"), "wavelength_um"].tolist()
        raise ValueError(f"Duplicate wavelength_um values found: {dupes[:5]}")
    if df.isna().any().any():
        raise ValueError("Missing values remain after cleaning.")
    return df


def _xlsx_col_to_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    out = 0
    for ch in letters:
        out = out * 26 + (ord(ch.upper()) - ord("A") + 1)
    return out - 1


def _read_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    try:
        data = zf.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(data)
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    strings = []
    for si in root.findall("x:si", ns):
        texts = [node.text or "" for node in si.findall(".//x:t", ns)]
        strings.append("".join(texts))
    return strings


def _sheet_xml_path(zf: zipfile.ZipFile, sheet_name: str) -> str:
    ns = {
        "x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
        "rel": "http://schemas.openxmlformats.org/package/2006/relationships",
    }
    workbook = ET.fromstring(zf.read("xl/workbook.xml"))
    target_rid = None
    available = []
    for sheet in workbook.findall(".//x:sheet", ns):
        name = sheet.attrib.get("name")
        available.append(name)
        if name == sheet_name:
            target_rid = sheet.attrib.get(f"{{{ns['r']}}}id")
            break
    if target_rid is None:
        raise ValueError(f"Sheet {sheet_name!r} not found. Available sheets: {available}")
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    for rel in rels.findall("rel:Relationship", ns):
        if rel.attrib.get("Id") == target_rid:
            target = rel.attrib["Target"]
            if target.startswith("/"):
                return target.lstrip("/")
            return "xl/" + target.lstrip("/")
    raise ValueError(f"No worksheet relationship found for {sheet_name!r}")


def _cell_value(cell, shared_strings: list[str]) -> object:
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        texts = [node.text or "" for node in cell.findall(".//x:t", ns)]
        return "".join(texts)
    value = cell.find("x:v", ns)
    if value is None:
        return None
    text = value.text
    if cell_type == "s":
        return shared_strings[int(text)]
    return text


def _read_xlsx_sheet_without_openpyxl(
    input_path: Path,
    sheet_name: str,
    skiprows: int,
) -> pd.DataFrame:
    """Minimal xlsx reader for first-column numeric tables when openpyxl is absent."""
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    rows_out = []
    with zipfile.ZipFile(input_path) as zf:
        shared = _read_shared_strings(zf)
        sheet_path = _sheet_xml_path(zf, sheet_name)
        root = ET.fromstring(zf.read(sheet_path))
        for row in root.findall(".//x:sheetData/x:row", ns):
            row_number = int(row.attrib.get("r", "0"))
            if row_number <= skiprows:
                continue
            values = [None, None, None]
            for cell in row.findall("x:c", ns):
                ref = cell.attrib.get("r", "")
                col = _xlsx_col_to_index(ref)
                if 0 <= col < 3:
                    values[col] = _cell_value(cell, shared)
            rows_out.append(values)
    return pd.DataFrame(rows_out)


def plot_nk(df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.2), sharex=True)
    for ax, col, color in [
        (axes[0], "n", "#4C78A8"),
        (axes[1], "k", "#F58518"),
    ]:
        ax.plot(df["wavelength_um"], df[col], color=color, lw=1.6)
        ax.axvspan(QUANT_LO, QUANT_HI, color="#E8EEF7", alpha=0.65)
        ax.set_ylabel(col)
        ax.grid(alpha=0.25)
    axes[1].set_xlabel("Wavelength (um)")
    fig.suptitle("Measured lossy wall-film n,k from Sheet2")
    save_figure(fig, out_path)
    plt.close(fig)


def write_report(df: pd.DataFrame, path: Path, args: argparse.Namespace) -> None:
    ensure_dir(path.parent)
    band = df[(df["wavelength_um"] >= QUANT_LO) & (df["wavelength_um"] <= QUANT_HI)]
    allowed = (
        len(band) >= 3
        and (df["wavelength_um"] > 0).all()
        and (df["n"] > 0).all()
        and (df["k"] >= 0).all()
        and not df.duplicated("wavelength_um").any()
        and not df.isna().any().any()
    )
    path.write_text(
        "# V02 Measured Wall-Film n,k Preparation\n\n"
        f"Overall result level: **{'CODE_PASS' if allowed else 'FAIL'}**\n\n"
        "## Inputs\n\n"
        f"- Workbook: `{args.input}`\n"
        f"- Sheet: `{args.sheet}`\n"
        "- Row handling: first two rows skipped; first three columns renamed to "
        "`wavelength_um`, `n`, `k`.\n\n"
        "## Validation Summary\n\n"
        f"- Total numeric data points: {len(df)}\n"
        f"- Data wavelength range: {df['wavelength_um'].min():.6g}-"
        f"{df['wavelength_um'].max():.6g} um\n"
        f"- Points in {QUANT_LO:.4f}-{QUANT_HI:.3f} um: {len(band)}\n"
        f"- n range: {df['n'].min():.6g}-{df['n'].max():.6g}\n"
        f"- k range: {df['k'].min():.6g}-{df['k'].max():.6g}\n"
        f"- Negative k present: {bool((df['k'] < 0).any())}\n"
        f"- Allowed for this narrowband quantitative wall-film simulation: "
        f"{'yes' if allowed else 'no'}\n\n"
        "## Notes\n\n"
        "- The material is named `measured_lossy_wall_film`; no chemical "
        "identity is inferred from the Excel data.\n"
        "- No extrapolation is allowed outside the measured wavelength range.\n",
        encoding="utf-8",
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare measured wall-film n,k Excel data.")
    parser.add_argument(
        "--input",
        type=Path,
        default=project_path("data", "raw", "measured_lossy_wall_film_nk.xlsx"),
    )
    parser.add_argument("--sheet", default="Sheet2")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    paths = _paths()
    for path in paths.values():
        ensure_dir(path.parent)
    try:
        df = load_and_clean_excel(args.input, args.sheet)
        df.to_csv(paths["csv"], index=False)
        plot_nk(df, paths["figure"])
        write_report(df, paths["report"], args)
    except Exception as exc:
        _write_failure_report(paths["report"], exc, args)
        raise


if __name__ == "__main__":
    main()
