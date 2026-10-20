from pathlib import Path
import re


REPORT_FILE = Path("structural_reports_all.txt")


def test_combined_report_has_correct_cao_mole_fraction_labels():
    text = REPORT_FILE.read_text(encoding="utf-8")

    labels = re.findall(
        r"^=+\s*x\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*=+$",
        text,
        flags=re.MULTILINE,
    )

    assert labels == ["0.10", "0.20", "0.30", "0.40"]
