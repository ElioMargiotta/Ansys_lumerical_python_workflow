"""
Loop over the wavelengths in WAVELENGTHS, patching `params.yaml`
(λmin = λmax = w), launching app.py, then restoring the file.

Requires:
    pip install "ruamel.yaml>=0.17"            # round-trip YAML
"""

from pathlib import Path
import subprocess, sys, os, shutil, tempfile
from ruamel.yaml import YAML                   # ⬅ round-trip loader

# ──────────────────────────────────── constants ──────────────────────────────
PARAM_FILE   = Path("configs/params.yaml")             # must exist
WAVELENGTHS  = [0.50e-6, 0.55e-6, 0.60e-6]     # metres  (edit as you like)
PY           = sys.executable                  # launch with *same* interpreter
ENV_UTF8     = {**os.environ, "PYTHONUTF8": "1"}  # force UTF-8 inside the child
YAML         = YAML()                          # one YAML() instance is enough
YAML.preserve_quotes = True                    # keep "1.5" etc. quoted

# ───────────────────────────────── helper ------------------------------------
def run_app() -> None:
    """Call app.py and fail loud if it exits non-zero."""
    completed = subprocess.run(
        [PY, "app.py"],
        capture_output=True,
        text=True,
        encoding="utf-8",          # decode what we captured
        env=ENV_UTF8,              # child prints UTF-8 → no 🧪 crash
        check=True                 # raise on error
    )
    # Optional: stream child’s stdout
    if completed.stdout:
        print(completed.stdout, end="")

# ─────────────────────────────── main sweep ----------------------------------
if __name__ == "__main__":
    original = PARAM_FILE.read_text(encoding="utf-8")   # exact bytes → full restore
    print(f"🗂  Loaded {PARAM_FILE}")

    try:
        for w in WAVELENGTHS:
            print(f"\n🔧  Setting λ = {w:.2e} m in {PARAM_FILE} …")
            # ‥parse original every iteration, so comments stay *untouched*
            data = YAML.load(original)
            data["Source"]["lambda_min"] = w
            data["Source"]["lambda_max"] = w

            with PARAM_FILE.open("w", encoding="utf-8") as fh:   # ← give YAML a stream
                YAML.dump(data, fh)
            print("🚀  Running simulation …")
            run_app()
            print("✅  Run finished.")

    finally:
        PARAM_FILE.write_text(original, encoding="utf-8")
        print(f"\n♻️  Restored pristine {PARAM_FILE}")
