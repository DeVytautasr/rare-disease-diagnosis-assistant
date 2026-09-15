#!/usr/bin/env bash
# Install the breakpoint-evidence assistant (MINIMAL tier).
#
# What this installs: Python packages only, into a virtual environment in this
# directory. It does NOT need, and does not install:
#   - a reference genome (no FASTA is ever opened)
#   - samtools or bcftools (pysam does that work in-process)
#   - delly (it produces candidate sets; this tool reads them)
#   - a GPU, or internet access at run time
#
# Optional extras, each independent:
#   IGV       -> picture panels        (install IGV separately; set IGV_PATH)
#   ollama    -> local model chat      (install ollama separately)
#   API key   -> Claude chat backend   (pip install -r requirements-api.txt)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
PY="${PYTHON:-python3}"
VENV="${VENV:-.venv}"

echo "==> checking Python"
if ! command -v "$PY" >/dev/null 2>&1; then
  echo "ERROR: $PY not found. Install Python 3.10 or newer, or set PYTHON=/path/to/python3." >&2
  exit 1
fi
"$PY" - <<'PYCHK'
import sys
if sys.version_info < (3, 10):
    sys.exit(f"ERROR: Python 3.10+ required (fastmcp and requests need it); found {sys.version.split()[0]}")
print(f"    Python {sys.version.split()[0]} OK")
PYCHK

echo "==> creating virtual environment in $VENV"
# A previous failed run can leave a half-built venv behind. Testing for the
# directory treated that as usable and then failed on a missing activate
# script; test for the thing we actually need.
if [ ! -f "$VENV/bin/activate" ]; then
  rm -rf "$VENV"
  if ! "$PY" -m venv "$VENV" 2>/tmp/venv_err.$$; then
    PYVER="$("$PY" -c 'import sys;print("%d.%d"%sys.version_info[:2])' 2>/dev/null || echo 3)"
    # Debian splits ensurepip out of the stdlib, so `python3 -m venv` fails on a
    # stock image. Telling the user to apt-install it is useless if they have no
    # root -- which is the common case on a managed laptop. A venv WITHOUT pip
    # still builds from the stdlib alone, and pip can be bootstrapped into it
    # over the network. Found by installing on a machine without sudo, not by
    # reading the docs.
    echo "    ensurepip unavailable; building the environment without pip and bootstrapping it"
    rm -rf "$VENV"
    if "$PY" -m venv --without-pip "$VENV" 2>>/tmp/venv_err.$$ \
       && "$VENV/bin/python" - <<'BOOT' 2>>/tmp/venv_err.$$
import ssl, sys, urllib.request, tempfile, os, runpy
url = "https://bootstrap.pypa.io/get-pip.py"
try:
    data = urllib.request.urlopen(url, timeout=120, context=ssl.create_default_context()).read()
except Exception as e:
    sys.exit(f"could not download {url}: {e}")
fd, path = tempfile.mkstemp(suffix=".py")
os.write(fd, data); os.close(fd)
sys.argv = [path, "--quiet"]
runpy.run_path(path, run_name="__main__")
BOOT
    then
      rm -f /tmp/venv_err.$$
      echo "    pip bootstrapped OK"
    else
      cat /tmp/venv_err.$$ >&2 || true
      rm -f /tmp/venv_err.$$
      rm -rf "$VENV"
      echo >&2
      echo "ERROR: could not create a virtual environment." >&2
      echo "  Either install the venv package (needs root):" >&2
      echo "       sudo apt install python${PYVER}-venv" >&2
      echo "       (if that package does not exist: sudo apt install python3-venv)" >&2
      echo "  or run this script on a machine with internet access, which lets it" >&2
      echo "  bootstrap pip by itself without root." >&2
      exit 1
    fi
  else
    rm -f /tmp/venv_err.$$
  fi
fi
if [ ! -f "$VENV/bin/activate" ]; then
  echo "ERROR: $VENV exists but has no bin/activate — the environment is incomplete." >&2
  echo "  Remove it and run again:   rm -rf $VENV && bash install.sh" >&2
  exit 1
fi
# shellcheck disable=SC1091
. "$VENV/bin/activate"

echo "==> upgrading pip"
python -m pip install --quiet --upgrade pip

echo "==> installing MINIMAL requirements"
python -m pip install -r requirements.txt

if [ "${WITH_API:-0}" = "1" ]; then
  echo "==> installing the optional Anthropic API backend"
  python -m pip install -r requirements-api.txt
fi

echo "==> verifying"
if python -m stage1_igv_assistant.ui --check; then
  echo
  echo "Install OK."
  echo "  Start it with:   . $VENV/bin/activate && python -m stage1_igv_assistant.ui"
  echo "  Then open:       http://127.0.0.1:8765"
  echo
  echo "  If the dropdowns are empty, tell it where your data is:"
  echo "    cp sv-assistant.conf.example sv-assistant.conf   # then edit data_dir"
  echo "  or point it at files directly:"
  echo "    python -m stage1_igv_assistant.ui \\"
  echo "        --dataset SAMPLE=/path/to/sample.bam \\"
  echo "        --candidates SAMPLE=/path/to/sample.bcf"
else
  echo "Install FAILED the check above. The banner says what is missing." >&2
  exit 1
fi
