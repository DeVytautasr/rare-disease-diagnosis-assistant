#!/usr/bin/env bash
# Installs IGV 2.17.4 (the desktop app used by stage1_igv_assistant's
# igv_screenshot/evidence_panel tools) to ~/IGV_2.17.4 -- the first path
# run_igv_screenshot auto-detects (see bam_tools.py). Safe to re-run: if
# IGV is already installed there, this skips straight to the java/DISPLAY
# checks below.
#
# Usage:
#   bash scripts/install_igv.sh
#
# To install to (or use) a different location, set IGV_PATH to the full
# path of an igv.sh before running the MCP server / tools -- it's checked
# before ~/IGV_2.17.4, ~/igv, and /opt/igv. This script itself always
# targets ~/IGV_2.17.4.

set -euo pipefail

IGV_VERSION="2.17.4"
IGV_URL="https://data.broadinstitute.org/igv/projects/downloads/2.17/IGV_${IGV_VERSION}.zip"
INSTALL_DIR="${HOME}/IGV_${IGV_VERSION}"
IGV_SH="${INSTALL_DIR}/igv.sh"

if [ -x "$IGV_SH" ]; then
    echo "IGV ${IGV_VERSION} already installed at ${IGV_SH} -- skipping download."
else
    if ! command -v curl >/dev/null 2>&1; then
        echo "ERROR: curl is required to download IGV but was not found on PATH." >&2
        exit 1
    fi
    if ! command -v python3 >/dev/null 2>&1; then
        echo "ERROR: python3 is required to extract the IGV archive but was not found on PATH." >&2
        exit 1
    fi

    echo "Installing IGV ${IGV_VERSION} to ${INSTALL_DIR} ..."
    TMP_DIR="$(mktemp -d)"
    trap 'rm -rf "$TMP_DIR"' EXIT

    ZIP_PATH="${TMP_DIR}/IGV_${IGV_VERSION}.zip"
    echo "Downloading ${IGV_URL} ..."
    if ! curl -fSL --retry 3 -o "$ZIP_PATH" "$IGV_URL"; then
        echo "ERROR: failed to download IGV from ${IGV_URL}" >&2
        exit 1
    fi

    echo "Extracting ..."
    if ! python3 -c "
import zipfile
with zipfile.ZipFile('${ZIP_PATH}') as z:
    z.extractall('${TMP_DIR}')
"; then
        echo "ERROR: failed to extract ${ZIP_PATH} -- the download may be corrupt." >&2
        exit 1
    fi

    EXTRACTED_DIR="${TMP_DIR}/IGV_${IGV_VERSION}"
    if [ ! -d "$EXTRACTED_DIR" ]; then
        echo "ERROR: expected ${EXTRACTED_DIR} after extraction but it's not there --" >&2
        echo "IGV's zip layout may have changed since this script was written." >&2
        exit 1
    fi

    rm -rf "$INSTALL_DIR"
    mv "$EXTRACTED_DIR" "$INSTALL_DIR"
    # The zip stores igv.sh as executable, but Python's zipfile doesn't
    # apply that permission bit on extraction -- set it explicitly.
    chmod +x "$IGV_SH"
    echo "Installed to ${INSTALL_DIR}"
fi

if [ ! -x "$IGV_SH" ]; then
    echo "ERROR: ${IGV_SH} not found or not executable after install." >&2
    exit 1
fi

echo
echo "Checking runtime requirements ..."

if command -v java >/dev/null 2>&1; then
    echo "  java: OK ($(command -v java))"
else
    echo "  WARNING: java not found on PATH. IGV requires it. This project's"
    echo "  'rda' conda environment provides java -- run 'conda activate rda'"
    echo "  (or 'conda run -n rda ...') before using IGV-dependent tools."
fi

if [ -n "${DISPLAY:-}" ]; then
    echo "  DISPLAY: OK (${DISPLAY})"
else
    echo "  WARNING: DISPLAY is not set. IGV's screenshot tools need a usable"
    echo "  X11 display -- on WSL2, WSLg normally provides one automatically"
    echo "  (DISPLAY=:0). On a headless machine with no display at all, run"
    echo "  IGV-dependent tools under 'xvfb-run' (install Xvfb separately;"
    echo "  this script does not install it)."
fi

echo
echo "IGV ${IGV_VERSION} ready at: ${IGV_SH}"
echo "This is run_igv_screenshot's default auto-detected path, so no extra"
echo "configuration is needed. To use a different IGV install instead, set"
echo "IGV_PATH=/path/to/igv.sh before running the MCP server / tools."
