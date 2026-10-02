"""Explicit configuration, with autodiscovery as the fallback.

Phase 11. Every path the tool needs used to be hardcoded to this developer's
home directory. Someone else's machine has a different layout, so each is now
resolvable in a fixed precedence:

    command-line flag  >  environment variable  >  config file  >  built-in default

and the built-in default is the old hardcoded value, so an existing install
behaves exactly as before if nothing is configured.

Nothing here raises. A missing config file, an unreadable one, or a path that
does not exist all resolve to a STATUS that the caller reports — the tool is
usable without any of these, just visibly less capable. That is the whole point:
the failure this module exists to prevent is the exclude template being absent
and the filter chain quietly running one step short.

Config file (stdlib configparser, no dependency):

    [paths]
    data_dir         = ~/public_data
    exclude_template = ~/reference/human.hg38.excl.tsv
    igv              = ~/IGV_2.17.4/igv.sh
    gene_table       = ~/public_data/annotation/genes_grch38.tsv.gz
    omim_mim2gene    = ~/public_data/annotation/mim2gene.txt
    gene_disorders   = ~/public_data/annotation/genes_to_disease.txt

    [datasets]
    DEMO = /srv/data/demo.bam

    [candidates]
    DEMO = /srv/data/demo.bcf

    [test_data]
    labels = DEMO

Labels registered explicitly (here or on the command line) are treated by the
interface as private data unless [test_data] lists them; autodiscovered files
under data_dir are test data.

Searched, first readable wins: $SV_CONFIG, ./sv-assistant.conf,
~/.config/sv-assistant/config.ini
"""
import configparser
import os

ENV = {"data_dir": "SV_DATA_DIR",
       "exclude_template": "SV_EXCLUDE_TEMPLATE",
       "igv": "IGV_PATH",
       "gene_table": "SV_GENE_TABLE",
       "omim_mim2gene": "SV_OMIM_MIM2GENE",
       "gene_disorders": "SV_GENE_DISORDERS"}

DEFAULTS = {"data_dir": "~/public_data",
            "exclude_template": "~/reference/human.hg38.excl.tsv",
            "igv": "",          # empty: bam_tools' own IGV search runs
            # Phase 26, all optional (see scripts/make_gene_table.py)
            "gene_table": "~/public_data/annotation/genes_grch38.tsv.gz",
            "omim_mim2gene": "~/public_data/annotation/mim2gene.txt",
            "gene_disorders": "~/public_data/annotation/genes_to_disease.txt"}

CONFIG_CANDIDATES = [os.environ.get("SV_CONFIG", ""),
                     "sv-assistant.conf",
                     os.path.expanduser("~/.config/sv-assistant/config.ini")]


def _read():
    cp = configparser.ConfigParser()
    cp.optionxform = str                      # preserve dataset LABEL case
    for c in CONFIG_CANDIDATES:
        if c and os.path.isfile(c):
            try:
                cp.read(c)
                return cp, os.path.abspath(c)
            except Exception as e:
                return configparser.ConfigParser(), f"{c} (unreadable: {e})"
    return cp, None


_CP, CONFIG_FILE = _read()


def _expand(p, relative_to=None):
    """Expand ~ and $VARS. A RELATIVE path is resolved against the directory the
    config file lives in, not the working directory — a demo bundle must work
    wherever it is unzipped, and `cd` must not change what the tool loads."""
    if not p:
        return p
    out = os.path.expanduser(os.path.expandvars(p))
    if not os.path.isabs(out) and relative_to:
        out = os.path.normpath(os.path.join(relative_to, out))
    return out


def _config_dir():
    return os.path.dirname(CONFIG_FILE) if CONFIG_FILE and os.path.isfile(CONFIG_FILE) else None


def get(key, cli=None):
    """Resolve one path. Returns (value, source)."""
    if cli:
        return _expand(cli), "command line"
    env = os.environ.get(ENV.get(key, ""), "")
    if env:
        return _expand(env), f"${ENV[key]}"
    if _CP.has_option("paths", key):
        return (_expand(_CP.get("paths", key), _config_dir()),
                f"config file {CONFIG_FILE}")
    return _expand(DEFAULTS.get(key, "")), "built-in default"


def status(key, cli=None, kind="file"):
    """Resolve and check. Never raises; the caller reports."""
    value, source = get(key, cli)
    if not value:
        return {"key": key, "path": None, "source": source, "found": False,
                "reason": "not configured"}
    ok = os.path.isdir(value) if kind == "dir" else os.path.isfile(value)
    return {"key": key, "path": value, "source": source, "found": bool(ok),
            "reason": None if ok else f"no such {kind}: {value}"}


def registered(section):
    """Explicit LABEL=PATH entries from the config file, expanded."""
    if not _CP.has_section(section):
        return {}
    return {k: _expand(v, _config_dir()) for k, v in _CP.items(section) if v}


def test_data_labels():
    """Labels the config file declares to be test data: [test_data] labels = A, B.
    The interface treats every other explicitly registered label as private, so a
    cloud model never receives it without the person's confirmation."""
    if not _CP.has_option("test_data", "labels"):
        return set()
    return {x.strip() for x in _CP.get("test_data", "labels").split(",") if x.strip()}
