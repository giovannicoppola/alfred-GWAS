#!/usr/bin/env python3

"""Shared configuration and helpers for the alfred-GWAS scripts."""

import json
import os
import sqlite3
import sys
import zipfile
from urllib.parse import quote


HERE = os.path.dirname(os.path.abspath(__file__))
WF_DATA = os.getenv('alfred_workflow_data') or HERE
INDEX_DB = os.path.join(WF_DATA, 'index.db')

# Alfred shows a fixed-height list, so serialising thousands of items only costs
# time. Broad searches are capped and the user is told to narrow them.
MAX_RESULTS = 200

os.makedirs(WF_DATA, exist_ok=True)


def log(s, *args):
    if args:
        s = s % args
    print(s, file=sys.stderr)


def logF(log_message, file_name):
    with open(file_name, "a") as f:
        f.write(log_message + "\n")


def alfredItems(items, variables=None):
    """Print a feedback payload for Alfred."""
    result = {"items": items, "variables": variables or {}}
    print(json.dumps(result))


def alfredError(title, subtitle):
    """Print a single non-actionable error item for Alfred."""
    alfredItems([{
        "title": title,
        "subtitle": subtitle,
        "valid": False,
        "icon": {"path": "icons/Warning.png"}
    }])


def checkDatabase():

    DB_ZIPPED = os.path.join(HERE, 'index.db.zip')

    if os.path.exists(DB_ZIPPED):  # there is a zipped database: distribution version
        log("found distribution database, extracting")
        with zipfile.ZipFile(DB_ZIPPED, "r") as zip_ref:
            zip_ref.extractall(WF_DATA)
        os.remove(DB_ZIPPED)


def readOnly():
    """Open the index read-only, so a missing file is reported rather than created.

    The path is percent-encoded: it is read as a URI, and an unescaped '?' or
    '#' in the user's home directory would otherwise truncate it.
    """
    db = sqlite3.connect(f"file:{quote(INDEX_DB)}?mode=ro", uri=True)
    db.row_factory = sqlite3.Row
    return db


def fetchColophon():
    """Return the catalog version string, or None if the database is unusable."""
    try:
        conn = readOnly()
        try:
            rs = conn.execute("SELECT colophon FROM colophon").fetchone()
        finally:
            conn.close()
    except sqlite3.Error:
        return None

    return rs[0] if rs else None


def cappedQuery(db, columns, table, where, params, orderBy):
    """Fetch at most MAX_RESULTS rows, with the true total when more matched.

    One extra row reveals whether the result was cut short, so the second
    COUNT(*) is only paid on the broad searches that actually need it.
    """
    rows = db.execute(
        f"SELECT {columns} FROM {table} WHERE {where} {orderBy} LIMIT ?",
        (*params, MAX_RESULTS + 1)).fetchall()

    if len(rows) <= MAX_RESULTS:
        return rows, len(rows)

    total = db.execute(f"SELECT COUNT(*) FROM {table} WHERE {where}", params).fetchone()[0]
    return rows[:MAX_RESULTS], total


def truncationItem(shown, total):
    """A closing item telling the user the list was cut short."""
    return {
        "title": f"… and {total - shown:,} more",
        "subtitle": f"showing the first {shown:,} of {total:,} – type more to narrow the search",
        "valid": False,
        "icon": {"path": "icon.png"}
    }


def requireDatabase():
    """Tell the user to rebuild instead of crashing when the database is missing."""
    if COLOPHON is None:
        alfredError("GWAS database not available",
                    "Run ::rebuild to download and build the catalog")
        sys.exit(0)


checkDatabase()
COLOPHON = fetchColophon()
GWAS_REF = f"ref: GWAS catalog, {COLOPHON}" if COLOPHON else "ref: GWAS catalog"
