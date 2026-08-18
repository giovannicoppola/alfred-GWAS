#!/usr/bin/env python3

"""Shared configuration and helpers for the alfred-GWAS scripts."""

import json
import os
import sqlite3
import sys
import zipfile


HERE = os.path.dirname(os.path.abspath(__file__))
WF_DATA = os.getenv('alfred_workflow_data') or HERE
INDEX_DB = os.path.join(WF_DATA, 'index.db')

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


def fetchColophon():
    """Return the catalog version string, or None if the database is unusable.

    Opened read-only so that a missing file is reported instead of being
    created as an empty database.
    """
    try:
        conn = sqlite3.connect(f"file:{INDEX_DB}?mode=ro", uri=True)
        try:
            rs = conn.execute("SELECT colophon FROM colophon").fetchone()
        finally:
            conn.close()
    except sqlite3.Error:
        return None

    return rs[0] if rs else None


def requireDatabase():
    """Tell the user to rebuild instead of crashing when the database is missing."""
    if COLOPHON is None:
        alfredError("GWAS database not available",
                    "Run ::rebuild to download and build the catalog")
        sys.exit(0)


checkDatabase()
COLOPHON = fetchColophon()
GWAS_REF = f"ref: GWAS catalog, {COLOPHON}" if COLOPHON else "ref: GWAS catalog"
