#!/usr/bin/env python3


import os
import sys
import zipfile
import sqlite3
import tempfile



WF_DATA = os.getenv('alfred_workflow_data') or os.path.join(tempfile.gettempdir(), 'alfred-gwas')
INDEX_DB = os.path.join(WF_DATA, 'index.db')

os.makedirs(WF_DATA, exist_ok=True)

def log(s, *args):
    if args:
        s = s % args
    print(s, file=sys.stderr)

    
def logF(log_message, file_name):
    with open(file_name, "a") as f:
        f.write(log_message + "\n")

def checkDatabase():

    DB_ZIPPED = 'index.db.zip'

    if os.path.exists(DB_ZIPPED):  # there is a zipped database: distribution version
        log ("found distribution database, extracting")
        with zipfile.ZipFile(DB_ZIPPED, "r") as zip_ref:
            zip_ref.extractall(WF_DATA)
        os.remove (DB_ZIPPED)


def databaseReady():
    """True once ::rebuild has built the index. Checks for a table, not just
    the file: an earlier connect can leave an empty index.db behind."""
    if not os.path.exists(INDEX_DB) or os.path.getsize(INDEX_DB) == 0:
        return False
    try:
        conn = sqlite3.connect(INDEX_DB)
        found = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' "
                             "AND name='GeneTrait'").fetchone()
        conn.close()
        return found is not None
    except sqlite3.Error:
        return False


# shown by the gene and trait searches until the database has been built
NOT_BUILT = {"items": [{
    "title": "GWAS database not built yet",
    "subtitle": "Run the ::rebuild keyword to download the GWAS Catalog and build it",
    "valid": False,
}]}


def fetchColophon():

    if not os.path.exists(INDEX_DB):  # don't create an empty index.db
        return "unknown"
    # Importing the gene annotation table from the gene lookup DB
    conn = sqlite3.connect(INDEX_DB)
    cursor = conn.cursor()
    cursor.execute(f"SELECT * FROM colophon")
    rs = cursor.fetchone()
    conn.close()
    if not rs:
        return "unknown"
    return rs[0]



try:
    checkDatabase()
    colophon = fetchColophon()
except Exception:
    colophon = "unknown"
GWAS_REF = f"ref: GWAS catalog, {colophon}"