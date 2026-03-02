#!/usr/bin/env python3


import os
import sys
import zipfile
import sqlite3



WF_DATA = os.getenv('alfred_workflow_data')
INDEX_DB = WF_DATA + '/index.db'

if not os.path.exists(WF_DATA):
    os.makedirs(WF_DATA)

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
    elif os.path.exists('index.db'):  # there is a new version, possibly rebuilt via script: replace the version in DATA
        os.rename('index.db', INDEX_DB)


def fetchColophon():
    
    # Importing the gene annotation table from the gene lookup DB
    conn = sqlite3.connect(INDEX_DB)
    cursor = conn.cursor()
    cursor.execute(f"SELECT * FROM colophon")
    rs = cursor.fetchone()
    conn.close()
    return rs [0]



checkDatabase()
colophon = fetchColophon()
GWAS_REF = f"ref: GWAS catalog, {colophon}"