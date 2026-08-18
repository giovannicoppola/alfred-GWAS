#!/usr/bin/env python3

### Gene QUERY

#### Tuesday, June 7, 2022, 5:37 PM
# Partly cloudy ⛅️  🌡️+71°F (feels +71°F, 63%) 🌬️↑17mph 🌓 Tue Jun  7 17:37:33 2022
# W23Q2 – 158 ➡️ 206 – 85 ❇️ 279

## revision for v0.3
#NYP – Clear ☀️   🌡️+55°F (feels +54°F, 26%) 🌬️→6mph 🌓 Mon Mar 27 00:55:03 2023
#W13Q1 – 86 ➡️ 278 – 320 ❇️ 44

# checking summaries from the GWAS catalog, by gene


import os
import sqlite3
import sys
import traceback

from config import INDEX_DB, alfredError, alfredItems, requireDatabase

MYSOURCE = os.getenv('mySource', '')
MYENTRY_Q = os.getenv('myENTRY_Q', '')
MYGENE = os.getenv('currentGeneID', '')
MYARG = sys.argv[1] if len(sys.argv) > 1 else ''


def queryGenes():

    if MYSOURCE == "GWG" and MYARG == '':
        MYINPUT = MYENTRY_Q

    elif MYSOURCE == "geneMasterSearch":
        MYINPUT = MYGENE

    else:
        MYINPUT = MYARG

    MYINPUT = (MYINPUT or '').strip()

    orderS = 'ORDER BY nTraits DESC'

    # flags are matched as whole words, so '--p' never fires on '--paper'
    words = MYINPUT.split()
    if "--p" in words:
        orderS = "ORDER BY nPapers DESC"
        MYINPUT = ' '.join(word for word in words if word != "--p")

    MYQUERY = "%" + MYINPUT + "%"

    db = sqlite3.connect(INDEX_DB)
    db.row_factory = sqlite3.Row

    # COALESCE keeps genes with no annotation searchable by their Ensembl id
    rs = db.execute(f"""SELECT *
        FROM geneCounts
        WHERE COALESCE(searchField, gene) LIKE ? {orderS}""", (MYQUERY,)).fetchall()
    db.close()

    if not rs:
        alfredError("No matches", "Try a different query")
        return

    items = []
    myResLen = len(rs)

    for countR, r in enumerate(rs, start=1):

        geneID = r['gene']
        FeatureName = r['GeneName'] or geneID

        PapersCount = r['nPapers']
        paperString = "paper" if (PapersCount == 1) else "papers"

        TraitCount = r['nTraits']
        traitString = "trait" if (TraitCount == 1) else "traits"

        subtitleString = (f"{countR}/{myResLen}"
                          f" – associated with {TraitCount:,} {traitString}, "
                          f"from {PapersCount:,} {paperString}")

        titleString = f"{FeatureName}: {TraitCount:,} {traitString} ({PapersCount:,} {paperString})"

        #### COMPILING OUTPUT
        items.append({
            "title": titleString,
            "subtitle": subtitleString,

            "arg": FeatureName,

            "variables": {
                "currentTrait": FeatureName,
                "currentGenes": geneID,
                "currentTITLE": titleString,
                "myENTRY_Q": MYINPUT,
                "mySource": "GWG"
            },

            "icon": {
                "path": ""
            }
        })

    alfredItems(items)


def main():
    requireDatabase()
    queryGenes()


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        traceback.print_exc(file=sys.stderr)
        alfredError(f"Error: {type(err).__name__}", str(err))
