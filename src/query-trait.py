#!/usr/bin/env python3
# -*- coding: utf-8 -*-


### TRAIT QUERY

#### Tuesday, May 31, 2022, 1:20 PM
#Partly cloudy ⛅️  🌡️+76°F (feels +80°F, 69%) 🌬️→12mph 🌑 Tue May 31 07:24:28 2022
#W22Q2 – 151 ➡️ 213 – 20 ❇️ 345

#Version 0.3
#NYP Light rain 🌦   🌡️+45°F (feels +39°F, 76%) 🌬️↙16mph 🌓 Mon Mar 27 22:21:15 2023
#W13Q1 – 86 ➡️ 278 – 320 ❇️ 44

# GWAS summaries of studied traits


import os
import sys
import traceback

from config import (alfredError, alfredItems, cappedQuery, readOnly,
                    requireDatabase, truncationItem)

MYENTRY_Q = os.getenv('myENTRY_Q', '')  # this is breadcrumbs to enable the 'back' feature
MYSOURCE = os.getenv('mySource', '')
MYARG = sys.argv[1] if len(sys.argv) > 1 else ''


def queryTraits():

    if MYSOURCE == "traitGene" and MYARG == '':
        MYINPUT = MYENTRY_Q

    else:
        MYINPUT = MYARG

    MYINPUT = (MYINPUT or '').strip()
    MYQUERY = "%" + MYINPUT + "%"

    db = readOnly()
    rs, total = cappedQuery(
        db, "MAPPED_TRAIT, ImplicatedGenes_count, papers_count",
        "traitCounts", "MAPPED_TRAIT LIKE ?", [MYQUERY],
        "ORDER BY papers_count DESC, ImplicatedGenes_count DESC")
    db.close()

    if not rs:
        alfredError("No matches", "Try a different query")
        return

    items = []

    for countR, r in enumerate(rs, start=1):

        FeatureName = r['MAPPED_TRAIT']

        PapersCount = r['papers_count']
        paperString = "paper" if (PapersCount == 1) else "papers"

        GeneCount = r['ImplicatedGenes_count']
        geneString = "gene" if (GeneCount == 1) else "genes"

        subtitleString = (f"{countR}/{total:,}"
                          f" – {GeneCount} {geneString}, from {PapersCount} {paperString}")

        bigText = f"{FeatureName} – {GeneCount:,} {geneString}, from {PapersCount:,} {paperString}"
        titleString = f"{FeatureName}: {GeneCount:,} {geneString}, {PapersCount:,} {paperString}"

        #### COMPILING OUTPUT
        items.append({
            "title": titleString,
            "subtitle": subtitleString,

            "arg": FeatureName,

            "variables": {
                "currentTrait": FeatureName,
                "myENTRY_Q": MYINPUT,
                "currentTITLE": titleString,
                "bigText": bigText,
                "mySource": "GWT"
            },

            "icon": {
                "path": ""
            }
        })

    if total > len(rs):
        items.append(truncationItem(len(rs), total))

    alfredItems(items)


def main():
    requireDatabase()
    queryTraits()


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        traceback.print_exc(file=sys.stderr)
        alfredError(f"Error: {type(err).__name__}", str(err))
