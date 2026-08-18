#!/usr/bin/env python3

### ITEMS-QUERY
# showing individual items (i.e. papers, genes, loci)

#### Tuesday, May 31, 2022, 5:05 PM
# Partly cloudy ⛅️  🌡️+76°F (feels +80°F, 69%) 🌬️→12mph 🌑 Tue May 31 07:24:28 2022
# W22Q2 – 151 ➡️ 213 – 20 ❇️ 345

import os
import sqlite3
import sys
import traceback

from config import INDEX_DB, GWAS_REF, alfredError, alfredItems, requireDatabase


MYSOURCE = os.getenv('mySource', '')
MYENTRY = sys.argv[1] if len(sys.argv) > 1 else ''
MYTITLE = os.getenv('currentTITLE', '')
MYENTRY_Q = os.getenv('myENTRY_Q', '')

# Search modifiers, recognised as whole words and never passed to the search itself
FLAGS = ('--es', '--p')

# Map incoming source to the source needed to reproduce the previous step (for back from papers)
BACK_SOURCE_MAP = {
    "traitGene": "GWT",   # trait→gene→papers: back to gene list (needs GWT)
    "geneTrait": "GWG",   # gene→trait→papers: back to trait list (needs GWG)
    "genePap": "GWG",     # gene→papers: back to trait list (needs GWG)
}


def connect():
    db = sqlite3.connect(INDEX_DB)
    db.row_factory = sqlite3.Row
    return db


def hasFlag(text, flag):
    """True when flag appears as its own word (so '--p' never matches '--paper')."""
    return flag in text.split()


def stripFlags(text):
    return ' '.join(word for word in text.split() if word not in FLAGS)


def fmt(value, spec='.2f'):
    """Format a possibly NULL/NaN numeric column, returning '' when there is no value."""
    if value is None:
        return ''
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ''
    if number != number:  # NaN
        return ''
    return format(number, spec)


def orBetaBlock(low_value, high_value):
    """Render the OR/beta range, tolerating either end being absent."""
    low, high = fmt(low_value), fmt(high_value)
    if not low and not high:
        return 'NA'
    if low == high or not low or not high:
        return high or low
    return f'{low}–{high}'


def noMatches(subtitle="Try a different query"):
    alfredItems([{
        "title": "No matches",
        "subtitle": subtitle,
        "valid": False,
        "icon": {"path": "icons/Warning.png"}
    }], {"myTextOutput": ""})


def emit(items, lines):
    """Print the item list along with the plain-text version used by the copy action."""
    lines = [f"**{MYTITLE}**"] + lines + [GWAS_REF]
    alfredItems(items, {"myTextOutput": "\n".join(lines)})


def showGenes():
    """Genes associated with the selected trait."""
    MYTRAIT = os.getenv('currentTrait', '')

    # `* 1` forces a numeric sort on databases built by older versions, where
    # these columns were stored as text.
    if hasFlag(MYENTRY, '--es'):
        orderS = "ORDER BY OR_Bmax IS NULL, OR_Bmax * 1 DESC"
    else:
        orderS = "ORDER BY PapCount * 1 DESC, pMax * 1 DESC"

    db = connect()
    rs = db.execute(f"SELECT * FROM GeneTrait WHERE trait = ? {orderS}", (MYTRAIT,)).fetchall()
    db.close()

    if not rs:
        noMatches(f"No genes recorded for {MYTRAIT}")
        return

    items, lines = [], []
    myResLen = len(rs)

    for countR, r in enumerate(rs, start=1):
        GeneName = r['GeneName'] or r['gene']
        Trait = r['trait']
        papList = r['PapList']
        locus = r['locus']
        keyCount = r['KeyCount']
        KeyList = r['KeyList']

        PapCount = r['PapCount']
        paperString = "paper" if (PapCount == 1) else "papers"

        pMax = fmt(r['pMax']) or 'NA'
        OR_B_block = orBetaBlock(r['OR_Bmin'], r['OR_Bmax'])

        itemString = (f"{GeneName}: {PapCount} {paperString} ({papList}), "
                      f"pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.)")
        myBIGFONT = (f"{MYTRAIT}-{GeneName} ({locus}): {PapCount} {paperString} ({papList}), "
                     f"pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.)")

        lines.append(f"\t{countR}. {itemString}")

        #### COMPILING OUTPUT
        items.append({
            "title": itemString,
            "subtitle": f"{countR}/{myResLen} {Trait} {locus} - ⬆️ for GTEx",
            "quicklookurl": f"https://gtexportal.org/home/gene/{GeneName}",
            "arg": "",
            "variables": {
                "currentTITLE": itemString,
                "mySource": "traitGene",
                "myAction": "",
                "myBIGFONT": myBIGFONT,
                "myKEYlist": KeyList,
                "myENTRY_Q": MYENTRY_Q
            },
            "icon": {
                "path": ""
            }
        })

    emit(items, lines)


def showTraits():
    """Traits associated with the selected gene."""
    MYGENE = os.getenv('currentTrait', '')  # gene name, i.e. the breadcrumb Alfred pre-fills
    if MYSOURCE == "geneMasterSearch":
        MYGENES = os.getenv('currentGeneID')
    else:
        MYGENES = os.getenv('currentGenes')

    if not MYGENES:
        alfredError("No gene selected", "Start again from the gene search")
        return

    if hasFlag(MYENTRY, '--es'):
        orderS = "ORDER BY OR_Bmax IS NULL, OR_Bmax * 1 DESC"
    else:
        orderS = "ORDER BY PapCount * 1 DESC, pMax * 1 DESC"

    # allow search refinement: drop the breadcrumb and any flags from the typed query
    MYSTRING = MYENTRY.replace(MYGENE, '') if MYGENE else MYENTRY
    MYSTRING = stripFlags(MYSTRING).strip()

    where, params = "gene = ?", [MYGENES]
    if MYSTRING:
        where += " AND trait LIKE ?"
        params.append(f"%{MYSTRING}%")

    db = connect()
    rs = db.execute(f"SELECT * FROM GeneTrait WHERE {where} {orderS}", params).fetchall()
    db.close()

    if not rs:
        noMatches(f"No traits for {MYGENE or MYGENES} matching '{MYSTRING}'"
                  if MYSTRING else f"No traits recorded for {MYGENE or MYGENES}")
        return

    items, lines = [], []
    myResLen = len(rs)

    for countR, r in enumerate(rs, start=1):
        GeneName = r['GeneName'] or r['gene']
        Trait = r['trait']
        PapCount = r['PapCount']
        paperString = "paper" if (PapCount == 1) else "papers"
        KeyList = r['KeyList']
        keyCount = r['KeyCount']
        papList = r['PapList']

        pMax = fmt(r['pMax']) or 'NA'
        OR_B_block = orBetaBlock(r['OR_Bmin'], r['OR_Bmax'])

        itemString = (f"{Trait}: {PapCount} {paperString} ({papList}), "
                      f"pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.)")
        myBIGFONT = (f"**{Trait}**-{GeneName}: {PapCount} {paperString} ({papList}), "
                     f"pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.) – {GWAS_REF}")

        lines.append(f"\t{countR}. {itemString}")

        #### COMPILING OUTPUT
        items.append({
            "title": itemString,
            "subtitle": f"{countR}/{myResLen} {Trait}",
            "variables": {
                "mySource": "geneTrait",
                "myAction": "",
                "myBIGFONT": myBIGFONT,
                "myKEYlist": KeyList,
                "myENTRY_Q": MYENTRY_Q
            },
            "arg": "",
            "icon": {
                "path": ""
            }
        })

    emit(items, lines)


def showPapers():
    """Individual associations supporting the selected gene-trait pair."""
    MYKEYS = []
    for key in (os.getenv('myKEYlist') or '').split(','):
        key = key.strip()
        if key:
            MYKEYS.append(int(key) if key.isdigit() else key)

    if not MYKEYS:
        noMatches("No associations to show")
        return

    if hasFlag(MYENTRY, '--es'):
        orderS = 'ORDER BY "OR or BETA" IS NULL, "OR or BETA" * 1 DESC'
    else:
        orderS = "ORDER BY DATE DESC"

    placeholders = ','.join(['?'] * len(MYKEYS))
    # select by name: the associations table mirrors the catalog file, whose
    # column order is not ours to rely on
    sql = (f'SELECT PUBMEDID, DATE, STUDY, "DISEASE/TRAIT", MAPPED_GENE, '
           f'PVALUE_MLOG, "OR or BETA" FROM associations '
           f'WHERE key IN ({placeholders}) {orderS}')

    db = connect()
    rs = db.execute(sql, MYKEYS).fetchall()
    db.close()

    if not rs:
        noMatches("No associations to show")
        return

    items, lines = [], []
    myResLen = len(rs)

    for countR, r in enumerate(rs, start=1):
        Study = r['STUDY']
        StudyDate = r['DATE'] or ''
        Trait = r['DISEASE/TRAIT']
        pubmedID = r['PUBMEDID']
        mappedGene = r['MAPPED_GENE']

        OR = fmt(r['OR or BETA'])
        pVal = fmt(r['PVALUE_MLOG'], '.1f') or 'NA'

        myTitle = f"{Trait}-({mappedGene}), {OR} ({pVal})"
        mySubTitle = f"{countR}/{myResLen}–{Study} ({pubmedID}, {StudyDate[0:4]})"

        lines.append(f"{countR}. {myTitle} {pubmedID}, {StudyDate[0:4]}")

        #### COMPILING OUTPUT
        items.append({
            "subtitle": mySubTitle,
            "title": myTitle,
            "variables": {
                "mySource": BACK_SOURCE_MAP.get(MYSOURCE, MYSOURCE),
                "myAction": "openPubMed",
                "myPUBMED": pubmedID,
                "myBIGFONT": mySubTitle,
                "myENTRY_Q": MYENTRY_Q
            },
            "arg": "",
            "icon": {
                "path": ""
            }
        })

    emit(items, lines)


HANDLERS = {
    "GWT": showGenes,          # source is the GWAS trait search
    "GWG": showTraits,         # source is the GWAS gene search
    "geneMasterSearch": showTraits,
    "geneTrait": showPapers,   # to get papers after gene > trait
    "traitGene": showPapers,   # to get papers after trait > gene
    "genePap": showPapers,
}


def main():
    requireDatabase()

    handler = HANDLERS.get(MYSOURCE)
    if handler is None:
        alfredError("Nothing to show", f"Unknown source: {MYSOURCE or '(none)'}")
        return

    handler()


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        traceback.print_exc(file=sys.stderr)
        alfredError(f"Error: {type(err).__name__}", str(err))
