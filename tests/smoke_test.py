#!/usr/bin/env python3
"""Smoke tests for the alfred-GWAS query scripts.

Builds a small fixture database with the same schema the real build script
produces, then runs every entry point the workflow can reach and checks that
each one prints valid Alfred JSON and exits cleanly.

Run with:  python3 tests/smoke_test.py
No third-party dependencies (the query scripts are stdlib-only by design).
"""

import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), 'src')

# the associations table mirrors the catalog file, so the fixture mirrors it too
ASSOCIATION_COLUMNS = [
    'DATE ADDED TO CATALOG', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'JOURNAL', 'LINK',
    'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE', 'REPLICATION SAMPLE SIZE',
    'REGION', 'CHR_ID', 'CHR_POS', 'REPORTED GENE(S)', 'MAPPED_GENE',
    'UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS',
    'UPSTREAM_GENE_DISTANCE', 'DOWNSTREAM_GENE_DISTANCE', 'STRONGEST SNP-RISK ALLELE',
    'SNPS', 'MERGED', 'SNP_ID_CURRENT', 'CONTEXT', 'INTERGENIC',
    'RISK ALLELE FREQUENCY', 'P-VALUE', 'PVALUE_MLOG', 'P-VALUE (TEXT)', 'OR or BETA',
    '95% CI (TEXT)', 'PLATFORM [SNPS PASSING QC]', 'CNV', 'MAPPED_TRAIT',
    'MAPPED_TRAIT_URI', 'STUDY ACCESSION', 'GENOTYPING TECHNOLOGY', 'key',
    'ImplicatedGenes',
]

APOSTROPHE_TRAIT = "Crohn's disease"


def buildFixture(path):
    """A database exercising the awkward cases: apostrophes, NULLs, missing annotation."""
    db = sqlite3.connect(path)

    db.execute("CREATE TABLE colophon (colophon TEXT, timestamp TEXT)")
    db.execute("INSERT INTO colophon VALUES ('e115_r2026-02-16_full', '2026-03-02 at 10:00')")

    db.execute("""CREATE TABLE geneCounts
                  (gene TEXT, nTraits INTEGER, GeneName TEXT, searchField TEXT,
                   PUBMEDID TEXT, nPapers INTEGER)""")
    db.executemany("INSERT INTO geneCounts VALUES (?,?,?,?,?,?)", [
        ('ENSG00000001', 2, 'NOD2', 'NOD2,CARD15,ENSG00000001', '123,456', 2),
        # a gene with no annotation row, as ~40% of the catalog is
        ('ENSG00000002', 1, None, None, '789', 1),
    ])

    db.execute("""CREATE TABLE traitCounts
                  (MAPPED_TRAIT TEXT, ImplicatedGenes TEXT, ImplicatedGenes_count INTEGER,
                   PubmedIDs TEXT, papers_count INTEGER)""")
    db.executemany("INSERT INTO traitCounts VALUES (?,?,?,?,?)", [
        (APOSTROPHE_TRAIT, 'ENSG00000001', 1, '123', 1),
        ('glucose measurement', 'ENSG00000002', 1, '789', 1),
    ])

    db.execute("""CREATE TABLE GeneTrait
                  (trait TEXT, gene TEXT, OR_Bmax REAL, OR_Bmin REAL, KeyList TEXT,
                   pMax REAL, PapList TEXT, locus TEXT, KeyCount INTEGER,
                   PapCount INTEGER, GeneName TEXT, searchField TEXT)""")
    db.executemany("INSERT INTO GeneTrait VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", [
        (APOSTROPHE_TRAIT, 'ENSG00000001', 1.4, 1.1, '1,2', 12.5, '123,456',
         '16q12', 2, 2, 'NOD2', 'NOD2,ENSG00000001'),
        # NULL effect size and NULL p-value, on a gene with no annotation
        ('glucose measurement', 'ENSG00000002', None, None, '3', None, '789',
         '1p31', 1, 1, None, None),
    ])

    columns = ','.join(f'"{c}" TEXT' for c in ASSOCIATION_COLUMNS)
    db.execute(f"CREATE TABLE associations ({columns})")

    def association(key, pmlog, orbeta, date):
        row = [''] * len(ASSOCIATION_COLUMNS)
        row[1] = '12345678'
        row[3] = date
        row[6] = 'Some Study'
        row[7] = APOSTROPHE_TRAIT
        row[14] = 'NOD2'
        row[28] = pmlog
        row[30] = orbeta
        row[38] = key
        return row

    db.executemany(
        "INSERT INTO associations VALUES (%s)" % ','.join(['?'] * len(ASSOCIATION_COLUMNS)),
        [association(1, 12.5, 1.4, '2020-01-01'),
         association(2, 9.0, 1.1, '2021-05-05'),
         # one row with no p-value and no effect size, as the catalog contains
         association(3, None, None, '2019-01-01'),
         # oldest row but the largest effect: date order and effect-size order
         # disagree, so a broken ORDER BY cannot pass by luck
         association(4, 7.0, 2.5, '2018-03-03')])

    db.commit()
    db.close()


def run(script, argument='', **env):
    """Run a query script the way Alfred does and return its parsed output."""
    environment = dict(os.environ)
    environment['alfred_workflow_data'] = WF_DATA
    environment.update({k: v for k, v in env.items() if v is not None})
    for key in ('mySource', 'currentTrait', 'currentGenes', 'currentGeneID',
                'currentTITLE', 'myKEYlist', 'myENTRY_Q'):
        if key not in env:
            environment.pop(key, None)

    proc = subprocess.run([sys.executable, script, argument], cwd=SRC,
                          capture_output=True, text=True, env=environment)
    return proc


FAILURES = []


def check(name, proc, expect_items=None, title_contains=None, absent=None,
          first_title_contains=None):
    problems = []

    if proc.returncode != 0:
        problems.append(f"exit code {proc.returncode}")

    payload = None
    try:
        payload = json.loads(proc.stdout)
    except (ValueError, TypeError):
        problems.append(f"stdout is not JSON: {proc.stdout[:120]!r}")

    if payload is not None:
        items = payload.get('items')
        if not isinstance(items, list) or not items:
            problems.append("no items in payload")
        else:
            if expect_items is not None and len(items) != expect_items:
                problems.append(f"expected {expect_items} items, got {len(items)}")
            titles = ' | '.join(str(i.get('title', '')) for i in items)
            if title_contains and title_contains not in titles:
                problems.append(f"expected {title_contains!r} in titles: {titles[:160]!r}")
            if absent and absent in titles:
                problems.append(f"unexpected {absent!r} in titles: {titles[:160]!r}")
            first = str(items[0].get('title', ''))
            if first_title_contains and first_title_contains not in first:
                problems.append(f"expected {first_title_contains!r} first, got {first!r}")

    if 'Traceback' in proc.stderr:
        problems.append(f"traceback on stderr: {proc.stderr.strip().splitlines()[-1]}")

    if problems:
        FAILURES.append((name, problems))
        print(f"FAIL  {name}")
        for problem in problems:
            print(f"        {problem}")
    else:
        print(f"ok    {name}")


def main():
    global WF_DATA
    WF_DATA = tempfile.mkdtemp(prefix='alfred-gwas-test-')
    try:
        buildFixture(os.path.join(WF_DATA, 'index.db'))

        # --- gene search
        check("gene search finds an annotated gene",
              run('query-gene.py', 'NOD2'), expect_items=1, title_contains='NOD2')
        check("gene search finds an unannotated gene by its Ensembl id",
              run('query-gene.py', 'ENSG00000002'), expect_items=1,
              title_contains='ENSG00000002')
        check("gene search survives an apostrophe in the query",
              run('query-gene.py', "Crohn's"), title_contains='No matches')
        check("gene search treats a bare --p as a sort flag, not a search term",
              run('query-gene.py', '--p'), expect_items=2)
        check("gene search handles no argument",
              run('query-gene.py', ''), expect_items=2)

        # --- trait search
        check("trait search finds a trait containing an apostrophe",
              run('query-trait.py', "Crohn"), expect_items=1,
              title_contains=APOSTROPHE_TRAIT)
        check("trait search reports no matches cleanly",
              run('query-trait.py', 'zzzz'), title_contains='No matches')

        # --- trait -> gene
        check("trait drill-down works for a trait containing an apostrophe",
              run('query-items.py', '', mySource='GWT', currentTrait=APOSTROPHE_TRAIT,
                  currentTITLE=APOSTROPHE_TRAIT), expect_items=1, title_contains='NOD2')
        check("trait drill-down renders NULL p-value and NULL effect size",
              run('query-items.py', '', mySource='GWT', currentTrait='glucose measurement',
                  currentTITLE='glucose measurement'), expect_items=1, absent='None')
        check("trait drill-down sorts by effect size with --es",
              run('query-items.py', '--es', mySource='GWT', currentTrait=APOSTROPHE_TRAIT,
                  currentTITLE=APOSTROPHE_TRAIT), expect_items=1)

        # --- gene -> trait
        check("gene drill-down lists traits",
              run('query-items.py', '', mySource='GWG', currentTrait='NOD2',
                  currentGenes='ENSG00000001', currentTITLE='NOD2'),
              expect_items=1, title_contains=APOSTROPHE_TRAIT)
        check("gene drill-down refines by substring",
              run('query-items.py', 'Crohn', mySource='GWG', currentTrait='NOD2',
                  currentGenes='ENSG00000001', currentTITLE='NOD2'), expect_items=1)
        check("gene drill-down reports an unmatched refinement instead of crashing",
              run('query-items.py', 'zzzz', mySource='GWG', currentTrait='NOD2',
                  currentGenes='ENSG00000001', currentTITLE='NOD2'),
              title_contains='No matches')
        check("gene drill-down treats --es as a flag, not a search term",
              run('query-items.py', '--es', mySource='GWG', currentTrait='NOD2',
                  currentGenes='ENSG00000001', currentTITLE='NOD2'), expect_items=1)
        check("gene drill-down survives an apostrophe typed as a refinement",
              run('query-items.py', "Crohn's", mySource='GWG', currentTrait='NOD2',
                  currentGenes='ENSG00000001', currentTITLE='NOD2'), expect_items=1)
        check("gene drill-down honours currentGeneID from the master gene search",
              run('query-items.py', '', mySource='geneMasterSearch',
                  currentGeneID='ENSG00000001', currentTrait='NOD2', currentTITLE='NOD2'),
              expect_items=1, title_contains=APOSTROPHE_TRAIT)

        # --- papers
        check("papers list renders a row with no p-value and no effect size",
              run('query-items.py', '', mySource='geneTrait', myKEYlist='1,2,3',
                  currentTITLE='NOD2'), expect_items=3, absent='None')
        check("papers list sorts by date by default",
              run('query-items.py', '', mySource='geneTrait', myKEYlist='1,2,4',
                  currentTITLE='NOD2'), expect_items=3, first_title_contains='1.10')
        check("papers list sorts by effect size with --es",
              run('query-items.py', '--es', mySource='geneTrait', myKEYlist='1,2,4',
                  currentTITLE='NOD2'), expect_items=3, first_title_contains='2.50')
        check("papers list handles an empty key list",
              run('query-items.py', '', mySource='geneTrait', myKEYlist='',
                  currentTITLE='NOD2'), title_contains='No matches')

        # --- missing database
        empty = tempfile.mkdtemp(prefix='alfred-gwas-empty-')
        try:
            previous, globals()['WF_DATA'] = WF_DATA, empty
            check("a missing database asks the user to rebuild",
                  run('query-gene.py', 'NOD2'), title_contains='not available')
            check("a missing database does not break the drill-down either",
                  run('query-items.py', '', mySource='GWT', currentTrait='x'),
                  title_contains='not available')
            globals()['WF_DATA'] = previous
        finally:
            shutil.rmtree(empty, ignore_errors=True)

    finally:
        shutil.rmtree(WF_DATA, ignore_errors=True)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == '__main__':
    sys.exit(main())
