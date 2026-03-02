#!/usr/bin/env python3

# Saturday, December 12, 2020: GWAS lookup
# Version 2: Thursday, May 26, 2022, 3:31 PM (current script length: 21 min)
# trying to create a more efficient version from NYP, Saturday, March 25, 2023 (duration: sec)
# Light rain 🌦   🌡️+43°F (feels +35°F, 73%) 🌬️↖9mph 🌒 Sat Mar 25 13:23:29 2023
# W12Q1 – 84 ➡️ 280 – 318 ❇️ 46


# either install docopt (`pip3 install docopt`), or have the docopt.py file bundled in the same directory

"""GWAScat-builder

Process the GWAS catalog so that it can be used in an Alfred Workflow
Main steps:
1. Read and import into a database the master GWAS list (associations_table() function, `associations` table in the database)
1a. logging basic stats
2. for each unique mapped trait:
    2a. count number of implicated papers
    2b. count number of implicated genes



Usage:
    GWAS-cat [<file_name>]
    GWAS-cat (-h|--help)
    GWAS-cat (-V|--version)

Options:
    -V, --version                   show version number and exit
    -h, --help                      show this message and exit

    -q, --quiet                     show error messages and above
    -v, --verbose                   show info messages and above
    -d, --debug                     show debug messages

If no <file_name> is provided, the script downloads the latest associations
file directly from the GWAS Catalog (EBI).
"""


__version__ = "version 0.4, March 2026"
__author__ = "giovanni"


import sqlite3
from time import time
import pandas as pd
import numpy as np
import os
import datetime
import logging
import urllib.request
import sys
import zipfile
import re


DEFAULT_LOG_LEVEL = logging.WARNING #log level

# initializing the log file
timeStr = datetime.datetime.now().strftime ("%Y-%m-%d-%H-%M")
timeSta = datetime.datetime.now().strftime ("%Y-%m-%d at %H:%M")
myTimeStamp = f"Script run on {timeSta}"

WF_DATA = os.getenv('alfred_workflow_data', '.')
LOG_DIR = os.path.join(WF_DATA, 'logs')
os.makedirs(LOG_DIR, exist_ok=True)
LOG_FILE = os.path.join(LOG_DIR, f"{timeStr}_log.md")

ANNOTATION_DB = 'geneNames.db'
INDEX_DB = os.path.join(WF_DATA, 'index.db')
GWAS_DOWNLOAD_URL = 'https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false'


def get_server_filename():
    """Query the GWAS catalog server for the current release filename."""
    req = urllib.request.Request(GWAS_DOWNLOAD_URL, method='HEAD')
    with urllib.request.urlopen(req) as response:
        content_disp = response.headers.get('Content-Disposition', '')
    # e.g. "attachement; filename=gwas_catalog_v1.0.2-associations_e115_r2026-02-16_full.zip"
    match = re.search(r'filename=(.+)', content_disp)
    if match:
        return match.group(1).strip()
    return None


def download_gwas_catalog():
    """Download the GWAS catalog associations file from EBI.

    Returns (tsv_path, already_existed) tuple.
    """
    server_filename = get_server_filename()
    # derive a TSV name from the zip name (e.g. gwas_catalog_v1.0.2-associations_e115_r2026-02-16_full.tsv)
    tsv_name = server_filename.replace('.zip', '.tsv') if server_filename else f"associations_{timeStr}.tsv"
    dest_tsv = os.path.join(WF_DATA, tsv_name)

    if os.path.exists(dest_tsv):
        print(f"File already exists: {dest_tsv}", file=sys.stderr)
        return dest_tsv, True

    print(f"Downloading GWAS catalog from EBI...", file=sys.stderr)
    print(f"\tServer file: {server_filename}", file=sys.stderr)

    # download the zip
    dest_zip = os.path.join(WF_DATA, server_filename)
    urllib.request.urlretrieve(GWAS_DOWNLOAD_URL, dest_zip)

    # extract the TSV from the zip
    with zipfile.ZipFile(dest_zip, 'r') as zf:
        tsv_inside = zf.namelist()[0]
        zf.extract(tsv_inside, WF_DATA)
        extracted_path = os.path.join(WF_DATA, tsv_inside)
        os.rename(extracted_path, dest_tsv)

    # clean up the zip
    os.remove(dest_zip)

    print(f"\tDownload complete: {dest_tsv}", file=sys.stderr)
    return dest_tsv, False


def log(s, *args):
    if args:
        s = s % args
    print(s, file=sys.stderr)


def logF(log_message, file_name):
    with open(file_name, "a") as f:
        f.write(log_message + "\n")

def makeColophon(DATA_FILE):

    match = re.search(r'associations_(.+?)\.tsv', DATA_FILE)

    if match:
        colophon = match.group(1)

    else:
        colophon = timeStr
    return colophon


def associations_table(INDEX_DB, DATA_FILE):
    """Create the master associations table from the master file"""

    # reading the master file in:
    print (f"\t 1. reading the master associations file in", end = "...", file=sys.stderr)
    myData = pd.read_csv(DATA_FILE, sep='\t', header=0, low_memory=False,
                         dtype={'PVALUE_MLOG': 'float64', 'OR or BETA': 'float64',
                                'PUBMEDID': str, 'CHR_ID': str, 'CHR_POS': str,
                                'UPSTREAM_GENE_ID': str, 'DOWNSTREAM_GENE_ID': str,
                                'SNP_GENE_IDS': str})
    # fill NaN in string columns only (keep NaN in float columns for proper aggregation)
    str_cols = myData.select_dtypes(include='object').columns
    myData[str_cols] = myData[str_cols].fillna('')
    print ("done", file=sys.stderr)
    print (f"\t 2. creating the sqlite database", end = "...", file=sys.stderr)
    con = sqlite3.connect(INDEX_DB)
    cursor = con.cursor()

    myData['key'] = range(1, len(myData.index)+1)

    # Pre-compute ImplicatedGenes column (vectorized, used by all downstream steps)
    # Concatenate gene columns with comma, then clean up leading/trailing/double commas
    gene_cols = ['UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS']
    raw = myData[gene_cols[0]].str.strip().str.cat(
        [myData[c].str.strip() for c in gene_cols[1:]], sep=','
    )
    # Remove empty segments from concatenation (e.g. ",," becomes ",")
    myData['ImplicatedGenes'] = raw.str.replace(r',+', ',', regex=True).str.strip(',')

    print ("done", file=sys.stderr)

    colophon = makeColophon (DATA_FILE)

    # Write colophon first (small, fast)
    cursor.execute(f"DROP TABLE IF EXISTS colophon")
    cursor.execute(f"CREATE TABLE colophon (colophon TEXT, timestamp TEXT)")
    cursor.execute(f"INSERT INTO colophon (colophon, timestamp) VALUES (?, ?)", (colophon, timeSta))
    con.commit()
    con.close ()

    return myData


def unique_set_of_comma_separated_values(s):
    return ','.join(set(s.split(',')))

def count_comma_separated_elements(s):
    if s:
        return len(s.split(','))
    else:
        return 0

def createTraitCounts (myDataFrame):

    # Importing the gene annotation table from the gene lookup DB
    connAnn = sqlite3.connect(ANNOTATION_DB)
    geneAnnotation = pd.read_sql_query("SELECT GeneName, EnsemblGeneId, searchField FROM geneAnnotation", connAnn)
    connAnn.close()

    con = sqlite3.connect(INDEX_DB)
    cursor = con.cursor()

    cursor.execute(f"DROP TABLE IF EXISTS geneAnnotation")
    geneAnnotation.to_sql("geneAnnotation", con, index=False)

    # ImplicatedGenes is pre-computed in associations_table(); explode once and reuse
    exploded = myDataFrame[['key', 'MAPPED_TRAIT', 'PUBMEDID', 'OR or BETA',
                            'PVALUE_MLOG', 'REGION', 'ImplicatedGenes']].copy()
    exploded = exploded[exploded['ImplicatedGenes'] != '']
    exploded['ImplicatedGenes'] = exploded['ImplicatedGenes'].str.split(',')
    exploded = exploded.explode('ImplicatedGenes')
    exploded['ImplicatedGenes'] = exploded['ImplicatedGenes'].str.strip()
    exploded = exploded[exploded['ImplicatedGenes'] != '']


    ### 1. Pubmed + Gene Summary by trait
    print (f"\n\t1. Pubmed and Gene Summary by trait", end = "...", file=sys.stderr)

    pub_agg = myDataFrame.groupby('MAPPED_TRAIT')['PUBMEDID'].agg(
        PubmedIDs=lambda x: ', '.join(x.astype(str).unique()),
        papers_count='nunique'
    ).reset_index()

    # Gene summary by trait — use pre-exploded data
    gene_by_trait = exploded.groupby('MAPPED_TRAIT')['ImplicatedGenes'].agg(
        ImplicatedGenes=lambda x: ','.join(set(x))
    ).reset_index()
    gene_by_trait['ImplicatedGenes_count'] = gene_by_trait['ImplicatedGenes'].str.count(',') + 1

    merged_data = gene_by_trait.merge(pub_agg, on='MAPPED_TRAIT')

    cursor.execute(f"DROP TABLE IF EXISTS traitCounts")
    merged_data.to_sql("traitCounts", con, index=False)

    print (f"done", file=sys.stderr)

    ### 2. geneCounts table
    print (f"\n\t2. geneCounts table", end = "...", file=sys.stderr)

    geneCounts = exploded.groupby('ImplicatedGenes')['MAPPED_TRAIT'].nunique().reset_index()
    geneCounts.columns = ['gene', 'nTraits']
    geneCounts = geneCounts[geneCounts['gene'] != '']

    geneCounts = geneCounts.merge(geneAnnotation, left_on='gene', right_on='EnsemblGeneId', how='left')
    geneCounts = geneCounts.drop('EnsemblGeneId', axis=1)

    print (f"done", file=sys.stderr)

    ### 3. genes implicated in Pubmed (paperCounts + complete geneCounts)
    print (f"\n\t3. genes implicated in Pubmed", end = "...", file=sys.stderr)

    # Paper-level gene summary
    pub_genes = exploded.groupby('PUBMEDID')['ImplicatedGenes'].agg(
        lambda x: ','.join(set(x))
    ).reset_index()
    pub_genes.columns = ['PUBMEDID', 'ImplicatedGenes']
    pub_genes['ImplicatedGenes'] = pub_genes['ImplicatedGenes'].apply(unique_set_of_comma_separated_values)
    pub_genes['ImplicatedGenes_count'] = pub_genes['ImplicatedGenes'].str.count(',') + 1

    cursor.execute(f"DROP TABLE IF EXISTS paperCounts")
    pub_genes.to_sql("paperCounts", con, index=False)

    # Complete geneCounts with paper counts per gene
    countPub = exploded.groupby('ImplicatedGenes')['PUBMEDID'].agg(
        PUBMEDID=lambda x: ','.join(set(x.astype(str)))
    ).reset_index()
    countPub['nPapers'] = countPub['PUBMEDID'].str.count(',') + 1

    geneCounts = geneCounts.merge(countPub, left_on='gene', right_on='ImplicatedGenes', how='left')
    geneCounts = geneCounts.drop('ImplicatedGenes', axis=1)

    cursor.execute(f"DROP TABLE IF EXISTS geneCounts")
    geneCounts.to_sql("geneCounts", con, index=False)

    print (f"done", file=sys.stderr)

    ### 4. gene-trait table (no more redundant groupby on unique key)
    print (f"\n\t4. gene-trait table", end = "...", file=sys.stderr)

    # exploded already has one gene per row with all needed columns
    pair_data = exploded[['key', 'MAPPED_TRAIT', 'PUBMEDID', 'OR or BETA',
                          'PVALUE_MLOG', 'REGION', 'ImplicatedGenes']].copy()
    pair_data['key'] = pair_data['key'].astype(str)
    pair_data['PUBMEDID'] = pair_data['PUBMEDID'].astype(str)

    pair_dataAGG = pair_data.groupby(['MAPPED_TRAIT', 'ImplicatedGenes']).agg(
        OR_Bmax=('OR or BETA', 'max'),
        OR_Bmin=('OR or BETA', 'min'),
        KeyList=('key', ','.join),
        pMax=('PVALUE_MLOG', 'max'),
        PapList=('PUBMEDID', lambda x: ','.join(set(x.dropna()))),
        locus=('REGION', 'first')
    ).reset_index()
    pair_dataAGG.columns = ['trait', 'gene', 'OR_Bmax', 'OR_Bmin', 'KeyList', 'pMax', 'PapList', 'locus']

    pair_dataAGG['KeyCount'] = pair_dataAGG['KeyList'].str.count(',') + 1
    pair_dataAGG['PapCount'] = pair_dataAGG['PapList'].str.count(',') + 1

    pair_dataAGG = pair_dataAGG.merge(geneAnnotation, left_on='gene', right_on='EnsemblGeneId', how='left')
    pair_dataAGG = pair_dataAGG.drop('EnsemblGeneId', axis=1)

    cursor.execute(f"DROP TABLE IF EXISTS GeneTrait")
    pair_dataAGG.to_sql("GeneTrait", con, index=False)

    print (f"done", file=sys.stderr)

    ### 5. Write the associations table (deferred to end to not block computation)
    print (f"\n\t5. writing associations table", end = "...", file=sys.stderr)

    cursor.execute(f"DROP TABLE IF EXISTS associations")
    myDataFrame.to_sql("associations", con, index=False)

    print (f"done", file=sys.stderr)

    ### 6. Create indexes for query performance
    print (f"\n\t6. creating indexes", end = "...", file=sys.stderr)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_traitCounts_trait ON traitCounts(MAPPED_TRAIT)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_geneCounts_gene ON geneCounts(gene)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_geneCounts_name ON geneCounts(GeneName)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_geneCounts_search ON geneCounts(searchField)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_GeneTrait_trait ON GeneTrait(trait)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_GeneTrait_gene ON GeneTrait(gene)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_geneAnnotation_ensembl ON geneAnnotation(EnsemblGeneId)")
    con.commit()

    print (f"done", file=sys.stderr)


def main(args=None):
    """Run program."""


    from docopt import docopt
    main_start_time = time()
    args = docopt(__doc__, version=__version__)
    # Access the file_name argument, or download if not provided
    DATA_FILE = args['<file_name>']
    if DATA_FILE is None:
        DATA_FILE, already_existed = download_gwas_catalog()

    if args.get('--verbose'):
        logging.basicConfig(level=logging.INFO)
    elif args.get('--quiet'):
        logging.basicConfig(level=logging.ERROR)
    elif args.get('--debug'):
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=DEFAULT_LOG_LEVEL)

    myCurrentLogLevel = logging.getLevelName(logging.getLogger().getEffectiveLevel())
    logging.info (args)
    logging.info("Set log level to {}".format(myCurrentLogLevel))

    try:

        if os.path.exists(INDEX_DB):
           os.remove(INDEX_DB)


        logF (f"# GWAS catalog build\n{myTimeStamp}", file_name =  LOG_FILE)
        logF (f"file used: {DATA_FILE}\n", file_name =  LOG_FILE)



    # importing the master association table into the database

        print (f"1. creating the master associations table", file=sys.stderr)
        myData = associations_table(INDEX_DB,DATA_FILE)

    # outputting basic stats
        n = myData['PUBMEDID'].nunique()
        print (f"number of unique studies: {n:,}", file=sys.stderr)
        logF (f"- number of unique studies: {n:,}", file_name =  LOG_FILE)

        n = myData['DISEASE/TRAIT'].nunique()
        print (f"number of unique disease/traits: {n:,}", file=sys.stderr)
        logF (f"- number of unique disease/traits: {n:,}", file_name =  LOG_FILE)

        n = myData['MAPPED_TRAIT'].nunique()
        print (f"number of unique mapped traits: {n:,}", file=sys.stderr)
        logF (f"- number of unique mapped traits: {n:,}", file_name =  LOG_FILE)



    # generate trait counts using 1) the master dataset, 2) the list of unique traits, 3) the association table source dictionary
        print (f"\n2. generating trait and gene counts ...", file=sys.stderr)

        createTraitCounts (myData)




        main_timeElapsed = time() - main_start_time
        print(f"\nTotal script duration: {round (main_timeElapsed,2)} seconds", file=sys.stderr)
        logF (f"\nTotal script duration: {round (main_timeElapsed,2)} seconds, {main_timeElapsed/60:.1f} minutes", file_name =  LOG_FILE)

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)


if __name__ == '__main__':
    main ()
