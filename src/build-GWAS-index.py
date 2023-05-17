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
    GWAS-cat <file_name>
    GWAS-cat (-h|--help)
    GWAS-cat (-V|--version)

Options:
    -V, --version                   show version number and exit
    -h, --help                      show this message and exit

    -q, --quiet                     show error messages and above
    -v, --verbose                   show info messages and above
    -d, --debug                     show debug messages
"""


__version__ = "version 0.3, March 2023"
__author__ = "giovanni"


import sqlite3
from time import time
import pandas as pd
import numpy as np
import os
import datetime
import logging

import re


DEFAULT_LOG_LEVEL = logging.WARNING #log level

# initializing the log file    
timeStr = datetime.datetime.now().strftime ("%Y-%m-%d-%H-%M")
timeSta = datetime.datetime.now().strftime ("%Y-%m-%d at %H:%M")
myTimeStamp = f"Script run on {timeSta}"
LOG_FILE = f"logs/{timeStr}_log.md"
ANNOTATION_DB = 'geneNames.db'
INDEX_DB = 'index.db'


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
    print (f"\t 1. reading the master associations file in", end = "...")
    myData = pd.read_csv(DATA_FILE, sep='\t', header=0, dtype = "str")
    print ("done")
    print (f"\t 2. creating the sqlite database", end = "...")
    con = sqlite3.connect(INDEX_DB)
    cursor = con.cursor()
  

    myData['key'] = range(1, len(myData.index)+1)
    myData['PVALUE_MLOG'] = myData['PVALUE_MLOG'].astype('float')
    myData['OR or BETA'] = myData['OR or BETA'].astype('float')
    myData = myData.fillna('')

    cursor.execute(f"DROP TABLE IF EXISTS associations")
    myData.to_sql("associations", con, index=False)
    

    print ("done")

    colophon = makeColophon (DATA_FILE)
    cursor.execute(f"DROP TABLE IF EXISTS colophon")
    
    # Create the table
    cursor.execute(f"CREATE TABLE colophon (colophon TEXT, timestamp TEXT)")

    # Insert the values into the table
    cursor.execute(f"INSERT INTO colophon (colophon, timestamp) VALUES (?, ?)", (colophon, timeSta))

    # Commit the changes and close the connection
    con.commit()
    con.close ()

    return myData
    



def unique_set_of_comma_separated_values(s):
    return ','.join(pd.unique(s.split(',')))

def count_comma_separated_elements(s):
    if s:
        return len(s.split(','))
    else:
        return 0

def createTraitCounts (myDataFrame):
    
    # pbar = tqdm(total=5) # 5 steps
    # pbar.update(1)  
    
    # Importing the gene annotation table from the gene lookup DB
    connAnn = sqlite3.connect(ANNOTATION_DB)

    # Read the SQL table into a Pandas dataframe
    geneAnnotation = pd.read_sql_query("SELECT GeneName, EnsemblGeneId, searchField FROM geneAnnotation", connAnn)
    # Close the database connection
    connAnn.close()

    con = sqlite3.connect(INDEX_DB)
    cursor = con.cursor()
    
    #adding ensembl ID to the search 

    # define a function to concatenate column A and B if both are non-empty, otherwise use value of A
    # def concat_or_a(row):
    #     if pd.notna(row['searchField']):
    #         return row['searchField'] + ',' + row['EnsemblGeneId']
    #     else:
    #         return row['EnsemblGeneId']

    # # apply the function to each row of the dataframe and create a new column C with the result
    # geneAnnotation['searchField'] = geneAnnotation.apply(concat_or_a, axis=1)



    

    cursor.execute(f"DROP TABLE IF EXISTS geneAnnotation")
    geneAnnotation.to_sql("geneAnnotation", con, index=False)
    

    




    ### 1. Pubmed data summary by trait
    print (f"\n\t1. Pubmed and Gene Summary by trait", end = "...")
    #counting the number of PUBMED IDs per trait
    count_data = myDataFrame.groupby('MAPPED_TRAIT')['PUBMEDID'].nunique().reset_index()
    #merging PUBMED ID numbers per trait
    pub_data = myDataFrame.groupby('MAPPED_TRAIT')['PUBMEDID'].apply(lambda x: ', '.join(x.astype(str).unique())).reset_index()

    pub_data = pub_data.merge(count_data, on='MAPPED_TRAIT')

    
    # Renaming the columns
    pub_data.columns = ['MAPPED_TRAIT', 'PubmedIDs','papers_count']

    print ("done")
    
    ## 2. implicated genes summary
    
    #pbar.update(1) #step 2
    print (f"\n\t2. implicated genes summary", end = "...")
    # Group the data by 'MAPPED_TRAIT' and combine the non-null values in 'UPSTREAM', 'DOWNSTREAM', and 'SNP' columns
    myDataFrame = myDataFrame.replace('', np.nan)
    gene_data = myDataFrame.groupby('MAPPED_TRAIT').agg({'UPSTREAM_GENE_ID': lambda x: ','.join(set(x.dropna())),
                                                       'DOWNSTREAM_GENE_ID': lambda x: ','.join(set(x.dropna())),
                                                       'SNP_GENE_IDS': lambda x: ','.join(set(x.dropna()))
                                                       })
    # joining in a single field
    gene_data['ImplicatedGenes'] = gene_data.apply(lambda row: ','.join([str(x) for x in row.values if x != '']), axis=1)
    #counting
    gene_data['ImplicatedGenes_count'] = gene_data['ImplicatedGenes'].apply(count_comma_separated_elements)

    gene_data = gene_data.drop(['UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID','SNP_GENE_IDS'], axis=1)
    merged_data = gene_data.merge(pub_data,on='MAPPED_TRAIT')
    
    
    
    cursor.execute(f"DROP TABLE IF EXISTS traitCounts")
    merged_data.to_sql("traitCounts", con, index=False)
    
    print (f"done")

    ## geneCounts table
    #pbar.update(1) #step 4
    print (f"\n\t3. geneCounts table", end = "...")
    
    #print (merged_data.head(10))
    # Split the 'ImplicatedGenes' column into multiple rows
    merged_data['ImplicatedGenes'] = merged_data['ImplicatedGenes'].str.split(',')
    df = merged_data.explode('ImplicatedGenes')
    
    # Count the number of 'MAPPED_TRAIT' associated with each gene
    df['ImplicatedGenes'] = df['ImplicatedGenes'].str.strip()
    geneCounts = df.groupby(['ImplicatedGenes']).agg({
                                                    'MAPPED_TRAIT': lambda x: x.nunique()
                                                    
                                                       })
    #geneCounts['nPapers'] = geneCounts['PubmedIDs'].apply(count_comma_separated_elements)
    
    
    geneCounts = geneCounts.reset_index()
    geneCounts.columns.values[0] = 'gene'
    geneCounts.columns.values[1] = 'nTraits'
    #geneCounts.columns.values[2] = 'paperList'
    
    geneCounts = geneCounts.replace('', np.nan)
    geneCounts = geneCounts.dropna(subset=['gene'])

    geneCounts = pd.merge(geneCounts, geneAnnotation, left_on='gene',right_on='EnsemblGeneId', how='left' )
    geneCounts = geneCounts.drop('EnsemblGeneId', axis=1)
    #print (geneCounts.head(20))
    
    print (f"done")
    

    
    ## Pubmed
    #pbar.update(1) #step 3
    print (f"\n\t4. genes implicated in Pubmed", end = "...")
    gene_pub = myDataFrame.groupby('PUBMEDID').agg({'UPSTREAM_GENE_ID': lambda x: ','.join(set(x.dropna())),
                                                       'DOWNSTREAM_GENE_ID': lambda x: ','.join(set(x.dropna())),
                                                       'SNP_GENE_IDS': lambda x: ','.join(set(x.dropna()))
                                                       })
    gene_pub['ImplicatedGenes'] = gene_pub.apply(lambda row: ','.join([str(x) for x in row.values if x != '']), axis=1)
    
    gene_pub['ImplicatedGenes'] = gene_pub['ImplicatedGenes'].apply(unique_set_of_comma_separated_values)
    gene_pub['ImplicatedGenes_count'] = gene_pub['ImplicatedGenes'].apply(count_comma_separated_elements)
    
    gene_pub = gene_pub.drop(['UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID','SNP_GENE_IDS'], axis=1).reset_index()
    
    #print (gene_pub.head(20))
    cursor.execute(f"DROP TABLE IF EXISTS paperCounts")
    gene_pub.to_sql("paperCounts", con, index=False)
    print (f"done")


    # completing geneCounts table with paper counts
    gene_pub['ImplicatedGenes'] = gene_pub['ImplicatedGenes'].str.split(',')
    gene_pub = gene_pub.explode('ImplicatedGenes')
    gene_pub['ImplicatedGenes'] = gene_pub['ImplicatedGenes'].str.strip() 

    gene_pub['PUBMEDID'] = gene_pub['PUBMEDID'].astype(str)
    countPub = gene_pub.groupby(['ImplicatedGenes']).agg({
                                                       
                                                       'PUBMEDID': lambda x: ','.join(set(x.dropna()))
                                                       })
    countPub['nPapers'] = countPub['PUBMEDID'].apply(count_comma_separated_elements)
    
    geneCounts = geneCounts.merge(countPub,left_on='gene',right_on='ImplicatedGenes',how='left')
    
    cursor.execute(f"DROP TABLE IF EXISTS geneCounts")
    geneCounts.to_sql("geneCounts", con, index=False)
    


    #######
    ### locus-Trait pairs
    ######
    ## 3 gene-trait table
    #pbar.update(1) #step 5
    print (f"\n\t5. gene-trait table", end = "...")
    pair_data = myDataFrame.groupby('key').agg({'UPSTREAM_GENE_ID': lambda x: ','.join(set(x.dropna())),
                                                       'DOWNSTREAM_GENE_ID': lambda x: ','.join(set(x.dropna())),
                                                       'SNP_GENE_IDS': lambda x: ','.join(set(x.dropna())),
                                                       'OR or BETA': lambda x: x.iloc[0], #it is ok to take the first record, because grouping on 'key' will return one value only
                                                       'MAPPED_TRAIT': lambda x: x.iloc[0],
                                                       'PVALUE_MLOG': lambda x: x.iloc[0],
                                                       'PUBMEDID': lambda x: x.iloc[0],
                                                       'REGION': lambda x: x.iloc[0]
                                                       })
    pair_data['ImplicatedGenes'] = pair_data.loc[:, ['UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS']].apply(lambda row: ','.join([str(x) for x in row.values if x != '']), axis=1)
    pair_data = pair_data.drop(['UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID','SNP_GENE_IDS'], axis=1).reset_index()
    pair_data = pair_data.replace('', np.nan)
    pair_data = pair_data.dropna(subset=['ImplicatedGenes'])
    
    # Split the 'ImplicatedGenes' column into multiple rows
    #print (pair_data.head(50))
    pair_data['ImplicatedGenes'] = pair_data['ImplicatedGenes'].str.split(',')
    pair_data = pair_data.explode('ImplicatedGenes')
    pair_data['ImplicatedGenes'] = pair_data['ImplicatedGenes'].str.strip()
    #print (pair_data.head(50))
    pair_data['key'] = pair_data['key'].astype(str)
    pair_data['PUBMEDID'] = pair_data['PUBMEDID'].astype(str)
    pair_dataAGG = pair_data.groupby(['MAPPED_TRAIT','ImplicatedGenes']).agg({
                                                       'OR or BETA': ['max','min'],
                                                       'key': ','.join,
                                                       'PVALUE_MLOG': lambda x: x.max(),
                                                       'PUBMEDID': lambda x: ','.join(set(x.dropna())),
                                                       'REGION': lambda x: x.iloc[0]
                                                       })
    
    # renaming columns
    pair_dataAGG.columns = pair_dataAGG.columns.droplevel(0)
    pair_dataAGG.columns.values[0] = 'OR_Bmax'
    pair_dataAGG.columns.values[1] = 'OR_Bmin'
    pair_dataAGG.columns.values[2] = 'KeyList'
    pair_dataAGG.columns.values[3] = 'pMax'
    pair_dataAGG.columns.values[4] = 'PapList'
    pair_dataAGG.columns.values[5] = 'locus'
    pair_dataAGG = pair_dataAGG.reset_index()
    pair_dataAGG.columns.values[0] = 'trait'
    pair_dataAGG.columns.values[1] = 'gene'

    #counting keys and papers
    pair_dataAGG['KeyCount'] = pair_dataAGG['KeyList'].apply(count_comma_separated_elements)
    pair_dataAGG['PapCount'] = pair_dataAGG['PapList'].apply(count_comma_separated_elements)
    
    pair_dataAGG = pd.merge(pair_dataAGG, geneAnnotation, left_on='gene',right_on='EnsemblGeneId', how='left' )
    pair_dataAGG = pair_dataAGG.drop('EnsemblGeneId', axis=1)
    

    cursor.execute(f"DROP TABLE IF EXISTS GeneTrait")
    pair_dataAGG.to_sql("GeneTrait", con, index=False)
    
    print (f"done")


    
    #pbar.close()    

    
    


def main(args=None):
    """Run program."""


    from docopt import docopt
    main_start_time = time()
    args = docopt(__doc__, version=__version__)
    # Access the file_name argument
    DATA_FILE = args['<file_name>']

    if args.get('--verbose'):
        logging.basicConfig(level=logging.INFO)
        #log.setLevel(logging.INFO)
    elif args.get('--quiet'):
        logging.basicConfig(level=logging.ERROR)
        #log.setLevel(logging.ERROR)
    elif args.get('--debug'):
        logging.basicConfig(level=logging.DEBUG)
        #log.setLevel(logging.DEBUG)
    else:
        logging.basicConfig(level=DEFAULT_LOG_LEVEL)
        #log.setLevel(DEFAULT_LOG_LEVEL)

    myCurrentLogLevel = logging.getLevelName(logging.getLogger().getEffectiveLevel())
    logging.info (args)
    logging.info("Set log level to {}".format(myCurrentLogLevel)) 
    

    if os.path.exists(INDEX_DB):
       os.remove(INDEX_DB)

    
    logF (f"# GWAS catalog build\n{myTimeStamp}", file_name =  LOG_FILE)
    logF (f"file used: {DATA_FILE}\n", file_name =  LOG_FILE)

    
    
# importing the master association table into the database
    
    print (f"1. creating the master associations table")
    myData = associations_table(INDEX_DB,DATA_FILE)
    
# outputting basic stats
    n = myData['PUBMEDID'].nunique()
    print (f"number of unique studies: {n:,}")
    logF (f"- number of unique studies: {n:,}", file_name =  LOG_FILE)
          
    n = myData['DISEASE/TRAIT'].nunique()
    print (f"number of unique disease/traits: {n:,}")
    logF (f"- number of unique disease/traits: {n:,}", file_name =  LOG_FILE)
    
    n = myData['MAPPED_TRAIT'].nunique()
    print (f"number of unique mapped traits: {n:,}")
    logF (f"- number of unique mapped traits: {n:,}", file_name =  LOG_FILE)
    
    
    
# generate trait counts using 1) the master dataset, 2) the list of unique traits, 3) the association table source dictionary
    print (f"\n2. generating trait and gene counts ...")
    
    createTraitCounts (myData)
    
    
    
    
    main_timeElapsed = time() - main_start_time
    print(f"\nTotal script duration: {round (main_timeElapsed,2)} seconds")
    logF (f"\nTotal script duration: {round (main_timeElapsed,2)} seconds, {main_timeElapsed/60:.1f} minutes", file_name =  LOG_FILE)



    

   
    

if __name__ == '__main__':
    main ()

