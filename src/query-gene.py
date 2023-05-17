#!/usr/bin/env python3

### Gene QUERY

#### Tuesday, June 7, 2022, 5:37 PM
# Partly cloudy ⛅️  🌡️+71°F (feels +71°F, 63%) 🌬️↑17mph 🌓 Tue Jun  7 17:37:33 2022
# W23Q2 – 158 ➡️ 206 – 85 ❇️ 279

## revision for v0.3
#NYP – Clear ☀️   🌡️+55°F (feels +54°F, 26%) 🌬️→6mph 🌓 Mon Mar 27 00:55:03 2023
#W13Q1 – 86 ➡️ 278 – 320 ❇️ 44

# checking summaries from the GWAS catalog, by gene 


import sqlite3
import json
import sys
import os

from config import INDEX_DB, log
MYSOURCE = os.getenv('mySource')
MYENTRY_Q = os.getenv('myENTRY_Q')
MYGENE = os.getenv('currentGeneID')

def queryGenes ():
    result = {"items": [], "variables":{}}


    db = sqlite3.connect(INDEX_DB)
    cursor = db.cursor()

    if MYSOURCE in ["GWG"] and sys.argv[1] == '':
        MYINPUT= MYENTRY_Q

    elif MYSOURCE in ["geneMasterSearch"]:
        MYINPUT= MYGENE

    else:
        MYINPUT= sys.argv[1].strip()


    orderS = 'ORDER BY nTraits DESC'


    if "--p" in MYINPUT:
        orderS = "ORDER BY nPapers DESC "
        MYINPUT = MYINPUT.replace (' --p','')

    MYQUERY = "%" + MYINPUT + "%"



    try:
        cursor.execute(f"""SELECT *
        FROM geneCounts
        WHERE searchField LIKE '{MYQUERY}' {orderS}""")
        
        rs = cursor.fetchall()
    

        

    except sqlite3.OperationalError as err:
        result= {"items": [{
        "title": "Error: " + str(err),
        "subtitle": "Some error",
        "arg": "",
        "icon": {

                "path": "icons/Warning.png"
            }
        }]}
        print (json.dumps(result))
        raise err

    if (rs):
        myResLen = str(len (rs))
        countR=1
        
        for r in rs:
            
            
            FeatureName = r[2]
            geneID = r[0]
            
            PapersCount = r[5]
            paperString = "paper" if (PapersCount == 1) else "papers"

            TraitCount = r[1]
            traitString = "trait" if (TraitCount == 1) else "traits"
            
            
            subtitleString = (str(countR)+"/"+myResLen +
            f" – associated with {TraitCount:,} {traitString}, from {PapersCount:,} {paperString}"
                
                
            )
            
            titleString = f"{FeatureName}: {TraitCount:,} {traitString} ({PapersCount:,} {paperString})"
            
        
        
            #### COMPILING OUTPUT    
            result["items"].append({
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
            countR += 1  

                

        print (json.dumps(result))

    if MYQUERY and not rs:
        resultErr= {"items": [{
            "title": "No matches",
            "subtitle": "Try a different query",
            "arg": "",
            "icon": {
                "path": "icons/Warning.png"
                }
            
                }]}
        print (json.dumps(resultErr))



def main():
    
    queryGenes ()

if __name__ == '__main__':
    main ()


