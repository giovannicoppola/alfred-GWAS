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


import sqlite3
import json
import sys
import os


from config import INDEX_DB, log, databaseReady, NOT_BUILT
MYENTRY_Q = os.getenv('myENTRY_Q') # this is breadcrumbs to enable the 'back' feature
MYSOURCE = os.getenv('mySource')

def queryTraits ():
    if not databaseReady():
        print(json.dumps(NOT_BUILT))
        return
    db = sqlite3.connect(INDEX_DB)
    cursor = db.cursor()

    
    argv1 = sys.argv[1] if len(sys.argv) > 1 else ""
    if MYSOURCE in ["traitGene"] and argv1 == '':
        MYINPUT= MYENTRY_Q or ""

    else:
        MYINPUT= argv1


    

    MYQUERY= "%" + MYINPUT + "%"

    result = {"items": [], "variables":{}}

    try:
        cursor.execute(f"""SELECT *
            FROM traitCounts
            WHERE MAPPED_TRAIT LIKE ? 
            ORDER BY papers_count DESC, ImplicatedGenes_count DESC""", (MYQUERY,))
        
        rs = cursor.fetchall()
    
    except sqlite3.OperationalError as err:
        resultErr = {"items": [{
        "title": "Error: " + str(err),
        "subtitle": "Some error",
        "arg": "",
        "icon": {

                "path": "icons/Warning.png"
            }
        }]}
        print (json.dumps(resultErr))
        raise err

    if (rs):
        myResLen = str(len (rs))
        
        countR=1
        
        for r in rs:
            
            
            FeatureName = r[0]
            
            PapersCount = r[4]
            paperString = "paper" if (PapersCount == 1) else "papers"

            
            ImplicatedGenes = r[1]
            GeneCount = r[2]
            geneString = "gene" if (GeneCount == 1) else "genes"


            subtitleString = (str(countR)+"/"+myResLen +
            f" – {GeneCount} {geneString}, from {PapersCount} {paperString}"
            
                
            )
            
            
            bigText = f"{FeatureName} – {GeneCount:,} {geneString}, from {PapersCount:,} {paperString}"
            titleString = f"{FeatureName}: {GeneCount:,} {geneString}, {PapersCount:,} {paperString}"
            
        
        
            #### COMPILING OUTPUT    
            result["items"].append({
            "title": titleString,
            "subtitle": subtitleString,
                
            "arg": FeatureName,
            
            "variables": {
            "currentTrait": FeatureName,
            "myENTRY_Q": MYINPUT,
            "currentGenes": ImplicatedGenes,
            "currentTITLE": titleString,
            "bigText": bigText,
            "mySource": "GWT"
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
    

    queryTraits ()

if __name__ == '__main__':
    main ()
