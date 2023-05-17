#!/usr/bin/env python3

    
### ITEMS-QUERY 
# showing individual items (i.e. papers, genes, loci)

#### Tuesday, May 31, 2022, 5:05 PM
# Partly cloudy ⛅️  🌡️+76°F (feels +80°F, 69%) 🌬️→12mph 🌑 Tue May 31 07:24:28 2022
# W22Q2 – 151 ➡️ 213 – 20 ❇️ 345

import sqlite3
import json
import sys
import os
import webbrowser
from config import INDEX_DB, log, GWAS_REF



db = sqlite3.connect(INDEX_DB)
cursor = db.cursor()
result = {"items": [], "variables":{}}

MYSOURCE = os.getenv('mySource')
MYENTRY = sys.argv[1]
MYTITLE = os.getenv('currentTITLE')

def showGenes ():
    MYTRAIT = os.getenv('currentTrait')
    
    #MYGENES = os.getenv('currentGenes')
    orderS = " ORDER BY PapCount*1 DESC, pMax*1 DESC"
    
    
    if "--es" in MYENTRY:
        orderS = " ORDER BY OR_Bmax*1 DESC"
    
    
    sql = f"""SELECT * FROM GeneTrait 
            WHERE trait = '{MYTRAIT}' 
            {orderS}"""
            
    cursor.execute(sql)
    rs = cursor.fetchall()

    
    if (rs):
        myResLen = len (rs)
        countR=1
    
    for r in rs:
        
        GeneName = r[10] 
        Trait = r[0] 
        #GeneLocus = r[3]
        
        papList = papListAll = r[6]
        
        
        
        locus = r[7]
        keyCount = r[8]
        KeyList = r[4] 

        PapCount = r[9]
        paperString = "paper" if (PapCount == 1) else "papers"
        

        pMax = "{:.2f}".format(float(r[5]))
        
        if r[3]:
            OR_BetaMin = f"{r[3]:.2f}"
        else:
            OR_BetaMin = ''
    
        if r[2]:
            OR_BetaMax = f"{r[2]:.2f}"
        else:
            OR_BetaMax = ''
    
        
        if keyCount == 1:
            OR_B_block = OR_BetaMax
        else:
            OR_B_block = f'{OR_BetaMin}–{OR_BetaMax}'

        itemString = f"{GeneName}: {PapCount} {paperString} ({papList}), pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.)"
        myBIGFONT = f"{MYTRAIT}-{GeneName} ({locus}): {PapCount} {paperString} ({papListAll}), pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.)"
        if countR == 1:
            myTextOutput = f"**{MYTITLE}**\n"   
            
        myTextOutput = (myTextOutput 
            + f"\t{countR}"
            + ". "
            + itemString
            + "\n"
                
                )    
        if countR == myResLen:
            myTextOutput = f"{myTextOutput}{GWAS_REF}" 
        
    #### COMPILING OUTPUT    
        result["items"].append({
        "title": itemString,
        "subtitle": f"{countR}/{myResLen} {Trait} {locus} - ⬆️ for GTEx",
        "quicklookurl": f"https://gtexportal.org/home/gene/{GeneName}",
        "arg": "",
        "variables": {
        "currentTITLE": itemString,
        "mySource": "traitGene",
        "myBIGFONT": myBIGFONT, 
        "myKEYlist": KeyList
        },
        
        "icon": {   
        
        "path": ""
    }
        })
        countR += 1  

    
            
    result['variables'] = {"myTextOutput": myTextOutput}   
    print (json.dumps(result))

def showTraits ():
    MYGENE = os.getenv('currentTrait')
    if "mySource" == "geneMasterSearch":
        MYGENES = os.getenv("currentGeneID")
    else:
        MYGENES = os.getenv('currentGenes')

    
    orderS = " ORDER BY PapCount*1 DESC, pMax DESC"
    
    MYSTRING = sys.argv[1].replace(MYGENE,'').strip() #to allow search refinement
    
    if MYSTRING:
        SQL_SUBSTRING = f" AND trait LIKE '%{MYSTRING}%'"
    else:
        SQL_SUBSTRING = ''

    if "--es" in MYENTRY:
        orderS = " ORDER BY OR_Bmax*1 DESC"
        #MYENTRY = MYENTRY.replace ('--es','')
    
    db = sqlite3.connect(INDEX_DB)
    db.row_factory = sqlite3.Row
    cursor = db.cursor()
    sql = f"SELECT * FROM GeneTrait WHERE gene = '{MYGENES}' {SQL_SUBSTRING} {orderS}"
    
    cursor.execute(sql)
    rs = db.execute(sql).fetchall()
        

    
    if (rs):
        myResLen = len (rs)
        countR=1
    
    for r in rs:
        
        #GeneName = r[10] 
        GeneName = r['GeneName']

        Trait = r['trait'] 
        
        PapCount = r['PapCount']
        paperString = "paper" if (PapCount == 1) else "papers"
        KeyList = r['KeyList'] 
        keyCount = r['KeyCount']
        papList = papListAll = r['papList']
        


        pMax = "{:.2f}".format(float(r['pMax']))
        
        if r['OR_Bmin']:
            OR_BetaMin = "{:.2f}".format(float (r['OR_Bmin']))
        else:
            OR_BetaMin = ''
        if r['OR_Bmax']:
            OR_BetaMax = "{:.2f}".format(float (r['OR_Bmax']))
        else:
            OR_BetaMax = ''
        
        if keyCount == 1:
            OR_B_block = OR_BetaMax
        else:
            OR_B_block = f'{OR_BetaMin}–{OR_BetaMax}'

        if countR == 1:
            myTextOutput = f"**{MYTITLE}**\n"
        
        myTextOutput = (myTextOutput 
            + f"\t{countR}. "
            + f"{Trait}: {PapCount} {paperString} ({papList}), pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.)\n"  
                
                )
        if countR == myResLen:
            myTextOutput = f"{myTextOutput}{GWAS_REF}" 

        myBIGFONT = f"**{Trait}**-{GeneName}: {PapCount} {paperString} ({papListAll}), pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.) – {GWAS_REF}"
    
    #### COMPILING OUTPUT    
        result["items"].append({
        "title": f"{Trait}: {PapCount} {paperString} ({papList}), pMax: {pMax}, OR/B: {OR_B_block} ({keyCount} assoc.)",
        "subtitle": f"{countR}/{myResLen} {Trait}",
        "variables": {
        "mySource": "geneTrait",
        "myBIGFONT": myBIGFONT, 
        "myKEYlist": KeyList
        },    
        "arg": "",
        "icon": {  
        
        "path": ""
    }
        })
        countR += 1  

    
            
    result['variables'] = {"myTextOutput": myTextOutput}           
    print (json.dumps(result))



def showPapers (): 
    orderS = " ORDER BY DATE DESC"
    MYTITLE = os.getenv('currentTITLE')
    MYKEYS = os.getenv('myKEYlist')
    MYKEYS = MYKEYS.split(',')

    if "--es" in MYENTRY:
        orderS = " ORDER BY 'OR or BETA'*1 DESC"
        #MYENTRY = MYENTRY.replace ('--es','')
    
    
    sql = "SELECT * FROM associations WHERE key IN ({seq})".format(seq=','.join(['?']*len(MYKEYS))) + orderS
    cursor.execute(sql, MYKEYS)
    rs = cursor.fetchall()


    if (rs):
        myResLen = len (rs)
        countR=1
    
    for r in rs:
        
        Study = r[6] 
        StudyDate = r[3] 
        Trait = r[7] 
        pubmedID = r[1]
        mappedGene = r[14] 
        if r[30]:
            OR =  f'{float(r[30]):.2f}'
        else:
            OR = ''
        
        pVal = float(r[28])
        myTitle = f"{Trait}-({mappedGene}), {OR} ({pVal:.1f})"
        mySubTitle = f"{countR}/{myResLen}–{Study} ({pubmedID}, {StudyDate[0:4]})"

        if countR == 1:
            myTextOutput = f"**{MYTITLE}** \n"  
         
        myTextOutput = (myTextOutput 
            + str(countR)
            + ". "
            + myTitle + " "
            + pubmedID + ", " + StudyDate[0:4]
            + "\n"
            
            )
        if countR == myResLen:
            myTextOutput = f"{myTextOutput}{GWAS_REF}" 


    #### COMPILING OUTPUT    
        result["items"].append({
        "subtitle": mySubTitle,
        "title": myTitle,
        
         "variables": {
        "mySource": "papers",
        "myPUBMED": pubmedID,
        "myBIGFONT": mySubTitle
        },            
        "arg": "",
        "icon": {   
        
        "path": ""
    }
        }) 
        countR += 1  


    
    result['variables'] = {"myTextOutput": myTextOutput}        
                
    print (json.dumps(result))

    


def main():


    if MYSOURCE == "GWT": #source is the GWAS trait search
        showGenes ()


    elif MYSOURCE == "GWG": #source is the GWAS gene search
        showTraits ()

    elif MYSOURCE == "geneTrait": # to get papers after gene > trait
        showPapers ()
        
    elif MYSOURCE == "traitGene": #to get papers after trait >gene
        showPapers ()
        
    elif MYSOURCE == "genePap":
        showPapers ()

    elif MYSOURCE == "papers":
        MYPUBMED = os.getenv('myPUBMED')
        webbrowser.open(f'https://pubmed.ncbi.nlm.nih.gov/{MYPUBMED}', new=2)


if __name__ == '__main__':
    main ()




