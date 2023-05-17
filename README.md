# alfred-GWAS 🧬
An [Alfred](https://www.alfredapp.com/) workflow to browse the GWAS catalog



<a href="https://github.com/giovannicoppola/alfred-GWAS/releases/latest/">
<img alt="Downloads"
src="https://img.shields.io/github/downloads/giovannicoppola/alfred-GWAS/total?color=purple&label=Downloads"><br/>
</a>

![](images/alfred-GWAS.gif)

<!-- MarkdownTOC autolink="true" bracket="round" depth="3" autoanchor="true" -->

- [Motivation](#motivation)
- [Setting up](#setting-up)
- [Basic Usage](#usage)
- [Known Issues](#known-issues)
- [Acknowledgments](#acknowledgments)
- [Changelog](#changelog)
- [Feedback](#feedback)

<!-- /MarkdownTOC -->


<h1 id="motivation">Motivation ✅</h1>

- to obtain, without opening websites, a summarized version of all genes associated with a particular trait in the GWAS catalog, or all traits asssociated with a gene or locus.
- quickly access papers supporting a gene-trait association



<h1 id="setting-up">Setting up ⚙️</h1>

### Needed
- [Alfred 5](https://www.alfredapp.com/) with Powerpack license


<h1 id="usage">Basic Usage 📖</h1>

## Querying by gene 🧬
- launch with keyword (default: `gwg`) or custom hotkey, enter a search string. Results will include a gene name, the number of associated traits, and the number of papers.
- output is sorted by number of traits by default. Adding `--p` to the search string will sort by number of papers. 
- once a gene is selected, associated traits are shown, including for each trait the range of OR (or beta), minimum p-value, number of papers supporting the association, and number of associated loci. It is possible to refine the list of associated traits by entering an additional search string. 

	- `ctrl-enter` will show the currently selected gene-trait pair in large font, and copy to clipboard
	- `cmd-enter` will copy the entire gene-trait list to clipboard
	- `cmd-option-enter` will go back to the previous list
- once a gene-trait pair is selected, individual associations are shown, with OR min p-value, and source publication. 
 	- `ctrl-enter` will show the currently selected association in large font, and copy to clipboard
	- `cmd-enter` will copy the entire set of associations to clipboard
- once an individual association is selected, `enter` will open the supporting publication in Pubmed. 	

	
## Querying by trait 👤
- launch with keyword (default: `gwt`) or custom hotkey. Results will include a gene name, the number of associated traits, and the number of papers. 
- once a gene is selected, associated traits are shown, including for each trait the range of OR (or beta), minimum p-value, number of papers supporting the association, and number of associated loci.
- output is sorted by number of supporting papers. Adding `--es` to the search string will sort by largest reported effect size (OR or beta). 
	- `ctrl-enter` will show the currently selected gene-trait pair in large font, and copy to clipboard
	- `cmd-enter` will copy the entire gene-trait list to clipboard
	- `cmd-option-enter` will go back to the previous list  
- once a gene-trait pair is selected, individual associations are shown, with OR min p-value, and source publication. 
 	- `ctrl-enter` will show the currently selected association in large font, and copy to clipboard
	- `cmd-enter` will copy the entire set of associations to clipboard
- once an individual association is selected, `enter` will open the supporting publication in Pubmed.  



## rebuilding database (optional) 🛠️
- `alfred-GWAS` comes with the database derived from version `e109_r2023-05-07` of the GWAS catalog (v1.0.2, all associations) available [here](https://www.ebi.ac.uk/gwas/docs/file-downloads).
- the script `build-GWAS-index.py` can be used to rebuild/update the database. Run with `python3 build-GWAS-index.py "path-to-GWAS-file.tsv"`, or select the file in Finder and use the Universal Action. 


<h1 id="known-issues">Limitations & known issues ⚠️</h1>

- None for now, but I have not done extensive testing, let me know if you see anything!



<h1 id="acknowledgments">Acknowledgments 😀</h1>

- Icons from [Flaticon](https://www.flaticon.com/)
	
	
<h1 id="changelog">Changelog 🧰</h1>

- 05-17-2023: version 0.4
- 03-25-2023: version 0.3
- 05-31-2022: version 0.2
- 12-12-2020: version 0.1


<h1 id="feedback">Feedback 🧐</h1>

Feedback welcome! If you notice a bug, or have ideas for new features, please feel free to get in touch either here, or on the [Alfred](https://www.alfredforum.com) forum. 

