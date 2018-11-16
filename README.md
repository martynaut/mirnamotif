# mirnamotif

miRNAmotif is an application that enables the user to 
*  search for specific motifs in the terminal loop or in whole pre-miRNA 
among sequences in a widely used miRNA database (miRBase) and 
* predict enriched motifs within sequences provided by the user. 

In particular, miRNAmotif is a simple and user-friendly application 
that gives researchers the opportunity to analyze known pre-miRNAs sequences 
for the presence of specific motifs that can be recognized by protein regulators 
of miRNA biogenesis.

## Licence

This application is open-source under MIT License (see LICENSE file)

## Availability

### Webserver
Webserved is freely available at [http://mirnamotif.ibch.poznan.pl/](http://mirnamotif.ibch.poznan.pl/)

### CLI
#### Database preparation
Files with miRNAs, linking sequence, loops and pre-miRNAs are already build from `hairpin.fa`, `hsa.gff3`,
 `mmu.gff3`, `ath.gff3` and `mature.fa`
but when needed they can be rebuild with new `hairpin.fa`, `hsa.gff3`,
 `mmu.gff3`, `ath.gff3` and `mature.fa` files (files need to be in `mirnamotif` folder) running
`python ./mirnamotif/parser.py`
#### Running mirnamotif - search for known motif

To find premirnas with single motif using mirnamotif run `python ./mirnamotif/mirnamotif.py GGAG` 
and `python ./mirnamotif/mirnamotif.py GGAG -s AGG` for two motifs
search.

To see options for this function please see `python mirnamotif/mirnamotif.py --help`

Options are:
`python ./mirnamotif/mirnamotif.py GGAG [-s AGG] [-t loops|sh|linking|mature]
[-d f|fr] [-o True|False] [-db hsa|mmu|ath] [-fp file_prefix]`

e.g.
`python ./mirnamotif/mirnamotif.py GGAG -s AGG -t linking
-d fr -o False -db mmu -fp mmu_linking`

Results will be saved in `./mirnamotif/results/` localization with day and hour within filename.

#### Running mirnamotif - search for new motif

To be able to search for new motifs meme suite is needed. Please download 4.11.4 version 
of meme http://web.mit.edu/meme_v4.11.4/share/doc/download.html into the same folder as 
repository, unpack the archive file and follow the instructions from INSTALL file.

After meme installation you can look for new motifs using command
`python ./mirnamotif/mirnamotif_find.py ./static/example_motif.txt -f ./results/new_meme_results -db sh`

First argument is a file directory with names of precursors (e.g. `mmu-let-7k`) or sequences.
Optional arguments include folder to save results, orientation and whether sequences searched 
should be loops or whole pre-miRNAs.

To see options please see `python ./mirnamotif/mirnamotif_find.py --help`.

Options are:
`python ./mirnamotif/mirnamotif_find.py ./static/example_motif.txt [-d f|fr] 
[-db sh|sh_ath|sh_mmu|sh_hsa|loops|loops_ath|loops_hsa|loops_mmu|mature|mature_ath|mature_hsa|mature_mmu|linking|linking_ath|linking_mmu|linking_hsa]
[-f ./results/new_meme_results] [-n directory_to_negative_sequences]`

Please remember that `-n` overwrites `-db` option.

e.g.
`python ./mirnamotif/mirnamotif_find.py ./static/example_motif.txt -d fr -db linking_hsa -f ./results/new_meme_results`

## Authors
Martyna O. Urbanek-Trzeciak

Edyta Jaworska

Wlodzimierz J. Krzyzosiak

## Citation

TBD
