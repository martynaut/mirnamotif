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
Files with miRNAs, loops and pre-miRNAs are already build from `hairpin.fa` 
but when needed they can be rebuild with new `hairpin.fa` file running
`python ./mirnamotif/parser.py`
#### Running mirnamotif - search for known motif

To find premirnas with single motif using mirnamotif run `python ./mirnamotif/mirnamotif.py GGAG` 
and `python ./mirnamotif/mirnamotif.py GGAG -s AGG` for two motifs
search.

To see options for this function please see `python mirnamotif/mirnamotif.py --help`

Results will be saved in `./mirnamotif/results/` localization with day and hour within filename.

#### Running mirnamotif - search for new motif

To be able to search for new motifs meme suite is needed. Please download 4.11.4 version 
of meme http://web.mit.edu/meme_v4.11.4/share/doc/download.html into the same folder as 
repository, unpack the archive file and follow the instructions from INSTALL file.

After meme installation you can look for new motifs using command
`python ./mirnamotif/mirnamotif_find.py ./static/example_motif.txt -f ./results/new_meme_results -db sh`

First argument is a file directory with names of precursors (e.g. `mmu-let-7k`) ro sequences.
Optional arguments include folder to save results, orientation and whether sequences searched 
should be loops or whole pre-miRNAs.

To see options please see `python ./mirnamotif/mirnamotif_find.py --help`.



## Authors
Martyna O. Urbanek-Trzeciak

Edyta Jaworska

Wlodzimierz J. Krzyzosiak

## Citation

TBD
