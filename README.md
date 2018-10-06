# mirnamotif

miRNAmotif is an application that enables the user to 
*  search for specific motifs in the terminal loop or in whole pre-miRNA 
among sequences in a widely used miRNA database (miRBase) and 
* predict enriched motifs within sequences provided by the user. 

In particular, miRNAmotif is a simple and user-friendly application 
that gives researchers the opportunity to analyze known pre-miRNAs sequences 
for the presence of specific motifs that can be recognized by protein regulators 
of miRNA biogenesis.

## Availability

### Webserver
Webserved is freely available at [http://mirnamotif.ibch.poznan.pl/](http://mirnamotif.ibch.poznan.pl/)

### CLI
#### Database preparation
Files with miRNAs, loops and pre-miRNAs are already build from `hairpin.fa` 
but when needed they can be rebuild with new `hairpin.fa` file running
`python ./mirnamotif/parser.py`
#### Running mirnamotif


## Authors
Martyna O. Urbanek-Trzeciak

Edyta Jaworska

Wlodzimierz J. Krzyzosiak

## Citation

TBD