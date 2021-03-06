# sequence_recovery.py
This script compares transcript sequence alignments of orthologous zebra finch and chicken [Ensembl](http://www.ensembl.org) gene models to the zebra finch (*Taeniopygia guttata*) genome assembly [taeGut1](https://genome.ucsc.edu/cgi-bin/hgGateway?db=taeGut1). 
### Background
The project in this repository was designed to analyze the ion channel genes assessed in [Friedrich et al. 2019 BMC Genomics](https://doi.org/10.1186/s12864-019-5871-2). This publication is publically available and a detailed description of the data curation and analysis pipeline can be found in the **Methods** section under the subheading "**Assessing gene model completeness and expanding gene models with additional sequence**".

In short, after finding that many zebra finch gene models were incomplete relative to other species including chicken, my goal was to use alignments of chicken models to the zebra finch genome to identify sequence that may be missing from current zebra finch gene models. For each gene, all blocks of aligned chicken model sequence that did not overlap with the orthologous zebra finch model were considered recovered sequence, representing potential coding regions of the genome that could extend current zebra finch models. 

The product of this script is a text file in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) which allows for visualization of recovered sequence blocks using UCSC's [taeGut1](https://genome.ucsc.edu/cgi-bin/hgGateway?db=taeGut1) Genome Browser.

## Getting started
Script (sequence_recovery.py) and input files should be placed in the same folder when running this code. The output file (output_BED.txt) is what should be produced after running.

### Python library prerequisites
Please be sure the following are installed in your python environment:
- os
- subprocess
- pandas
- re

### Input files
- ENS_key.txt - a text file linking HGNC gene names to zebra finch and chicken gene models
- alignments-to-taeGut1.psl - a psl file containing alignments of chicken and zebra finch Ensembl gene models to taeGut1 generated using the [UCSC Genome Browser Blat](https://genome.ucsc.edu/index.html)
- gene_to_transcript.txt - a text file linking Ensembl gene ID's to Ensembl transcript ID's

### Output file
- output_BED.txt - a BED track of the chicken sequence blocks that aligned in taeGut1 and did not overlap any of the corresponding zebra finch model(s)

## Author
Sami Friedrich, PhD candidate at Oregon Health and Science University
- [LinkedIn](https://www.linkedin.com/in/sami-friedrich/)
 
## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/samifriedrich/seq-recovery/blob/master/sequence_recovery.py) file for details.


[![DOI](https://zenodo.org/badge/171545573.svg)](https://zenodo.org/badge/latestdoi/171545573)

