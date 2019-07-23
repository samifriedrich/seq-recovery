# seq-recovery
This script compares alignments of zebra finch and chicken [Ensembl](http://www.ensembl.org) gene models to the zebra finch (*Taeniopygia guttata*) genome assembly [taeGut1](https://genome.ucsc.edu/cgi-bin/hgGateway?db=taeGut1) with the goal of recovering sequence that may be missing from zebra finch gene predictions (i.e. aligned chicken sequence blocks that do not overlap with any of the corresponding zebra finch model).

The project in this repository was designed to analyze zebra finch ion channel gene models in Friedrich et al. 2019 BMC Genomics, and a more detailed description of the data and analysis can be found in the **Methods** section of this publication under the subheading "Assessing gene model completeness and expanding gene models with additional sequence".

## Getting started
Script (sequence_recovery.py) and input files should be placed in the same folder when running this code. The output file (output_BED.txt) is what should be produced after running.

### Python library prerequisites
Please be sure the following are already installed in your python environment:
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
 - Sami Friedrich, PhD candidate at Oregon Health and Science University
 
## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/samifriedrich/seq-recovery/blob/master/sequence_recovery.py) file for details.


[![DOI](https://zenodo.org/badge/171545573.svg)](https://zenodo.org/badge/latestdoi/171545573)

