## README

This repository contains a number of *.py files used to process output data 
from bioinformatics processes (eggnog, Interproscan). It also contains code 
used for processing orthologous sequence (OG) data for species spanning the 
eukaryote tree of life (orthogroup_analysis).

The main functions necessary for analysing OG data are contained within group.py. 
Uniquely shared OG data can be converted to a heatmap (og_to_heatmap.py), unique OGs 
belonging to sets of 3+ eukaryote groups can be scanned (setquery.py), and shared 
genome content of 3 eukaryote groups can be visualised with venn3_euk_groups.py. 

covid_sequence_analysis contains code used to determine the significance of amino 
acid frequency changes along alignments of the COVID spike protein to highlight potential 
mutations.
