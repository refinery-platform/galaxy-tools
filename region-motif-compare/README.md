# Region-Motif-Compare Tools
Version 1.1 Released 2014  
Park Laboratory  
Center for Biomedical Informatics  
Harvard University  

Contact:  
Jeremy Liu (jeremy.liu@yale.edu)  
Nils Gehlenborg (nils@hms.harvard.edu)

## Overview
### Structure
The tool suite consists of:  
1. Two Rscripts: region_motif_compare.r and region_motif_intersect.r
2. Two Xml Files: region_motif_compare.xml and region_motif_intersect.xml
3. Motif Database Directory: region_motif_db
4. Dependency Library Directory: region_motif_lib
5. Galaxy Workflows: Galaxy-Workflow-Region_Motif_Count_Comparison_Test_Motifs.ga

### Description
1. region_motif_intersect.r: takes one bed file of regions as input. Then it calculates
the number of intersections of the regions and the motifs. region_motifs_intersect.r
outputs a tab separated values (tsv) file of motif names and intersection counts.

2. region_motif_compare.r: takes as input two tsv files of motifs / regions intersection
counts. These generally originate from running region_motif_intersect.r on two sets
of different regions with the same query motif database. Based on the counts, 
region_motif_compare.r then determines the enrichment (or depletion) of certain
motifs across the two regions. This is done by a correcting for the size and gc
content of the region, and applying a Poisson test to the counts. 
Then, region_motif_compare.r outputs the most significant enriched or depleted
motifs as a tsv. In addition, the tool outputs a diagnostic plot containing
graphical representations of the motif counts, gc correction curves, and significant 
motifs that distinguish the two regions (selected via p value).

3. Motif positions are stored in region_motif_db as compressed, index tabix files.  
Dependencies (i.e. plotting.r) are stored in region_motif_lib.

## Installation
Directions for installing the region motif tools into your local galaxy distribution

1) In ~/galaxy-dist/tools/ make a directory called "my_tools" and cd into it
2) Clone the relevant files from the github into the directory
3) In ~/galaxy-dist/tool_conf.xml add the following lines:
  <section id="mTools" name="My Tools">
    <tool file="my_tools/region_motif_intersect.xml" />
    <tool file="my_tools/region_motif_compare.xml" />
  </section>
4) In ~/galaxy-dist/tools/my_tools/region_motif_lib/ you will have to run the 
following gcc commands to compile the shared library.
```
  gcc  -I/usr/local/include -fPIC  -g -O2 -c regions.cpp -o regions.o
  gcc  -shared -o regions.so regions.o -L/usr/lib64/R/lib -lR -lstdc++
```
5) In region_motif_compare.r and region_motif_intersect.r you may have to
adjust the file paths of commonDir or workingDir if it complains at you.

## Running the Tools
### Running from Galaxy and Refinery

### Running as Command Line Tools

## Motif Tabix File Creation
