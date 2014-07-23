# Region-Motif-Compare Tools
### Information
Version 1.1

Released 2014
Park Laboratory
Center for Biomedical Informatics
Harvard University

### Contact
Jeremy Liu (jeremy.liu@yale.edu)
Nils Gehlenb0rg (nils@hms.harvard.edu)

## Overview
### Structure

### Description



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
### Running as Command Line Tools

### Running from Galaxy and Refinery

## Motif Tabix File Creation
