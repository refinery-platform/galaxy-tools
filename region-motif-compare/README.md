# Region-Motif-Compare Tools
Version 1.2.0 Released 2015  
Park Laboratory  
Center for Biomedical Informatics  
Harvard University  

Contact  
Jeremy Liu (jeremy.liu@yale.edu)  
Nils Gehlenborg (nils@hms.harvard.edu)

## Overview
### Structure
This tool suite consists of:

1. Two Rscripts: region_motif_compare.r and region_motif_intersect.r
2. Two Xml Files: region_motif_compare.xml and region_motif_intersect.xml
3. Motif Database Configuration Files: motif_databases.loc.sample and 
tool_data_table_conf.xml.sample
4. Dependency Configuration Files: plotting.r, repository_dependencies.xml, 
and tool_dependencies.xml
5. Test Data: Located in test-data, this contains test BED files from the ENCODE
project that can be used to verify correct tool functioning.
6. Galaxy Workflows: Files with suffix ".ga" that can be imported into the local
Galaxy instance after installation of the tool.

### Description
1. **region_motif_intersect.r** (1 bed -> 1 tsv): 
Takes one bed file of regions as input. Then it calculates
the number of intersections of the regions and the motifs. region_motifs_intersect.r
outputs a tab separated values (tsv) file of motif names and intersection counts.  
**Important Note:** region_motif_intersect.r makes no assumptions about the nature
of the input regions. For example, if overlapping regions are inputted, motifs that
intersect the overlap will be double counted. Thus, it is recommended that regions
be merged before using this tool, using the merge tool in the Galaxy toolshed.

2. **region_motif_compare.r** (2 tsv -> 2 tsv & 1 png): 
Takes as input two tsv files of motifs / regions intersection
counts. These generally originate from running region_motif_intersect.r on two sets
of different regions with the same query motif database. Based on the counts, 
region_motif_compare.r then determines the enrichment (or depletion) of certain
motifs across the two regions. This is done by a correcting for the size and gc
content of the region, and applying a Poisson test to the counts. 
Then, region_motif_compare.r outputs the most significant enriched or depleted
motifs as a tsv. In addition, the tool outputs a diagnostic plot containing
graphical representations of the motif counts, gc correction curves, and significant 
motifs that distinguish the two regions (selected via p value).

3. **Motif Database**: To count the intersections of motif locations in regions, 
the tools rely on databases of motif positions as compressed, indexed tabix files.


### Pipeline Flowchart
Green: motif enrichment analysis   
Red: motif database preparation

![Region Motif Comparison Pipeline](/region-motif-compare/doc/pipeline.png)

## Installation
Directions for installing the region-motif-compare tools into a personal computer
and a local Galaxy instance. This requires Galaxy administrative priviledges.

1. Follow the online directions to install a local instance of Galaxy (getgalaxy.org).
Optionally, follow the directions to install Refinery (refinery-platform.readthedocs.org).

2. Login as an administrator and access the Admin control panel. The tools can
be found in the Galaxy main tool shed under the repository: region_motif_enrichment
by clicking "Search and browse tool sheds" and searching the package name.
Install the tools via the Galaxy tool installer UI.

3. The tool suite requires the installation of motif databases, via a Galaxy
data manager. Galaxy should automatically inform you of this dependency and 
install the data manager. For information on how to run the data manager and
how to download the provided motif databases, see region-motif-data-manager LINK.

4. Install R and the Bioconductor R package Rsamtools for dealing with tabix files
    ```
    $ R
    > source("http://bioconductor.org/biocLite.R")
    > biocLite("Rsamtools")
    ````

5. Alternatively, you can download the tools here and install them manually
into a local Galaxy instance.

## Running the Tools
### Running from Galaxy
1. To run the tools as workflows, import the .ga workflows included in the github
via the Galaxy workflow user interface. Then, upload and select two input BED files.

2. To run the tools individually, select the tool from the tools toolbar, provide
a BED file (Region Motif Intersect) or two tsv files (Region Motif Compare), and
select a query database from the dropdown menu.

### Running from Refinery
1. Import the .ga workflows into a local Galaxy instance. These workflows have
already been annotated for Refinery.

2. Add the local Galaxy instance to the Refinery installation.
    ````
    python manage.py create_workflowengine <instance_id> "<group_name>"
    ````

3. Import the Galaxy workflows into Refinery.
    ````
    python manage.py import_workflows
    ````
4. Run the tools from the Refinery user interface.

### Running as Command Line Tools
You can also run the tools from the command line, an example of which is shown below.
More information is found in the headers of the r source files.
````
cd ~/galaxy-dist/tools/my_tools
R --slave --vanilla -f region_motif_intersect.r --args <motif_db_bgz> <motif_db_tbi> <input_bed> <output_tsv>
R --slave --vanilla -f region_motif_compare.r --args ~/galaxy-dist <motif_pwms_meme> <path_to_region1_counts> <path_to_region2_counts> <enriched_motifs_output_tsv> <depleted_motifs_output_tsv> <plots_png>
````

## Interpreting Results
### Motif Database and Result Notation
TF motif positions for hg19 and mm9 were curated from three databases:  
ENCODE TF motif database "Pouya" (http://compbio.mit.edu/encode-motifs/)  
JASPAR database "Jaspar" (http://jaspar.genereg.net/)  
DNA binding specificities of human transciption factors "Jolma" (http://www.ncbi.nlm.nih.gov/pubmed/23332764)  

For ENCODE TF motifs, the genomic locations were taken straight from the database.
In addition, position weight matrices (pwms) were obtained by averaging the 
sites in the genome for a motif. These are labeled with "\_8mer\_". 
Fake motifs were also generated, by shuffling the pwms of actual motifs and 
mapping to the genome and are labeled with "_8mer_C".

For JASPAR and Jolma motifs, mast was run to determine genomic locations from the
provided pwms. The motif alignmment thresholds were set to the top 5k, 20k, 100k, and
250k sites and the redundant maps removed with the top 30k sites have the same score. 
These are labeled with "_t5000" and likewise.


## Motif Tabix File Creation
Starting with a BED file of motif positions (minimal chr, start, end), follow 
below to generate a tabix file that can be placed in /galaxy-dist/tool-data/motifs/.

1. Download Tabix (http://sourceforge.net/projects/samtools/files/tabix/) and install.
Add `tabix` and `bgzip` binaries to your file path.
    ````
tar -xvjf tabix-0.2.6.tar.bz2
cd tabix-0.2.6
make
    ````

2. Construct bgzip files and index files.
    ````
cd ~/galaxy-dist/tools/my_tools/region_motif/db
(grep ^"#" jaspar_motifs.bed; grep -v ^"#" jaspar_motifs.bed | sort -k1,1 -k2,2n) | bgzip > jaspa_motifs.bed.bgz
tabix -p bed jaspar_motifs.bed.bgz   # this generates jaspar_motifs.bed.bgz.tbi
    ````

3. Move the bgzip and tbi files to /galaxy-dist/tool-data/motifs/. Additionally,
edit /galaxy-dist/tool-data/toolshed.g2.bx.psu.edu/repos/jeremyjliu/region_motif_data_manager/<revision number>/motif_database.loc to reflect this change.
