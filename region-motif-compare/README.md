# Region-Motif-Compare Tools
Version 1.1 Released 2014  
Park Laboratory  
Center for Biomedical Informatics  
Harvard University  

Contact  
Jeremy Liu (jeremy.liu@yale.edu)  
Nils Gehlenborg (nils@hms.harvard.edu)

## Overview
### Structure
The tool suite consists of:

1. Two Rscripts: region_motif_compare.r and region_motif_intersect.r
2. Two Xml Files: region_motif_compare.xml and region_motif_intersect.xml
3. Motif Database Directory: region_motif_db
4. Dependency Library Directory: region_motif_lib
5. Galaxy Workflows: Files with suffix ".ga" that can be imported into the local
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

3. **region_motif_db**: Contains motif positions as compressed, indexed tabix files.

4. **region_motif_lib**: Contains dependencies (i.e. plotting.r) for region_motif_compare.r

## Installation
Directions for installing the region-motif-compare tools into a personal computer
and a local Galaxy instance.

1. Follow the online directions to install a local instance of Galaxy (getgalaxy.org).
Optionally, follow the directions to install Refinery (refinery-platform.readthedocs.org)

2. Clone the github repository to your local computer
    ````
    git clone https://github.com/parklab/refinery-galaxy-tools.git
    cd refinery-galaxy-tools/region-motif-compare
    ````

3. Make a directory for the tools in Galaxy instance. This serves as a category
for the tool in the tools sidebar. You can also place the tools in an existing
or alternatively named directory, but remember to update tool_conf.xml to reflect this.
    ````
    cd ~/galaxy-dist/tools/
    mkdir my_tools
    cd my_tools
    ````

4. Copy over ".r" and ".xml" files, as well as `region_motif_db` and `region_motif_lib`
    ````
    cd refinery-galaxy-tools/region-motif-compare
    cp *.r ~/galaxy-dist/tools/my_tools
    cp *.xml ~/galaxy-dist/tools/my_tools
    cp -r region_motif_db ~/galaxy-dist/tools/my_tools
    cp -r region_motif_lib ~/galaxy-dist/tools/my_tools
    ````

5. Edit `~/galaxy-dist/tool_conf.xml` to reflect the addition of the new tools.
Add the following lines within the `<toolbox>` tags. If in Step 3 you copied
the tools to a different directory than `my_tools`, edit the code snippet
to reflect the correct path name.
    ````
    <section id="mTools" name="My Tools">  
        <tool file="my_tools/region_motif_intersect.xml" />  
        <tool file="my_tools/region_motif_compare.xml" />  
    </section>
    ````

6. Download the motif databases and place them into `region_motif_db`
    ````
    cd ~/galaxy-dist/tools/my_tools/region_motif_db
    wget ????/pouya_motifs.bed.bgz
    wget ????/pouya_motifs.bed.bgz.tbi
    wget ????/jaspar_jolma_motifs.bed.bgz
    wget ????/jaspar_jolma_motifs.bed.bgz.tbi
    wget ????/mm9_motifs.bed.bgz
    wget ????/mm9_motifs.bed.bgz.tbi
    ````

7. Install the Bioconductor R package Rsamtools for dealing with tabix files
    ```
    $ R
    > source("http://bioconductor.org/biocLite.R")
    > biocLite("Rsamtools")
    ````

8. If in Step 3 you copied the tools to an existing directory or an alternatively
named directory, you must edit the following file paths.  
    In `region_motif_intersect.r` and `region_motif_compare.r` edit `commonDir`:  
    ````
    # Replace this line
    commonDir = concat(workingDir, "/tools/my_tools")
    # With this edited line
    commonDir = concat(workingDir, "<relative_path_from_galaxy_root>/<tool_directory>")
    ````
    In addition, edit `region_motif_intersect.xml` and `region_motif_compare.xml` to
    reflect the path of the tools relative to the galaxy root directory.
    ````
    <command interpreter="bash">
        /usr/bin/R --slave --vanilla -f $GALAXY_ROOT_DIR/<path_to_tools>/region_motif_intersect.r --args $GALAXY_ROOT_DIR $db_type $in_bed $out_tab
    </command>
    ````
    ````
    <command interpreter="bash">
        /usr/bin/R --slave --vanilla -f $GALAXY_ROOT_DIR/<path_to_tools>/region_motif_compare.r --args $GALAXY_ROOT_DIR $db_type $in_tab_1 $in_tab_2 $out_enriched $out_depleted $out_plots
    </command>
    ````

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
R --slave --vanilla -f region_motif_intersect.r --args ~/galaxy-dist p <path_to_bed_file> <path_to_output_tsv>
R --slave --vanilla -f region_motif_compare.r --args ~/galaxy-dist p <path_to_region1_counts> <path_to_region2_counts> <enriched_motifs_output_tsv> <depleted_motifs_output_tsv> <plots_png>
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
below to generate a tabix file that can be placed in `region_motif_db` and
used by the tools. 

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

3. Add the path to `jaspar_motifs.bed.bgz` to the selection options for the variable
`motifDB` in `region_motif_intersect.r` and `region_motif_compare.r`. To enable
the new database in Galaxy, you will have to edit the xml files for both tools.
