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
5. Galaxy Workflows: Files with suffix ".ga" that can be imported into the local
Galaxy instance after installation of the tool.

### Description
1. **region_motif_intersect.r** (1 bed -> 1 tsv): 
Takes one bed file of regions as input. Then it calculates
the number of intersections of the regions and the motifs. region_motifs_intersect.r
outputs a tab separated values (tsv) file of motif names and intersection counts.

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
Add the following lines within the `<toolbox>` tags.
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

7. If in Step 3 you copied the tools to an existing directory or an alternatively
named directory, you must edit the following file paths.  
    In `region_motif_intersect.r` and `region_motif_compare.r` edit `commonDir`:  
    ````
    # Replace this line
    commonDir = concat(workingDir, "/tools/my_tools")
    # With this edited line
    commonDir = concat(workingDir, "<relative_path_from_galaxy_root>/<tool_directory>")
    ````



## Running the Tools
### Running from Galaxy and Refinery

### Running as Command Line Tools

## Motif Tabix File Creation
