# Setup

Create a directory called `refinery` in the `tools` directory of your Galaxy installation and copy the files from the `tools` directory of this repository into it. 
Then add the following to the `toolbox` section of the `tool_conf.xml` file of your Galaxy installation.
```xml
    <label id="simulation" text="Simulation" />
    <section id="SimulationToolSet" name="Simulation Tools">
        <tool file="refinery/Bowtie.xml" />
        <tool file="refinery/bwa_mapper.xml" />
        <tool file="refinery/cut.xml" />
        <tool file="refinery/fastq_groomer.xml" />
        <tool file="refinery/FastQC.xml" />
        <tool file="refinery/filter.xml" />
        <tool file="refinery/MACS2.xml" />
        <tool file="refinery/igv_count.xml" />
        <tool file="refinery/igv_sort.xml" />
        <tool file="refinery/sam_to_bam.xml" />
        <tool file="refinery/SPP.xml" />
    </section>
 ```
Restart your Galaxy instance and import the workflows from the `workflows` directory of this repository.
