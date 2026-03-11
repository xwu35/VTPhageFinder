#------------ SET UP THE DIRECTORIES
dir = dict()
dir["output"] = dict()

# WORKFLOW DIRs
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")
dir["db"] = os.path.join(workflow.basedir, "..", "db")

# OUTPUT DIRs
dir["output"]["base"] = RESULTS_DIR
dir["output"]["reads_processing"] = os.path.join(dir["output"]["base"], "reads_processing")
dir["output"]["assembly"] = os.path.join(dir["output"]["base"], "assembly")
dir["output"]["check_contig_contamination"] = os.path.join(dir["output"]["base"], "check_contig_contamination")
dir["output"]["intermediate"] = os.path.join(dir["output"]["reads_processing"], "intermediate")

#------------ SET UP THE OUTPUT
host_removed_reads=[]
no_host_contig=[]
for sample in SAMPLE:
    host_removed_reads.append(os.path.join(dir["output"]["reads_processing"], "filtered_reads", sample + "_R1.unmapped.fastq.gz"))
    host_removed_reads.append(os.path.join(dir["output"]["reads_processing"], "filtered_reads", sample + "_R2.unmapped.fastq.gz"))
    host_removed_reads.append(os.path.join(dir["output"]["reads_processing"], "filtered_reads", sample + ".unmapped.singletons.fastq.gz"))
    no_host_contig.append(os.path.join(dir["output"]["check_contig_contamination"], "blastn_phix", sample + "_blastn_phix.out"))
    no_host_contig.append(os.path.join(dir["output"]["check_contig_contamination"], "blastn_human", sample + "_blastn_human_ANI95_AF50.out"))
    no_host_contig.append(os.path.join(dir["output"]["check_contig_contamination"], "no_host_contig_sequences", sample + "_contigs_1kb_no-host.fasta"))

# checkv scripts
scripts_input = [
    os.path.join(dir["scripts"], "anicalc.py"),
    os.path.join(dir["scripts"], "aniclust.py")
]

# raw reads: read counts and fastqc
fastqc_input = [
    os.path.join(dir["output"]["reads_processing"], "reads_statistics", "raw_reads", "R1_stats.tsv"),
    os.path.join(dir["output"]["reads_processing"], "reads_statistics", "raw_reads", "R2_stats.tsv"),
    os.path.join(dir["output"]["reads_processing"], "fastqc", "multiqc_raw_reads", "multiqc_report.html")
]

# preprocess: removed adaptoers, phix, human and host contmaination
preprocess_input = scripts_input + [os.path.join(dir["output"]["reads_processing"], "fastqc", "multiqc_after_trimmomatic", "multiqc_report.html")] + host_removed_reads + [
    os.path.join(dir["output"]["reads_processing"], "reads_statistics", "number_of_reads_removed_at_each_step.txt"),
    os.path.join(dir["output"]["reads_processing"], "reads_statistics", "reads_composition_barplot.svg")
]

# assembly and removed host contigs
assembly_input = [
    os.path.join(dir["output"]["assembly"], "contig_evaluation", "combined_assembly_stats.txt"),
    os.path.join(dir["output"]["assembly"], "all_contigs_1kb.fasta")
    ] + no_host_contig

#---------- DOWNLOAD SCRIPTS
localrules:
    download_checkV_scripts

rule download_checkV_scripts:
    """
    download scripts for calculating ANI and clustering
    """
    output:
        anicalc=os.path.join(dir["scripts"], "anicalc.py"),
        aniclust=os.path.join(dir["scripts"], "aniclust.py")
    params:
        dir=dir["scripts"]
    shell:
        """
        # download the two scripts from checkV which cannot be installed via conda
        # only download them if the files don't exsit
        # script for calculating ANI and AF
        if [[ ! -e {output.anicalc} ]]; then 
            wget -P {params.dir} https://bitbucket.org/berkeleylab/checkv/raw/3f185b5841e8c109848cd0b001df7117fe795c50/scripts/anicalc.py && chmod +x {output.anicalc}
        else
            # in case the files were downloaded, but not executable 
            chmod +x {output.anicalc}
        fi

        # script for clustering contigs into vOTUs
        if [[ ! -e {output.aniclust} ]]; then 
            wget -P {params.dir} https://bitbucket.org/berkeleylab/checkv/raw/3f185b5841e8c109848cd0b001df7117fe795c50/scripts/aniclust.py && chmod +x {output.aniclust}
        else
            # in case the files were downloaded, but not executable
            chmod +x {output.aniclust}
        fi
        """
