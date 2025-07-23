rule spades_assembly:
    input:
        R1=os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}_R1.unmapped.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}_R2.unmapped.fastq.gz"),
        singleton=os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}.unmapped.singletons.fastq.gz")
    output:
        contigs=os.path.join(dir["output"]["assembly"], "intermediate", "{sample}", "contigs.fasta"),
        renamed=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs.fasta")
    params:
        outdir=directory(os.path.join(dir["output"]["assembly"], "intermediate", "{sample}")),
        setting=config["spades"]["setting"]
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "spades.yml")
    shell:
        """
        # the file size of the singleton fastq files are not zero since they were not created using touch 
        # so cannot use [ -s {input.singleton} ], will count the reads
        # only use singleton fastq file if the total number of lines is >= 4 (-ge 4: greater than or equal; 4 lines means 1 reads in fastq)
        if [[ $(zcat {input.singleton} | wc -l) -ge 4 ]]; then
            single="-s {input.singleton}"
        else
            single=""
        fi

        # run SPADES
        spades.py {params.setting} \
            -1 {input.R1} \
            -2 {input.R2} \
            $single -o {params.outdir} \
            -t {threads}

        # rename contigs using the sample name
        sed 's/>/>{wildcards.sample}_/' {output.contigs} > {output.renamed}
        """

rule extract_1kb_contigs:
    """keep contigs >= 1kb"""
    input:
        renamed=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs.fasta")
    output:
        contigs_1kb=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # if contig file is not empty, then select contigs >= 1kb
        if [[ -s {input.renamed} ]]; then
            seqkit seq -m 1000 {input.renamed} > {output.contigs_1kb}
        else
            # if samples have no assembled contigs, then create an empty file to avoid error
            echo "No contigs were assembled from this sample"
            touch {output.contigs_1kb}
        fi
        """

rule assembly_status:
    """Check the quality of contigs using Quast"""
    input:
        contigs_1kb=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        report=os.path.join(dir["output"]["assembly"], "contig_evaluation", "{sample}", "report.txt")
    threads: 
        config["resources"]["small_cpu"]
    params:
        dir=lambda w, output: os.path.dirname(output.report),
        length=config["quast"]["length"]
    conda:
        os.path.join(dir["env"], "quast.yml")
    shell:
        """
        # if contigs_1kb is not empty, then evaluate contig status
        if [[ -s {input.contigs_1kb} ]]; then
            quast.py {input.contigs_1kb} -m {params.length} -t {threads} -o {params.dir}
        else
            echo "No contigs >= 1kb were obtained from this sample"
            touch {output.report}
        fi
        """

rule combine_all_contigs_1kb:
    """
    combine contigs from all samples
    """
    input:
        expand(os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta"), sample=SAMPLE)
    output:
        os.path.join(dir["output"]["assembly"], "all_contigs_1kb.fasta")
    shell:
        """
        cat {input} > {output}
        """