localrules:
    rename_raw_reads,
    download_phix_genome,
    download_human_genome

#######################
### rename raw reads ##
#######################

rule rename_raw_reads:
    """rename raw reads to make sure the extension works for fastqc"""
    input:
        R1=lambda wildcards: R1_MAP[wildcards.sample],
        R2=lambda wildcards: R2_MAP[wildcards.sample]
    output:
        R1=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R1.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R2.fastq.gz")
    shell:
        """
        cp {input.R1} {output.R1} 
        cp {input.R2} {output.R2} 
        """

###################
### QC raw reads ##
###################

rule fastqc_raw_reads:
    """run fastqc on raw reads"""
    input:
        R1=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R1.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R2.fastq.gz")
    output:
        R1_zip=os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads", "{sample}_R1_fastqc.zip"),
        R2_zip=os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads", "{sample}_R2_fastqc.zip")
    params:
        dir=lambda w, output: os.path.dirname(output.R1_zip)
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        fastqc {input.R1} -o {params.dir} -t {threads}
        fastqc {input.R2} -o {params.dir} -t {threads}
        """

rule multiqc_raw_reads:
    """Aggregate fastqc results from raw reads"""
    input:
        R1_zip=expand(os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads", "{sample}_R1_fastqc.zip"), sample=SAMPLE),
        R2_zip=expand(os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads", "{sample}_R2_fastqc.zip"), sample=SAMPLE)
    output:
        html=os.path.join(dir["output"]["reads_processing"], "fastqc", "multiqc_raw_reads", "multiqc_report.html")
    params:
        in_dir=directory(os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads")),
        out_dir=lambda w, output: os.path.dirname(output.html)
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        multiqc {params.in_dir} --outdir {params.out_dir} 
        """

rule get_raw_reads_status:
    """check raw reads counts"""
    input:
        R1=expand(os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R1.fastq.gz"), sample=SAMPLE),
        R2=expand(os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R2.fastq.gz"), sample=SAMPLE)
    output:
        R1_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "raw_reads", "R1_stats.tsv"),
        R2_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "raw_reads", "R2_stats.tsv")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    threads:
        config["resources"]["small_cpu"]
    shell:
        """
        seqkit stats -j {threads} -To {output.R1_stats} {input.R1}
        seqkit stats -j {threads} -To {output.R2_stats} {input.R2}
        """

##########################
### reads preprocessing ##
##########################

rule trimmomatic:
    """remove adapters and low quality reads"""
    input:
        R1=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R1.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R2.fastq.gz")
    output:
        R1P=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"),
        R1U=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R1.unpaired.fastq.gz"),
        R2P=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz"),
        R2U=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R2.unpaired.fastq.gz")
    params:
        seqType = config["trimmomatic"]["seqType"],
        phred = config["trimmomatic"]["phred"],
        adapter = config["trimmomatic"]["adapter"],
        adapter_params = config["trimmomatic"]["adapter_params"],
        post_adapter_params = config["trimmomatic"]["post_adapter_params"],
    log:
        os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}.trimmomatic.log"),
    threads:
        config["resources"]["med_cpu"]
    conda:
        os.path.join(dir["env"], "trimmomatic.yml")
    shell:
        """
        trimmomatic {params.seqType} {params.phred} -threads {threads} \
        {input.R1} {input.R2} \
        {output.R1P} {output.R1U} {output.R2P} {output.R2U} \
        ILLUMINACLIP:{params.adapter}:{params.adapter_params} \
        {params.post_adapter_params} 2>{log} # –baseout
        """

rule fastqc_trimmed_reads:
    """run fastqc on trimmed reads after trimmomatic to check if adapters were removed"""
    input:
        R1P=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"),
        R2P=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz")
    output:
        R1_zip=os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic", "{sample}_R1.trimmomatic_fastqc.zip"),
        R2_zip=os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic", "{sample}_R2.trimmomatic_fastqc.zip")
    params:
        dir=lambda w, output: os.path.dirname(output.R1_zip)
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        fastqc {input.R1P} -o {params.dir} -t {threads}
        fastqc {input.R2P} -o {params.dir} -t {threads}
        """

rule multiqc_trimmed_reads:
    """aggregate fastqc results from trimmed reads after trimmomatic"""
    input:
        R1_zip=expand(os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic", "{sample}_R1.trimmomatic_fastqc.zip"), sample=SAMPLE),
        R2_zip=expand(os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic", "{sample}_R2.trimmomatic_fastqc.zip"), sample=SAMPLE)
    output:
        html=os.path.join(dir["output"]["reads_processing"], "fastqc", "multiqc_after_trimmomatic", "multiqc_report.html")
    params:
        in_dir=directory(os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic")),
        out_dir=lambda w, output: os.path.dirname(output.html)
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        multiqc {params.in_dir} --outdir {params.out_dir} 
        """

rule get_reads_status_after_trimmomatic:
    """check reads counts after removing adapters and low quality reads"""
    input:
        R1=expand(os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"), sample=SAMPLE),
        R2=expand(os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz"), sample=SAMPLE)
    output:
        R1_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_trimmomatic", "R1_stats.tsv"),
        R2_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_trimmomatic", "R2_stats.tsv")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        seqkit stats -j {threads} -To {output.R1_stats} {input.R1}
        seqkit stats -j {threads} -To {output.R2_stats} {input.R2}
        """

rule download_phix_genome:
    """download the phix genome (GCF_000819615.1_ViralProj14015 is the same as NC_001422.1, the one used in FpVT paper). Check the sequence header"""
    output:
        os.path.join("resources", "phix_genome", "phix.fasta")
    log:
        os.path.join("resources", "phix_genome", "phix.log")
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
    shell:
        """
        # if the phix.fasta file doesn't exist, then download
        if [[ ! -f "{output}" ]]; then
            curl -L -o {output}.gz -s {params.url} > {log} 2>&1
            gunzip {output}.gz
        fi
        """

if config["mapper"]=="bowtie2":
    rule bowtie2_build_phix_genome:
        """build bowtie2 index for phiX genome"""
        input:
            phix_fasta=os.path.join("resources", "phix_genome", "phix.fasta")
        output:
            multiext(os.path.join("resources", "phix_genome", "phix"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params: 
            prefix=lambda w, input: os.path.splitext(input.phix_fasta)[0]
        threads: 
            config["resources"]["small_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2-build \
                --threads {threads} \
                {input.phix_fasta} \
                {params.prefix} 
            """

    rule phix_genome_mapping_bowtie2:
        """map reads back to phiX genome using bowtie2"""
        input:
            bt_index=multiext(os.path.join("resources", "phix_genome", "phix"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
            R1P=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"),
            R2P=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}.sam"))
        params:
            setting=config["bowtie2"]["extra_settings"],
            prefix=os.path.join("resources", "phix_genome", "phix")
        log:
            os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}.bowtie2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2 {params.setting} -p {threads} \
            -x {params.prefix} -1 {input.R1P} -2 {input.R2P} \
            -S {output.sam} 2> {log} 
            """

elif config["mapper"]=="minimap2":
    rule phix_genome_mapping_minimap2:
        """map reads back to phix genome using minimap2"""
        input:
            phix_fasta=os.path.join("resources", "phix_genome", "phix.fasta"),
            R1=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"),
            R2=os.path.join(dir["output"]["intermediate"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}.sam"))
        params:
            setting=config["minimap2"]["settings"]
        log:
            os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}.minimap2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "coverm.yml")
        shell:
            """
            minimap2 -ax {params.setting} \
            -t {threads} \
            --secondary=no \
            {input.phix_fasta} \
            {input.R1} \
            {input.R2} > {output.sam} 2>{log}
            """

rule phix_remove_mapped_reads:
    """keep unmapped paried reads and sort by name"""
    input:
        sam=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}.sam")
    output:
        bam=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_unmapped_sorted.bam")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # keep pairs where both reads are unmapped and sort by name
        samtools view -bS -f 12 {input.sam} | samtools sort - -n -o {output.bam} -@ {threads}
        """

rule phix_extract_unmapped_paired_reads:
    """extract the unmapped paired reads"""
    input:
        bam=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_unmapped_sorted.bam")
    output:
        R1P=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_R1.phixfilt.fastq.gz"),
        R2P=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_R2.phixfilt.fastq.gz")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools fastq \
        -1 {output.R1P} \
        -2 {output.R2P} \
        -0 /dev/null \
        -s /dev/null \
        -N {input.bam} \
        -@ {threads} # -N: add /1 or /2 to the read names
        """

rule get_reads_status_after_phix_reads_removal:
    """check reads counts after removing phix contaminated reads"""
    input:
        R1=expand(os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_R1.phixfilt.fastq.gz"), sample=SAMPLE),
        R2=expand(os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_R2.phixfilt.fastq.gz"), sample=SAMPLE)
    output:
        R1_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_phix_reads_removal", "R1_stats.tsv"),
        R2_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_phix_reads_removal", "R2_stats.tsv")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        seqkit stats -j {threads} -To {output.R1_stats} {input.R1}
        seqkit stats -j {threads} -To {output.R2_stats} {input.R2}
        """

rule download_human_genome:
    """download the human genome"""
    output:
        os.path.join("resources", "human_genome", "GRCh38.fasta")
    log:
        os.path.join("resources", "human_genome", "GRCh38.log")
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    shell:
        """
        # if the GRCh38.fasta file doesn't exist, then download
        if [[ ! -f "{output}" ]]; then
            curl -L -o {output}.gz -s {params.url} > {log} 2>&1
            gunzip {output}.gz
        fi
        """

if config["mapper"]=="bowtie2":
    rule bowtie2_build_human_genome:
        """build bowtie2 index for human genome"""
        input:
            human_fasta=os.path.join("resources", "human_genome", "GRCh38.fasta")
        output:
            multiext(os.path.join("resources", "human_genome", "GRCh38"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params: 
            prefix=lambda w, input: os.path.splitext(input.human_fasta)[0]
        threads: 
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2-build \
                --threads {threads} \
                {input.human_fasta} \
                {params.prefix} 
            """

    rule human_genome_mapping_bowtie2:
        """map reads back to human genome using bowtie2"""
        input:
            bt_index=multiext(os.path.join("resources", "human_genome", "GRCh38"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
            R1=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_R1.phixfilt.fastq.gz"),
            R2=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_R2.phixfilt.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}.sam"))
        params:
            setting=config["bowtie2"]["extra_settings"],
            prefix=os.path.join("resources", "human_genome", "GRCh38")
        log:
            os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}.bowtie2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2 {params.setting} -p {threads} \
                -x {params.prefix} \
                -1 {input.R1} -2 {input.R2} \
                -S {output.sam} 2> {log} 
            """

elif config["mapper"]=="minimap2":
    rule human_genome_mapping_minimap2:
        """map reads back to human genome using minimap2"""
        input:
            R1=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_R1.phixfilt.fastq.gz"),
            R2=os.path.join(dir["output"]["intermediate"], "phix_filtered", "{sample}_R2.phixfilt.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}.sam"))
        params:
            setting=config["minimap2"]["settings"],
            human_fasta=os.path.join("resources", "human_genome", "GRCh38.fasta")
        log:
            os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}.minimap2.log")
        threads:
            config["resources"]["med_cpu"]
        resources:
            mem_mb=config["resources"]["small_mem"]
        conda:
            os.path.join(dir["env"], "coverm.yml")
        shell:
            """
            minimap2 -ax {params.setting} \
                -t {threads} \
                --secondary=no \
                {params.human_fasta} \
                {input.R1} \
                {input.R2} > {output.sam} 2>{log}
            """

rule human_remove_mapped_reads:
    """keep unmapped paried reads and sort by name"""
    input:
        sam=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}.sam")
    output:
        bam=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_unmapped_sorted.bam")
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # keep pairs where both reads are unmapped and sort by name
        samtools view -bS -f 12 {input.sam} | samtools sort - -n -o {output.bam} -@ {threads}
        """

rule human_extract_unmapped_paired_reads:
    """extract unmapped paried reads"""
    input:
        bam=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_unmapped_sorted.bam")
    output:
        R1=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_R1.humanfilt.fastq.gz"),
        R2=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_R2.humanfilt.fastq.gz")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools fastq \
            -1 {output.R1} \
            -2 {output.R2} \
            -0 /dev/null \
            -s /dev/null \
            -N {input.bam} \
            -@ {threads}
        """

rule get_reads_status_after_human_reads_removal:
    """check read counts after removing human contaminated reads"""
    input:
        R1=expand(os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_R1.humanfilt.fastq.gz"), sample=SAMPLE),
        R2=expand(os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_R2.humanfilt.fastq.gz"), sample=SAMPLE)
    output:
        R1_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_human_reads_removal", "R1_stats.tsv"),
        R2_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_human_reads_removal", "R2_stats.tsv")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    threads:
        config["resources"]["small_cpu"]
    shell:
        """
        seqkit stats -j {threads} -To {output.R1_stats} {input.R1}
        seqkit stats -j {threads} -To {output.R2_stats} {input.R2}
        """

if config["mapper"]=="bowtie2":
    rule bowtie2_build_host_genome:
        """build bowtie2 index for host genome"""
        input:
            ref=config["host_genome"]
        output:
            multiext(config["host_genome"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params: 
            prefix=config["host_genome"]
        threads: 
            config["resources"]["small_cpu"]
        conda:
           os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2-build \
                --threads {threads} \
                {input.ref} \
                {params.prefix} 
            """

    rule host_genome_mapping_bowtie2:
        """map reads back to host genome using bowtie2"""
        input:
            bt_index=multiext(config["host_genome"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
            R1=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_R1.humanfilt.fastq.gz"),
            R2=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_R2.humanfilt.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}.sam"))
        params:
            setting=config["bowtie2"]["extra_settings"],
            prefix=config["host_genome"]
        log:
            os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}.bowtie2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2 {params.setting} -p {threads} \
                -x {params.prefix} \
                -1 {input.R1} -2 {input.R2} \
                -S {output.sam} 2> {log} 
            """
            
elif config["mapper"]=="minimap2":
    rule host_genome_mapping_minimap2:
        """map reads back to host genome using minimap2"""
        input:
            R1=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_R1.humanfilt.fastq.gz"),
            R2=os.path.join(dir["output"]["intermediate"], "human_filtered", "{sample}_R2.humanfilt.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}.sam"))
        params:
            setting=config["minimap2"]["settings"],
            host_fasta=config["host_genome"]
        log:
            os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}.minmap2.log")
        threads:
            config["resources"]["med_cpu"]
        resources:
            mem_mb=config["resources"]["small_mem"]
        conda:
            os.path.join(dir["env"], "coverm.yml")
        shell:
            """
            minimap2 -ax {params.setting} \
                -t {threads} \
                --secondary=no {params.host_fasta} \
                {input.R1} {input.R2} > {output.sam} 2> {log}
            """

rule host_separate_mapped_and_umapped_reads:
    """separate mapped and umapped reads"""
    input:
        sam=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}.sam")
    output:
        mapped=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted.bam"),
        unmapped=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_unmapped_sorted.bam")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # separate mapped and unmapped reads and sort the bam file by reference coordinates 
        samtools view -bS -F 4 {input.sam} | samtools sort - -o {output.mapped} -@ {threads}
        samtools view -bS -f 4 {input.sam} | samtools sort - -o {output.unmapped} -@ {threads}
        """

rule host_separate_prophage_and_nonphage:
    """separate prophage and nonprophage regions"""
    input:
        mapped=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted.bam")
    output:
        prophage=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted_prophage.bam"),
        nonprophage=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted_non-prophage.bam")
    params:
        bedFile=config["host_prophage_region"]
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # index the bam file
        samtools index {input.mapped} 

        # separate mapped prophage and non-prophage reads
        # {output.prophage} will contain alignments that overlap with regions in the bed file
        samtools view {input.mapped} \
            -b -h -o {output.prophage} \
            -U {output.nonprophage} \
            -L {params.bedFile} \
            -@ {threads}
        """

rule host_filter_prophage_region:
    """remove reads mapped to the prophage region using the given percent identity"""
    input:
        prophage=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted_prophage.bam")
    output:
        prophageFiltered=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted_prophage_filtered.bam")
    threads:
        config["resources"]["med_cpu"]
    params:
        pid=config["prophage_identity"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """                    
        coverm filter -b {input.prophage} \
            -o {output.prophageFiltered} \
            --inverse \
            --min-read-percent-identity {params.pid} \
            --threads {threads}
        """

rule host_filter_nonphage_region:
    """remove reads mapped to the non-prophage region using the given percent identity"""
    input:
        nonprophage=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted_non-prophage.bam")
    output:
        nonprophageFiltered=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted_non-prophage_filtered.bam")
    threads:
        config["resources"]["med_cpu"]
    params:
        pid=config["nonprophage_identity"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """                  
        coverm filter -b {input.nonprophage} \
            -o {output.nonprophageFiltered} \
            --inverse \
            --min-read-percent-identity {params.pid} \
            --threads {threads}
        """

rule host_merge_bam_files:
    """merge unmapped, prophage filtered and non-prophage filtered bam files"""
    input:
        unmapped=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_unmapped_sorted.bam"),
        prophageFiltered=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted_prophage_filtered.bam"),
        nonprophageFiltered=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_mapped_sorted_non-prophage_filtered.bam")
    output:
        kept=temp(os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_kept.bam"))
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools merge -o {output.kept} -@ {threads} \
            {input.unmapped} \
            {input.prophageFiltered} \
            {input.nonprophageFiltered} 
        """

rule host_sort_merged_bam:
    """sort the merged bam files by name"""
    input:
        kept=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_kept.bam")
    output:
        keptSorted=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_kept_sorted.bam")
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # sort by name 
        samtools sort -n {input.kept} -o {output.keptSorted}
        """

rule host_extract_kept_reads:
    """extract the kept reads"""
    input:
        keptSorted=os.path.join(dir["output"]["intermediate"], "host_filtered", "{sample}_kept_sorted.bam")
    output:
        R1=os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}_R1.unmapped.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}_R2.unmapped.fastq.gz"),
        singleton=os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}.unmapped.singletons.fastq.gz")
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools fastq \
            -1 {output.R1} -2 {output.R2} \
            -s {output.singleton} \
            -N {input.keptSorted} 
        """

rule get_final_reads_status:
    """check read counts after removing host contamination"""
    input:
        R1=expand(os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}_R1.unmapped.fastq.gz"), sample=SAMPLE),
        R2=expand(os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}_R2.unmapped.fastq.gz"), sample=SAMPLE),
        singleton=expand(os.path.join(dir["output"]["reads_processing"], "filtered_reads", "{sample}.unmapped.singletons.fastq.gz"), sample=SAMPLE)
    output:
        R1_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_host_reads_removal", "R1_stats.tsv"),
        R2_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_host_reads_removal", "R2_stats.tsv"),
        singleton_stats=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_host_reads_removal", "singleton_stats.tsv")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        seqkit stats -j {threads} -To {output.R1_stats} {input.R1}
        seqkit stats -j {threads} -To {output.R2_stats} {input.R2}
        seqkit stats -j {threads} -To {output.singleton_stats} {input.singleton}
        """

rule combine_all_reads_status:
    """combine reads counts and calculate number of removed reads at each step"""
    input:
        raw_reads=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "raw_reads", "R1_stats.tsv"),
        after_trimmomatic=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_trimmomatic", "R1_stats.tsv"),
        after_phix=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_phix_reads_removal", "R1_stats.tsv"),
        after_human=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_human_reads_removal", "R1_stats.tsv"),
        after_host_paired=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_host_reads_removal", "R1_stats.tsv"),
        after_host_singleton=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "after_host_reads_removal", "singleton_stats.tsv")
    output:
        table=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "number_of_reads_removed_at_each_step.txt"),
        figure=os.path.join(dir["output"]["reads_processing"], "reads_statistics", "reads_composition_barplot.svg")
    params:
        script=os.path.join(dir["scripts"], "combine_read_counts_and_plot.R")
    conda:
        os.path.join(dir["env"], "R.yml")
    script:
        "{params.script}"