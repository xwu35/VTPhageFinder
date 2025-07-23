rule make_blastdb_phix:
    """create blast database for phiX genome"""
    input:
        os.path.join("resources", "phix_genome", "phix.fasta")
    output:
        multiext(os.path.join("resources", "phix_genome", "blastDB", "phix_nucleotide_db"), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto")
    params: 
        prefix=os.path.join("resources", "phix_genome", "blastDB", "phix_nucleotide_db"),
        dbtype=config["blastn"]["dbtype"]
    threads: 
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "alignment.yml")
    shell:
        """
        makeblastdb -in {input} -dbtype {params.dbtype} -out {params.prefix}
        """

rule blastn_phix:
    """check if there contigs still aligned to the phiX genome (see if the reads removal was sufficient)"""
    input:
        blastdb=multiext(os.path.join("resources", "phix_genome", "blastDB", "phix_nucleotide_db"), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"),
        contigs_1kb=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        blastn_phix_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_phix", "{sample}_blastn_phix.out")
    params:
        db_path=os.path.join("resources", "phix_genome", "blastDB", "phix_nucleotide_db"),
        evalue=config["blastn"]["evalue"],
        max_target_seqs=config["blastn"]["max_target_seqs"],
        outfmt=config["blastn"]["outfmt"]
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "alignment.yml")
    shell:
        """
        blastn -db {params.db_path} \
            -query {input.contigs_1kb} \
            -out {output.blastn_phix_output} \
            -outfmt {params.outfmt} \
            -max_target_seqs {params.max_target_seqs} \
            -num_threads {threads} \
            -evalue {params.evalue} 
        """

rule make_blastdb_human:
    """create blast database for human genome"""
    input:
        os.path.join("resources", "human_genome", "GRCh38.fasta")
    output:
        multiext(os.path.join("resources", "human_genome", "blastDB", "GRCh38_nucleotide_db"), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto")
    params: 
        prefix=os.path.join("resources", "human_genome", "blastDB", "GRCh38_nucleotide_db"),
        dbtype=config["blastn"]["dbtype"]
    threads: 
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "alignment.yml")
    shell:
        """
        makeblastdb -in {input} -dbtype {params.dbtype} -out {params.prefix}
        """

rule blastn_human:
    """check if there are contigs still aligned to the human genome (see if the reads removal was sufficient)"""
    input:
        blastdb=multiext(os.path.join("resources", "human_genome", "blastDB", "GRCh38_nucleotide_db"), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"),
        contigs_1kb=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        blastn_human_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_human", "{sample}_blastn_human.out")
    params:
        db_path=os.path.join("resources", "human_genome", "blastDB", "GRCh38_nucleotide_db"),
        evalue=config["blastn"]["evalue"],
        max_target_seqs=config["blastn"]["max_target_seqs"],
        outfmt=config["blastn"]["outfmt"]
    threads:
        config["resources"]["med_cpu"]
    conda:
        os.path.join(dir["env"], "alignment.yml")
    shell:
        """
        blastn -db {params.db_path} \
            -query {input.contigs_1kb} \
            -out {output.blastn_human_output} \
            -outfmt {params.outfmt} \
            -max_target_seqs {params.max_target_seqs} \
            -num_threads {threads} \
            -evalue {params.evalue} 
        """

rule blastn_human_ani_calculation:
    """
    calculate ANI from blastn
    """
    input:
        blastn_human_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_human", "{sample}_blastn_human.out")
    output:
        ani_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_human", "{sample}_blastn_human_anicalc.out")
    params:
        anicalc=os.path.join(dir["scripts"], "anicalc.py")
    conda:
        os.path.join(dir["env"], "biopython.yml")
    shell:
        """
        # if blastn output is not empty, then calculate ANI
        if [[ -s {input.blastn_human_output} ]]; then
            {params.anicalc} -i {input.blastn_human_output} -o {output.ani_output}
        else
            echo "No contigs aligned to the human genome. BlastN result is empty."
            touch {output.ani_output}
        fi
        """

rule blastn_human_contig_selection:
    """
    select the contigs aligned to the human genome with ANI >= 95% and AF >=50%
    """
    input:
        ani_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_human", "{sample}_blastn_human_anicalc.out")
    output:
        ani_selected=os.path.join(dir["output"]["check_contig_contamination"], "blastn_human", "{sample}_blastn_human_ANI95_AF50.out")
    shell:
        """
        # if ani_output is not empty, then select contigs based on ANI and AF
        if [[ -s {input.ani_output} ]]; then
            awk '{{OFS="\t"}}{{FS="\t"}}{{ if (NR==1) {{print $1, $2, $4, $5}} else if ($4>=95 && $5>=50) {{print $1, $2, $4, $5}}}}' {input.ani_output} > {output.ani_selected}
        else
            echo "No contigs aligned to the human genome. ANI File is empty"
            touch {output.ani_selected}
        fi
        """

rule make_blastdb_host:
    """create blast database for host genome"""
    input:
        ref=config["host_genome"]
    output:
        multiext(os.path.join(os.path.dirname(config["host_genome"]), "blastDB", os.path.basename(config["host_genome"])), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto")
    params: 
        # get the path to the host fasta file (resources/host/host.fasta -- get resources/host) and then add blastDB and host.fasta as the prefix
        # eventually: resource/host/blastDB/host.fasta 
        prefix=os.path.join(os.path.dirname(config["host_genome"]), "blastDB", os.path.basename(config["host_genome"])),
        dbtype=config["blastn"]["dbtype"]
    threads: 
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "alignment.yml")
    shell:
        """
        makeblastdb -in {input.ref} -dbtype {params.dbtype} -out {params.prefix}
        """

rule blastn_host:
    """blast against the host genome"""
    input:
        blastdb=multiext(os.path.join(os.path.dirname(config["host_genome"]), "blastDB", os.path.basename(config["host_genome"])), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"),
        contigs_1kb=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        blastn_host_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_host", "{sample}_blastn_host.out")
    threads:
        config["resources"]["med_cpu"]
    params:
        db_path=os.path.join(os.path.dirname(config["host_genome"]), "blastDB", os.path.basename(config["host_genome"])),
        evalue=config["blastn"]["evalue"],
        max_target_seqs=config["blastn"]["max_target_seqs"],
        outfmt=config["blastn"]["outfmt"]
    conda:
        os.path.join(dir["env"], "alignment.yml")
    shell:
        """
        blastn -db {params.db_path} \
            -query {input.contigs_1kb} \
            -out {output.blastn_host_output} \
            -outfmt {params.outfmt} \
            -max_target_seqs {params.max_target_seqs} \
            -num_threads {threads} \
            -evalue {params.evalue} 
        """

rule blastn_host_ani_calculation:
    """
    caculate ANI from blastn
    """
    input:
        blastn_host_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_host", "{sample}_blastn_host.out")
    output:
        ani_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_host", "{sample}_blastn_host_anicalc.out")
    params:
        anicalc=os.path.join(dir["scripts"], "anicalc.py")
    conda:
        os.path.join(dir["env"], "biopython.yml")
    shell:
        """
        # if blastn output is not empty, then calculate ANI
        if [[ -s {input.blastn_host_output} ]]; then
            {params.anicalc} -i {input.blastn_host_output} -o {output.ani_output}
        else
            echo "No contigs aligned to the host genome. BlastN result is empty."
            touch {output.ani_output}
        fi
        """

rule blastn_host_contig_selection:
    """
    select the contigs aligned to the host genome with ANI >= 95% and AF >=85%
    """
    input:
        ani_output=os.path.join(dir["output"]["check_contig_contamination"], "blastn_host", "{sample}_blastn_host_anicalc.out")
    output:
        ani_selected=os.path.join(dir["output"]["check_contig_contamination"], "blastn_host", "{sample}_blastn_host_anicalc_selected.out")
    shell:
        """
        # if ani_output is not empty, then select contigs based on ANI and AF
        if [[ -s {input.ani_output} ]]; then
            awk '{{OFS="\t"}}{{FS="\t"}}{{ if (NR==1) {{print $1, $2, $4, $5}} else if ($4>=95 && $5>=85) {{print $1, $2, $4, $5}}}}' {input.ani_output} > {output.ani_selected}
        else
            echo "No contigs aligned to the host genome. ANI File is empty"
            touch {output.ani_selected}
        fi
        """

rule get_host_contigs_list:
    """
    extract the contig names aligned to the host genome with ANI >= 95% and AF >=85%
    """
    input:
        ani_selected=os.path.join(dir["output"]["check_contig_contamination"], "blastn_host", "{sample}_blastn_host_anicalc_selected.out")
    output:
        host_contig_id=os.path.join(dir["output"]["check_contig_contamination"], "blastn_host", "{sample}_blastn_host_anicalc_selected_contigs_id.txt")
    shell:
        """
        # if ani_output is not empty, then get the contig names aligned to host genome with ANI >= 95% and AF >=85%
        if [[ -s {input.ani_selected} ]]; then
            cat {input.ani_selected} | tail -n +2 | cut -f1 > {output.host_contig_id}
        else
            echo "No contigs aligned to the host genome with ANI >= 95% and AF >=85%"
            touch {output.host_contig_id}
        fi
        """

rule extract_no_host_contig_sequences:
    """
    remove contigs aligned to the host genome with ANI >= 95% and AF >= 85% from the contig file
    """
    input:
        contigs_1kb=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta"),
        host_contig_id=os.path.join(dir["output"]["check_contig_contamination"], "blastn_host", "{sample}_blastn_host_anicalc_selected_contigs_id.txt")
    output:
        no_host_contigs=os.path.join(dir["output"]["check_contig_contamination"], "no_host_contig_sequences", "{sample}_contigs_1kb_no-host.fasta")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # if host_contig_id is not empty, then extract the non-host contigs
        if [[ -s {input.host_contig_id} ]]; then
            seqkit grep -v -i -f {input.host_contig_id} {input.contigs_1kb} > {output.no_host_contigs}
        else
            echo "No contigs aligned to the host genome with ANI >= 95% and AF >=85%. All contigs >= 1kb are considered non-host contigs."
            cp {input.contigs_1kb} {output.no_host_contigs}
        fi
        """

