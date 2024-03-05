#Data downloaded from: http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz

# Default rule to run entire workflow, only works once inputs/outputs correctly filled in all rules
rule all:
    input: "SRR1705851_VarScan_results.vcf"

# Download raw sequence data
rule download_data:
    output: "SRR1705851.fastq.gz"
    shell:
        "wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz -O {output}"

# Create an index for the reference genome for bwa
rule index_genome_bwa:
    input: "influenza_hemagglutinin_gene.fasta"
    output:
        expand("influenza_hemagglutinin_gene.fasta.{ext}", ext=['sa', 'amb', 'ann', 'pac', 'bwt'])
    shell:
        "bwa index {input}"

# Map the raw reads to the reference genome
rule map_reads:
    input:
        genome = "influenza_hemagglutinin_gene.fasta",
        reads = "SRR1705851.fastq.gz",
        idxfile = expand("influenza_hemagglutinin_gene.fasta.{ext}", ext=['sa', 'amb', 'ann', 'pac', 'bwt'])
    output: "SRR1705851.sam"
    shell:
        "bwa mem -t 4 {input.genome} {input.reads} > {output}"

# Convert .sam to .bam file
rule samtools_import:
    input:
        samfile="SRR1705851.sam"
    output:"SRR1705851.bam"
    shell:
        """
        samtools view -S -b {input.samfile} -o {output} 
        """
        # original command with samtools v1.9
        ## samtools import {input.index} {input.samfile} {output}
        # but it gave segmentation fault error with samtools v1.10, so here we use samtools view instead

# Sort the bam alignment file
rule samtools_sort:
    input: "SRR1705851.bam"
    output: "alignment_SRR1705851_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

# Create an index for the sorted bam file
# again, this is different from the above indexing steps
rule bwa_index_sorted:
    input: "alignment_SRR1705851_sorted.bam"
    output: "alignment_SRR1705851_sorted.bam.bai"
    shell: "samtools index {input}"

# Generate pileup file with samtools, then call variants with bcftools
# From samtools doc: 'Pileup format consists of TAB-separated lines, with each line representing the pileup of reads at a single genomic position.'
rule samtools_mpileup:
    input:
        index="influenza_hemagglutinin_gene.fasta",
        sorted="alignment_SRR1705851_sorted.bam",
    output:
        "SRR1705851_VarScan_results.vcf"
    shell:
        """
        samtools mpileup -f {input.index} {input.sorted} -d 0 > my.mpileup && \
        java -jar VarScan.v2.4.0.jar mpileup2snp my.mpileup --min-var-freq 0.95 --variants --output-vcf 1 > {output}
        """
