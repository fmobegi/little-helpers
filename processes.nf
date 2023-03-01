 // Estimate quality of raw reads
process fastqc {
    conda '/home/fshimmunology/miniconda3/envs/fastqc_env'
    tag "$name"
    publishDir "${params.outdir}/01-reports/01-fastqc", mode: 'copy'

    input:
        tuple val(name), path(raw_reads)

    output:
        tuple val(name), path("*.html"), emit: html
        tuple val(name), path("*.zip") , emit: zip

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        fastqc -t 4 -f fastq $raw_reads
        """
}

// Convert FASTQ to FASTA using BBMap
process fq2fa {
    tag "$name"
    publishDir "${params.outdir}/02-fastq2fasta", mode: 'copy'

    input:
        tuple val(name), path(raw_reads)

    output:
        tuple val(name), path("*.fasta"), emit: fasta

    script:
        """
        reformat.sh in=${raw_reads} out=${name}.fasta
        """
}

// Run Kraken1 and Kraken2
process kraken {
    tag "${name.toString()}"
    publishDir "${params.outdir}/03-kraken", mode: 'copy'
    cache true

    input:
        tuple val(name), path(fasta)
        //path( kraken1_db )
                val(stage_at)

    output:
        tuple val(name), path("*.kraken1.report"), emit: reportk1
        tuple val(name), path("*.kraken2.report"), emit: reportk2
        tuple val(name), path("*.kraken2.txt")

    when:
        task.ext.when == null || task.ext.when

    script:

        """
        kraken --threads 4 \\
                --preload \\
                --db ${params.kraken1db} ${fasta} | \\
                kraken-report --db ${params.kraken1db} - \\
                > ${name}.${stage_at}.kraken1.report

        kraken2 --db ${params.kraken2db} \\
                --use-names  \\
                --threads ${params.num_threads} \\
                --output  ${name}.${stage_at}.kraken2.txt \\
                --report  ${name}.${stage_at}.kraken2.report  \\
                --memory-mapping ${fasta}
        """
}

// Determine the Kingdom: TAXA
process determine_kingdom {
    tag "${name.toString()}"
    publishDir "${params.outdir}/04-taxa_predict", mode: 'copy'
    cache true

    input:
        tuple val(name), path(reportk1)
        tuple val(name), path(reportk2)

    output:
        tuple val(name), path("*.kingdom.txt"), emit: taxa
    when:
        task.ext.when == null || task.ext.when

    shell:
        """
    awk '\$1>=30 && (\$4=="K" || \$4=="D") {for(i=6; i<=NF; i++) printf(" %s", \$i)}' !{reportk1} > file1.txt
    awk '\$1>=30 && (\$4=="K" || \$4=="D") {for(i=6; i<=NF; i++) printf(" %s", \$i)}' !{reportk2} > file2.txt

    if grep -q "Bacteria" file1.txt && grep -q "Bacteria" file2.txt; then
        echo "Bacteria" > kingdom.txt
    elif grep -q "Fungi" file1.txt || grep -q "Fungi" file2.txt \\
    || grep -q "Eukaryota" file1.txt || grep -q "Eukaryota" file2.txt; then
        echo "Fungi" > kingdom.txt
    else
        echo "Unknown" > kingdom.txt
    fi

    mv kingdom.txt !{name}.kingdom.txt

    rm file1.txt file2.txt
        """
}



// Run BLASTN alignment to reference
process blastn_align {
    tag "${name.toString()}"
    publishDir "${params.outdir}/05-blastn", mode: 'copy'
    cache true

    input:
        tuple val(name), path(fasta)
        val( database_name )
        val( stage_at )

    output:
        tuple val(name), path("*.txt"),    emit: text
        tuple val(name), path("*.parsed"), emit: parsed

    when:
        task.ext.when == null || task.ext.when

    script:
        """
            export BLASTDB=/data/references/blast_db

        blastn \\
            ${params.blastoptions} \\
            -num_threads $params.num_threads \\
            -outfmt ${params.outfmt} \\
            -out ${name}.${stage_at}_blastn.txt \\
            -db $database_name \\
            -query $fasta

        awk '\$3>=90 && \$6/\$4>="0.75"  {print \$0}' ${name}.${stage_at}_blastn.txt| \\
        sort -nrk3 > "${name}.${stage_at}_blastn.parsed"
        """
}

// Run MetaPhLAn4
process metaphlan4 {

conda '/home/fshimmunology/miniconda3/envs/metaphlan4'
    label 'metaphlan4'
        tag {name}
    publishDir "${params.outdir}/06-metaphlan", mode: 'copy'

    input:
        tuple val(name), path(raw_reads)
                path(metaphlandb)

    output:
                tuple val(name), path("*.txt"), emit: profile

    """
        metaphlan ${raw_reads} \\
                --nproc ${params.num_threads} \\
                --input_type fastq \\
                --bowtie2db ${metaphlandb} \\
                -o ${name}.metaphlan4.txt
    """
}

// Process Run Canu2
process canu2 {
        conda '/home/fshimmunology/miniconda3/envs/canu2.2_env'
    tag {name}
    publishDir "${params.outdir}/07-canu2_assemblies/01-unpolished", mode: 'copy'
    label 'canu2'
    errorStrategy 'ignore'

        memory '30 GB'

    input:
        tuple val(name), path(raw_reads)

    output:
                tuple val(name), path("*.fasta"), emit: assembly
                tuple val(name), path("*.fasta.gz")
                tuple val(name), path("*.report")

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
                canu \
                        -p ${name} \
                        -d ${name} \
                        genomeSize=0.4m \
                        minInputCoverage=0.2 \
                        stopOnLowCoverage=0.2 \
                        -fast \
                        -nanopore ${raw_reads}

                cp ${name}/*contigs.fasta ${name}.contigs.fasta
                cp ${name}/*correctedReads.fasta.gz ${name}.correctedReads.fasta.gz
                cp ${name}/*.report ${name}.canu.report
        """
}

// polishing step 1 Racon
process racon {
        conda '/home/fshimmunology/miniconda3/envs/racon_env'
    tag {name}
    publishDir "${params.outdir}/07-canu2_assemblies/02-racon-polish", mode: 'copy'
        errorStrategy 'ignore'

    input:
                tuple val(name), path(raw_reads)
                tuple val(name), path(canu_contig)

    output:
                tuple val(name), path("*.fasta"), emit: racon

    """
    minimap2 \
    ${canu_contig} \
    ${raw_reads} > minimap.racon.paf

    racon -m 8 -x -6 -g -8 -w 500 -t 14 \
                --no-trimming \
                ${raw_reads} \
                minimap.racon.paf \
                ${canu_contig} > ${name}.contigs.racon.fasta
    """
}

// polishing step 2 Medaka
process medaka {
    conda '/home/fshimmunology/miniconda3/envs/medaka_env'
    tag {name}
    publishDir "${params.outdir}/07-canu2_assemblies/03-medaka-polish", mode: 'copy'
        errorStrategy 'ignore'

    input:
        tuple val(name), path(raw_reads)
                tuple val(name), path(racon_contig)

    output:
                tuple val(name), path("*.medaka.fasta"), emit: medaka

    """
    medaka_consensus \
                -d ${racon_contig} \
                -i ${raw_reads} \
                -o ${name}_medaka_output \
                -m r941_min_sup_g507

    seqkit sort -lr ${name}_medaka_output/consensus.fasta > ${name}.fasta
    seqkit replace -p '.+' -r '${name}_ctg_{nr}' --nr-width 2 ${name}.fasta > ${name}.contigs.racon.medaka.fasta
    """
}

// Process MultiQC report
process multiqc {
    tag "MultiQC"
    publishDir "${params.outdir}/01-reports/02-multiqc", mode: 'copy'
    label 'multiqc'
    errorStrategy 'ignore'

    input:
        path multiqc_files

    output:
        path "*multiqc_report.html", emit: report
        path "*_data"              , emit: data
        path "*_plots"             , optional:true, emit: plots

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        multiqc -f $args .
        """

}

// Render Rmarkdown script
process render_rmarkdown {
    tag {"Rmarkdown"}
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(project_dir)
        path(rmd)

    output:
        path("*.html")
        path("*.tsv")
        path("*.xlsx")

    script:
        output_f = "${name}_raw_read_accuracy.html"

        """
        #!/usr/bin/env Rscript -e
        ## rm(list=ls())

        rmarkdown::render(
            'render_id_the_bug_report.Rmd',
            outfile = $output_f,
            params = list(locus = $db_name),
            envir = parent.frame()
            )

        """
}
