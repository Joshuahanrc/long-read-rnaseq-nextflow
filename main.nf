#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Task 5: Long-read RNA-seq workflow
 * Adapted from Jonathan Göke's workshop workflow.
 */

params.samples    = null
params.refFa      = null
params.refGtf     = null
params.outdir     = 'task5_results'
params.bambu_mode = 'with_annotations'


process MINIMAP2_ALIGN {

    tag "${sample_id} (${protocol})"

    cpus 8
    memory '28 GB'
    maxForks 1

    input:
    tuple val(sample_id), val(protocol), path(reads)
    path refFa

    output:
    tuple val(sample_id), val(protocol),
          path("${sample_id}.sam"),
          emit: sam

    script:
    def minimap_options =
        protocol == 'directRNA' ?
        '-ax splice -uf -k14' :
        '-ax splice'

    """
    minimap2 \
        -t ${task.cpus} \
        ${minimap_options} \
        ${refFa} \
        ${reads} \
        > ${sample_id}.sam
    """
}


process SAM_TO_BAM {

    tag "${sample_id}"

    cpus 4
    memory '8 GB'

    publishDir "${params.outdir}/bam",
               mode: 'copy',
               overwrite: true

    input:
    tuple val(sample_id), val(protocol), path(reads_sam)

    output:
    tuple val(sample_id), val(protocol),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai"),
          emit: bam

    script:
    """
    samtools view \
        -@ ${task.cpus} \
        -b \
        -o ${sample_id}.unsorted.bam \
        ${reads_sam}

    samtools sort \
        -@ ${task.cpus} \
        -o ${sample_id}.sorted.bam \
        ${sample_id}.unsorted.bam

    samtools index \
        -@ ${task.cpus} \
        ${sample_id}.sorted.bam

    rm ${sample_id}.unsorted.bam
    """
}


process BAM_QC {

    tag "${sample_id}"

    cpus 2
    memory '4 GB'

    publishDir "${params.outdir}/qc",
               mode: 'copy',
               overwrite: true

    input:
    tuple val(sample_id), val(protocol),
          path(reads_bam), path(reads_bai)

    output:
    tuple val(sample_id),
          path("${sample_id}.quickcheck.txt"),
          path("${sample_id}.flagstat.txt"),
          path("${sample_id}.stats.txt"),
          emit: qc

    script:
    """
    samtools quickcheck -v ${reads_bam} \
        > ${sample_id}.quickcheck.txt 2>&1 \
        && echo "PASS" >> ${sample_id}.quickcheck.txt \
        || echo "FAIL" >> ${sample_id}.quickcheck.txt

    samtools flagstat \
        -@ ${task.cpus} \
        ${reads_bam} \
        > ${sample_id}.flagstat.txt

    samtools stats \
        -@ ${task.cpus} \
        ${reads_bam} \
        > ${sample_id}.stats.txt
    """
}


process BAMBU {

    tag "${bambu_mode}"

    cpus 2
    memory '24 GB'

    publishDir "${params.outdir}/bambu/${params.bambu_mode}",
               mode: 'copy',
               overwrite: true

    input:
    path refFa
    path refGtf
    path bam_files
    path bai_files
    val bambu_mode

    output:
    path "counts_transcript.txt"
    path "counts_gene.txt"
    path "extended_annotations.gtf"
    path "bambu_se.rds"

    script:
    if (bambu_mode == 'with_annotations') {

        bambu_call = """
        annotations <- prepareAnnotations("${refGtf}")

        se <- bambu(
            reads = bam_files,
            annotations = annotations,
            genome = "${refFa}",
            NDR = 1,
            ncore = ${task.cpus}
        )
        """

    } else if (bambu_mode == 'without_annotations') {

        bambu_call = """
        se <- bambu(
            reads = bam_files,
            genome = "${refFa}",
            NDR = 1,
            ncore = ${task.cpus}
        )
        """

    } else {
        error "bambu_mode must be with_annotations or without_annotations"
    }

    """
    #!/usr/bin/env Rscript --vanilla

    library(bambu)

    bam_files <- list.files(
        pattern = "\\\\.sorted\\\\.bam\$",
        full.names = TRUE
    )

    ${bambu_call}

    writeBambuOutput(se, path = "./")
    saveRDS(se, file = "bambu_se.rds")
    """
}


workflow {

    if (!params.samples) {
        error "Missing --samples"
    }

    if (!params.refFa) {
        error "Missing --refFa"
    }

    if (!params.refGtf) {
        error "Missing --refGtf"
    }

    if (!(params.bambu_mode in
          ['with_annotations', 'without_annotations'])) {
        error "Invalid --bambu_mode"
    }

    samples_ch = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                row.protocol,
                file(row.fastq, checkIfExists: true)
            )
        }

    ref_ch = Channel.value(
        file(params.refFa, checkIfExists: true)
    )

    gtf_ch = Channel.value(
        file(params.refGtf, checkIfExists: true)
    )

    mode_ch = Channel.value(params.bambu_mode)

    MINIMAP2_ALIGN(samples_ch, ref_ch)

    SAM_TO_BAM(MINIMAP2_ALIGN.out.sam)

    BAM_QC(SAM_TO_BAM.out.bam)

    bam_files_ch = SAM_TO_BAM.out.bam
        .map { sample_id, protocol, bam, bai -> bam }
        .collect()

    bai_files_ch = SAM_TO_BAM.out.bam
        .map { sample_id, protocol, bam, bai -> bai }
        .collect()

    BAMBU(
        ref_ch,
        gtf_ch,
        bam_files_ch,
        bai_files_ch,
        mode_ch
    )
}
