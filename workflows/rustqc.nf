QCENV = get_always('POSTQCENV')
QCBIN = get_always('POSTQCBIN')
QCPARAMS = get_always('rustqc_params_QC') ?: ''

MAPANNO = get_always('MAPPINGANNO')

// Map MONSDA strandedness to RustQC strandedness
def rustqc_stranded(stranded) {
    if (stranded == 'fr') return 'forward'
    else if (stranded == 'rf') return 'reverse'
    else return 'unstranded'
}

RUSTQC_STRANDED = rustqc_stranded(STRANDED ?: '')
RUSTQC_PAIRED = (PAIRED == 'paired') ? '-p' : ''

//RUSTQC on mapped BAMs

process rustqc_mapped{
    conda "$QCENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$QCENV"
    cpus THREADS
    cache 'lenient'
    label 'big_mem'

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("rustqc_summary.json") > 0)      "QC/${COMBO}/${CONDITION}/${file(filename).getParent().getName()}/rustqc_summary.json"
        else                                                    "QC/${COMBO}/${CONDITION}/${filename}"
    }

    input:
    path bam

    output:
    path "results/**", emit: rustqc_results
    path "results/rustqc_summary.json", emit: rustqc_json

    script:
    fn = file(bam).getSimpleName()
    anno = file("${workflow.workDir}/../${MAPANNO}")
    """
    rustqc rna $bam --gtf $anno -t ${task.cpus} $RUSTQC_PAIRED -s $RUSTQC_STRANDED --skip-dup-check -j results/rustqc_summary.json -o results/$fn $QCPARAMS
    """
}

workflow QC_MAPPING{
    take: collection

    main:

    rustqc_mapped(collection)

    emit:
    qc = rustqc_mapped.out.rustqc_results
}
