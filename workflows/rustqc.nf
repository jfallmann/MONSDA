RUSTQCENV = get_always('POSTQCENV')
RUSTQCBIN = get_always('POSTQCBIN')
RUSTQCPARAMS = get_always('rustqc_params_QC') ?: ''

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
    conda "$RUSTQCENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$RUSTQCENV"
    cpus THREADS
    cache 'lenient'
    label 'big_mem'

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        "QC/${COMBO}/${CONDITION}/"+filename.replaceFirst('^results/', '')
    }

    input:
    path bam

    output:
    path "results/**", emit: rustqc_results
    path "results/*/rustqc_summary.json", emit: rustqc_json

    script:
    fn = file(bam).getSimpleName()
    anno = file("${workflow.workDir}/../${MAPANNO}")
    """
    mkdir -p results/$fn && gzip -cdfq $anno > tmp_anno.gtf && $RUSTQCBIN rna $bam --gtf tmp_anno.gtf -t ${task.cpus} $RUSTQC_PAIRED -s $RUSTQC_STRANDED --skip-dup-check -j results/$fn/rustqc_summary.json -o results/$fn $RUSTQCPARAMS && rm -f tmp_anno.gtf
    """
}

workflow RUSTQC_MAPPING{
    take: collection

    main:

    rustqc_mapped(collection)

    emit:
    qc = rustqc_mapped.out.rustqc_results
}
