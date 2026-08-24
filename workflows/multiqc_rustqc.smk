if paired == 'paired':
    rule multiqc:
        input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                versions = "LOGS/versions.txt"
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC")
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
        conda:  "rustqc.yaml"
        container: "oras://jfallmann/monsda:rustqc"
        threads: 1
        params: qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); mkdir -p $OUT; SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; for d in QC/{wildcards.combo}/{wildcards.condition} QC/$(echo {wildcards.combo}|sed 's/rustqc/fastqc/g')/{wildcards.condition}; do case \" $SCAN \" in *\" $d \"*) continue;; esac; if [ -d \"$d\" ]; then SCAN=\"$SCAN $d\"; fi; done; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"

else:
    rule multiqc:
        input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                versions = "LOGS/versions.txt"
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC")
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
        conda:  "rustqc.yaml"
        container: "oras://jfallmann/monsda:rustqc"
        threads: 1
        params: qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); mkdir -p $OUT; SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; for d in QC/{wildcards.combo}/{wildcards.condition} QC/$(echo {wildcards.combo}|sed 's/rustqc/fastqc/g')/{wildcards.condition}; do case \" $SCAN \" in *\" $d \"*) continue;; esac; if [ -d \"$d\" ]; then SCAN=\"$SCAN $d\"; fi; done; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"
