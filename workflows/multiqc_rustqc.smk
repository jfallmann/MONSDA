if paired == 'paired':
    rule multiqc:
        input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                versions = "LOGS/versions.txt"
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                lst = "QC/Multi/{combo}/{condition}/qclist.txt"
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
        conda:  "rustqc.yaml"
        container: "oras://jfallmann/monsda:rustqc"
        threads: 1
        params: qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; FQ_COMBO=$(echo {wildcards.combo} | sed 's/rustqc/fastqc/g'); FQ_DIR=QC/${{FQ_COMBO}}/{wildcards.condition}; if [ -d \"${{FQ_DIR}}\" ]; then echo ${{FQ_DIR}} >> {output.tmp}; for f in ${{FQ_DIR}}/*;do ln -sfnr \"$f\" $COLLECT/ 2>> {log} || cp -f \"$f\" $COLLECT/ 2>> {log};done; fi; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"

else:
    rule multiqc:
        input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                versions = "LOGS/versions.txt"
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                lst = "QC/Multi/{combo}/{condition}/qclist.txt"
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
        conda:  "rustqc.yaml"
        container: "oras://jfallmann/monsda:rustqc"
        threads: 1
        params: qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; FQ_COMBO=$(echo {wildcards.combo} | sed 's/rustqc/fastqc/g'); FQ_DIR=QC/${{FQ_COMBO}}/{wildcards.condition}; if [ -d \"${{FQ_DIR}}\" ]; then echo ${{FQ_DIR}} >> {output.tmp}; for f in ${{FQ_DIR}}/*;do ln -sfnr \"$f\" $COLLECT/ 2>> {log} || cp -f \"$f\" $COLLECT/ 2>> {log};done; fi; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"
