rule qcall:
    input: expand("QC/Multi/{condition}/multiqc_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)))

if paired == 'paired':
    rule multiqc:
        input: expand("QC/{rawfile}_{read}_fastqc.zip", rawfile=list(SAMPLES), read=['R1','R2']),
               versions = "LOGS/versions.txt"
        output: html = report("QC/Multi/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/{condition}/tmp"),
                lst = "QC/Multi/{condition}/qclist.txt"
        log:    "LOGS/MULTIQC/multiqc/{condition}_multiqc.log"
        conda:  "qc.yaml"
        container: "oras://jfallmann/monsda:qc"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"

else:
    rule multiqc:
        input: expand("QC/{rawfile}_fastqc.zip", rawfile=list(SAMPLES)),
               versions = "LOGS/versions.txt"
        output: html = report("QC/Multi/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/{condition}/tmp"),
                lst = "QC/Multi/{condition}/qclist.txt"
        log:    "LOGS/MULTIQC/multiqc/{condition}_multiqc.log"
        conda:  "qc.yaml"
        container: "oras://jfallmann/monsda:qc"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"
