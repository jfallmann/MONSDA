# Mappers that do not emit SAM (e.g. piscem writing RAD) have no sam2bam step,
# so neither the BAMs nor any QC running on them are part of the workflow.
_hasbam = rules.has('sam2bam')
_mapbams = (
    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo)
    + expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo)
) if _hasbam else []

def _bamqc(rulename, **kwargs):
    """QC targets of a BAM consuming rule, empty for non-SAM mappers."""
    if not _hasbam or not rules.has(rulename):
        return []
    return expand(getattr(rules, rulename).output.js, file=samplecond(SAMPLES, config), combo=combo, **kwargs)

if paired == 'paired':
    rule multiqc:
        input:  _bamqc('rustqc_mapped'),
                _bamqc('rustqc_uniquemapped'),
                _mapbams,
                versions = "LOGS/versions.txt"
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC")
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
        conda:  "rustqc.yaml"
        container: "oras://jfallmann/monsda:rustqc"
        threads: 1
        params: qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', ""),
                preqcenv = config.get('PREQCENV', 'fastqc')
        shell:  "OUT=$(dirname {output.html}); mkdir -p $OUT; SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; for d in QC/{wildcards.combo}/{wildcards.condition} QC/$(echo {wildcards.combo}|sed 's/rustqc/{params.preqcenv}/g')/{wildcards.condition}; do case \" $SCAN \" in *\" $d \"*) continue;; esac; if [ -d \"$d\" ]; then SCAN=\"$SCAN $d\"; fi; done; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"

else:
    rule multiqc:
        input:  _bamqc('rustqc_mapped'),
                _bamqc('rustqc_uniquemapped'),
                _mapbams,
                versions = "LOGS/versions.txt"
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC")
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
        conda:  "rustqc.yaml"
        container: "oras://jfallmann/monsda:rustqc"
        threads: 1
        params: qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', ""),
                preqcenv = config.get('PREQCENV', 'fastqc')
        shell:  "OUT=$(dirname {output.html}); mkdir -p $OUT; SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; for d in QC/{wildcards.combo}/{wildcards.condition} QC/$(echo {wildcards.combo}|sed 's/rustqc/{params.preqcenv}/g')/{wildcards.condition}; do case \" $SCAN \" in *\" $d \"*) continue;; esac; if [ -d \"$d\" ]; then SCAN=\"$SCAN $d\"; fi; done; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"
