if rundedup:
    if paired == 'paired':
        if prededup:
            rule multiqc:
                input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_dedupmapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_uniquededup.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique"])
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/{condition}_multiqc.log"
                conda:  "rustqc.yaml"
                container: "oras://jfallmann/monsda:rustqc"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"
        else:
            rule multiqc:
                input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_dedupmapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_uniquededup.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique"])
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/{condition}_multiqc.log"
                conda:  "rustqc.yaml"
                container: "oras://jfallmann/monsda:rustqc"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"

    else:
        if prededup:
            rule multiqc:
                input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_dedupmapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_uniquededup.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique"])
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/{condition}_multiqc.log"
                conda:  "rustqc.yaml"
                container: "oras://jfallmann/monsda:rustqc"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"
        else:
            rule multiqc:
                input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_dedupmapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.rustqc_uniquededup.output.js, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique"])
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/{condition}_multiqc.log"
                conda:  "rustqc.yaml"
                container: "oras://jfallmann/monsda:rustqc"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"

else:
    if paired == 'paired':
        rule multiqc:
            input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo)
            output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                    lst = "QC/Multi/{combo}/{condition}/qclist.txt"
            log:    "LOGS/{combo}/{condition}_multiqc.log"
            conda:  "rustqc.yaml"
            container: "oras://jfallmann/monsda:rustqc"
            threads: 1
            params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
            shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"

    else:
        rule multiqc:
            input:  expand(rules.rustqc_mapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.rustqc_uniquemapped.output.js, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo)
            output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                    lst = "QC/Multi/{combo}/{condition}/qclist.txt"
            log:    "LOGS/{combo}/{condition}_multiqc.log"
            conda:  "rustqc.yaml"
            container: "oras://jfallmann/monsda:rustqc"
            threads: 1
            params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
            shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"
