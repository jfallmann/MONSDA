if rundedup:
    if paired == 'paired':
        if prededup:
            rule multiqc:
                input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo),
                        expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),
                        expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),
                        expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique", "sorted_dedup", "sorted_unique_dedup"])
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/{condition}_multiqc.log"
                conda:  "qc.yaml"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"
        else:
            rule multiqc:
                input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo),
                        expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),                    
                        expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique", "sorted_dedup", "sorted_unique_dedup"])
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/{condition}_multiqc.log"
                conda:  "qc.yaml"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"

    else:
        if prededup:
            rule multiqc:
                input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo),
                        expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo,type=["sorted", "sorted_unique", "sorted_dedup", "sorted_unique_dedup"])  
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/{condition}_multiqc.log"
                conda:  "qc.yaml"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"                    
        else:
            rule multiqc:
                input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo),
                        expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), combo=combo),                 
                        expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique", "sorted_dedup", "sorted_unique_dedup"])
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/{condition}_multiqc.log"
                conda:  "qc.yaml"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"

else:
    if paired == 'paired':
        rule multiqc:
            input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo),
                    expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),
                    expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo)
            output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                    lst = "QC/Multi/{combo}/{condition}/qclist.txt"
            log:    "LOGS/{combo}/{condition}_multiqc.log"
            conda:  "qc.yaml"
            threads: 1
            params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
            shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"

    else:
        rule multiqc:
            input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo),
                    expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo)
            output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                    lst = "QC/Multi/{combo}/{condition}/qclist.txt"
            log:    "LOGS/{combo}/{condition}_multiqc.log"
            conda:  "qc.yaml"
            threads: 1
            params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
            shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"
