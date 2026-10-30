DEDUPBIN, DEDUPENV = env_bin_from_config(config, 'DEDUP')

wlparams = tool_params(SAMPLES[0], None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('WHITELIST')
eparams = tool_params(SAMPLES[0], None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('EXTRACT')

wildcard_constraints:
    type = "sorted|sorted_unique"

if DEDUP_PAIRED == 'paired':
    if wlparams:
        rule whitelist:
            input:  r1 = lambda wildcards: raw_fastq_input(wildcards, 'R1'),
                    r2 = lambda wildcards: raw_fastq_input(wildcards, 'R2')
            output: wl = "DEDUP_FASTQ/{combo}/{file}_whitelist",
                    td = temp(directory("TMP/UMIWL/{combo}/{file}"))
            log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_whitelist.log"
            conda: ""+DEDUPENV+".yaml"
            container: "oras://jfallmann/monsda:"+DEDUPENV+""
            threads: 1
            params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('WHITELIST', ""),
                    dedup = DEDUPBIN
            shell:  "mkdir -p {output.td} && {params.dedup} whitelist {params.dpara} --temp-dir {output.td} --log={log} --stdin={input.r1} --read2-in={input.r2} --stdout={output.wl}"

        if eparams:
            rule extract:
                input:  r1 = lambda wildcards: raw_fastq_input(wildcards, 'R1'),
                        r2 = lambda wildcards: raw_fastq_input(wildcards, 'R2'),
                        wl = rules.whitelist.output.wl
                output: o1 = "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                        o2 = "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz",
                        td = temp(directory("TMP/UMIEX/{combo}/{file}"))
                log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_extract.log"
                conda: ""+DEDUPENV+".yaml"
                container: "oras://jfallmann/monsda:"+DEDUPENV+""
                threads: 1
                params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('EXTRACT', ""),
                        dedup = DEDUPBIN
                shell:  "mkdir -p {output.td} && {params.dedup} extract {params.dpara} --temp-dir {output.td} --log={log} --error-correct-cell --whitelist={input.wl} --stdin={input.r1} --read2-in={input.r2} --stdout={output.o1} --read2-out={output.o2}"
        else:
            rule extract:
                input:  r1 = lambda wildcards: raw_fastq_input(wildcards, 'R1'),
                        r2 = lambda wildcards: raw_fastq_input(wildcards, 'R2'),
                        wl = rules.whitelist.output.wl
                output: o1 = "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                        o2 = "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz",
                        td = temp(directory("TMP/UMIEX/{combo}/{file}"))
                log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_extract.log"
                threads: 1
                shell:  "mkdir -p {output.td} && echo 'EXTRACT unset, linking input fastq(s) as-is' > {log} && ln -sf $(readlink -f {input.r1}) {output.o1} 2>> {log} && ln -sf $(readlink -f {input.r2}) {output.o2} 2>> {log}"
    else:
        if eparams:
            rule extract:
                input:  r1 = lambda wildcards: raw_fastq_input(wildcards, 'R1'),
                        r2 = lambda wildcards: raw_fastq_input(wildcards, 'R2')
                output: o1 = "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                        o2 = "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz",
                        td = temp(directory("TMP/UMIEX/{combo}/{file}"))
                log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_extract.log"
                conda: ""+DEDUPENV+".yaml"
                container: "oras://jfallmann/monsda:"+DEDUPENV+""
                threads: 1
                params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('EXTRACT', ""),
                        dedup = DEDUPBIN
                shell:  "mkdir -p {output.td} && {params.dedup} extract {params.dpara} --temp-dir {output.td} --log={log} --stdin={input.r1} --read2-in={input.r2} --stdout={output.o1} --read2-out={output.o2}"
        else:
            rule extract:
                input:  r1 = lambda wildcards: raw_fastq_input(wildcards, 'R1'),
                        r2 = lambda wildcards: raw_fastq_input(wildcards, 'R2')
                output: o1 = "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                        o2 = "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz",
                        td = temp(directory("TMP/UMIEX/{combo}/{file}"))
                log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_extract.log"
                threads: 1
                shell:  "mkdir -p {output.td} && echo 'EXTRACT unset, linking input fastq(s) as-is' > {log} && ln -sf $(readlink -f {input.r1}) {output.o1} 2>> {log} && ln -sf $(readlink -f {input.r2}) {output.o2} 2>> {log}"

else:
    if wlparams:
        rule whitelist:
            input:  r1 = lambda wildcards: raw_fastq_input(wildcards)
            output: wl = "DEDUP_FASTQ/{combo}/{file}_whitelist",
                    td = temp(directory("TMP/UMIWL/{combo}/{file}"))
            log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_whitelist.log"
            conda: ""+DEDUPENV+".yaml"
            container: "oras://jfallmann/monsda:"+DEDUPENV+""
            threads: 1
            params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('WHITELIST', ""),
                    dedup = DEDUPBIN
            shell:  "mkdir -p {output.td} && {params.dedup} whitelist {params.dpara} --temp-dir {output.td} --log={log} --stdin={input.r1} --stdout={output.wl}"

        if eparams:
            rule extract:
                input:  r1 = lambda wildcards: raw_fastq_input(wildcards),
                        wl = rules.whitelist.output.wl
                output: o1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz",
                        td = temp(directory("TMP/UMIEX/{combo}/{file}"))
                log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_extract.log"
                conda: ""+DEDUPENV+".yaml"
                container: "oras://jfallmann/monsda:"+DEDUPENV+""
                threads: 1
                params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('EXTRACT', ""),
                        dedup = DEDUPBIN
                shell:  "mkdir -p {output.td} && {params.dedup} extract {params.dpara} --temp-dir {output.td} --log={log} --error-correct-cell --whitelist={input.wl} --stdin={input.r1} --stdout={output.o1}"
        else:
            rule extract:
                input:  r1 = lambda wildcards: raw_fastq_input(wildcards),
                        wl = rules.whitelist.output.wl
                output: o1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz",
                        td = temp(directory("TMP/UMIEX/{combo}/{file}"))
                log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_extract.log"
                threads: 1
                shell:  "mkdir -p {output.td} && echo 'EXTRACT unset, linking input fastq(s) as-is' > {log} && ln -sf $(readlink -f {input.r1}) {output.o1} 2>> {log}"

    else:
        if eparams:
            rule extract:
                input:  r1 = lambda wildcards: raw_fastq_input(wildcards)
                output: o1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz",
                        td = temp(directory("TMP/UMIEX/{combo}/{file}"))
                log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_extract.log"
                conda: ""+DEDUPENV+".yaml"
                container: "oras://jfallmann/monsda:"+DEDUPENV+""
                threads: 1
                params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('EXTRACT', ""),
                        dedup = DEDUPBIN
                shell:  "mkdir -p {output.td} && {params.dedup} extract {params.dpara} --temp-dir {output.td} --log={log} --stdin={input.r1} --stdout={output.o1}"
        else:
            rule extract:
                input:  r1 = lambda wildcards: raw_fastq_input(wildcards)
                output: o1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz",
                        td = temp(directory("TMP/UMIEX/{combo}/{file}"))
                log:   "LOGS/{combo}/{file}/DEDUP/umitools/dedup_extract.log"
                threads: 1
                shell:  "mkdir -p {output.td} && echo 'EXTRACT unset, linking input fastq(s) as-is' > {log} && ln -sf $(readlink -f {input.r1}) {output.o1} 2>> {log}"

if DEDUP_PAIRED == 'paired':
    rule dedupbam:
        input:  bam = "MAPPED/{combo}/{file}_mapped_{type}.bam"
        output: bam = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", category="DEDUP"),
                bai = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam.bai", category="DEDUP"),
                td = temp(directory("TMP/UMIDD/{combo}/{file}_{type}"))
        log:    "LOGS/{combo}/{file}/DEDUP/umitools/dedupbam_{type}.log"
        conda:  ""+DEDUPENV+".yaml"
        container: "oras://jfallmann/monsda:"+DEDUPENV+""
        threads: 1
        priority: 0               # This should be done after all mapping is done
        params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
                dedup = DEDUPBIN
        shell: "mkdir -p {output.td} && {params.dedup} dedup {params.dpara} --paired --temp-dir {output.td} --stdin={input.bam} --log={log} --stdout={output.bam} 2>> {log} && samtools index {output.bam} 2>> {log}"
else:
    rule dedupbam:  
        input:  bam = "MAPPED/{combo}/{file}_mapped_{type}.bam"
        output: bam = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", category="DEDUP"),
                bai = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam.bai", category="DEDUP"),
                td = temp(directory("TMP/UMIDD/{combo}/{file}_{type}"))
        log:    "LOGS/{combo}/{file}/DEDUP/umitools/dedupbam_{type}.log"
        conda:  ""+DEDUPENV+".yaml"
        container: "oras://jfallmann/monsda:"+DEDUPENV+""
        threads: 1
        priority: 0               # This should be done after all mapping is done
        params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
                dedup = DEDUPBIN
        shell: "mkdir -p {output.td} && {params.dedup} dedup {params.dpara} --temp-dir {output.td} --stdin={input.bam} --log={log} --stdout={output.bam} 2>> {log} && samtools index {output.bam} 2>> {log}"
