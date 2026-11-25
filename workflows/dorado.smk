from MONSDA.Params import basecall_signal_map

CALLERBIN, CALLERENV = env_bin_from_config(config,'BASECALL')

SIGNALS = basecall_signal_map(config, 'dorado')

def basecall_inputs(wildcards):
    return SIGNALS[wildcards.rawfile]

wildcard_constraints:
    rawfile = '|'.join(SAMPLES)

rule themall:
    input: expand("FASTQ/{rawfile}.fastq.gz", rawfile = SAMPLES)

rule call_base:
    input:  p5 = lambda wildcards: basecall_inputs(wildcards)
    output: fq = "FASTQ/{rawfile}.fastq.gz",
            bam = temp("FASTQ/{rawfile}.bam")
    log:    "LOGS/{rawfile}/BASECALL/dorado/dorado.log"
    conda:  ""+CALLERENV+".yaml"
    container: "oras://jfallmann/monsda:"+CALLERENV+""
    threads: MAXTHREAD
    params: caller = CALLERBIN,
            cpara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'BASECALL', CALLERENV)['OPTIONS'].get('BASECALL', ""),
            cmodel = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'BASECALL', CALLERENV)['OPTIONS'].get('MODEL', ""),
            fqdir = lambda wildcards, output: os.path.dirname(output.fq)
    shell:
        """
        set -euo pipefail
        taskdir=$(mktemp -d)
        trap 'rm -rf "$taskdir"' EXIT
        mkdir -p "$taskdir/inputs" "$taskdir/models"
        i=0
        for f in {input.p5:q}; do
            i=$((i+1))
            ln -s "$(readlink -f "$f")" "$taskdir/inputs/$(printf '%04d_%s' "$i" "$(basename "$f")")"
        done
        {params.caller} download --directory "$taskdir/models" --model {params.cmodel:q} &> {log:q}
        {params.caller} basecaller {params.cpara} "$taskdir/models/{params.cmodel}" "$taskdir/inputs" 2>> {log:q} 1> {output.bam:q}
        samtools view -h {output.bam:q} | samtools fastq -n - | pigz > {output.fq:q}
        """