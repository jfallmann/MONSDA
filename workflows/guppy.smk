from MONSDA.Params import basecall_signal_map

CALLERBIN, CALLERENV = env_bin_from_config(config,'BASECALL')

SIGNALS = basecall_signal_map(config, 'guppy')

def basecall_inputs(wildcards):
    return SIGNALS[wildcards.rawfile]

wildcard_constraints:
    rawfile = '|'.join(SAMPLES)

rule themall:
    input: expand("FASTQ/{rawfile}.fastq.gz", rawfile = SAMPLES)

rule call_base:
    input:  f5 = lambda wildcards: basecall_inputs(wildcards)
    output: fq = "FASTQ/{rawfile}.fastq.gz",
            summary = "FASTQ/{rawfile}_summary.txt",
            telemetry = "FASTQ/{rawfile}_telemetry.js"
    log:    "LOGS/{rawfile}/BASECALL/guppy/guppy.log"
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
        mkdir -p "$taskdir/inputs" "$taskdir/out"
        i=0
        for f in {input.f5:q}; do
            i=$((i+1))
            n=$(printf '%04d_%s' "$i" "$(basename "$f")")
            ln -s "$(readlink -f "$f")" "$taskdir/inputs/$n"
            printf '%s\\n' "$n"
        done > "$taskdir/inputs.list"
        {params.caller} {params.cpara} -c {params.cmodel:q} --compress_fastq -i "$taskdir/inputs" --input_file_list "$taskdir/inputs.list" -s "$taskdir/out" 2> {log:q}
        cat "$taskdir"/out/pass/fastq_runid_*.fastq.gz > {output.fq:q}
        cat "$taskdir"/out/*.log >> {log:q}
        mv -f "$taskdir/out/sequencing_summary.txt" {output.summary:q}
        mv -f "$taskdir/out/sequencing_telemetry.js" {output.telemetry:q}
        """