# Fixtures for tests/tool_audit/audit_tool_flags.py

## tiny.bam

Minimal coordinate-sorted BAM (3 reads, 2 alignment positions on a fake
`chr1`, one duplicated position) used for the `umicollapse` `minimal_run`
check, since `umicollapse --help` exits non-zero with "Missing input
file!" before printing usage and cannot be flag-checked via help text.

Regenerate with (inside the `apptainer` conda env):

    cat > tiny.sam << 'EOF'
    @HD	VN:1.6	SO:coordinate
    @SQ	SN:chr1	LN:1000
    r1_ACGTACGT	0	chr1	1	60	4M	*	0	0	ACGT	IIII
    r2_ACGTACGA	0	chr1	1	60	4M	*	0	0	ACGT	IIII
    r3_TTTTACGT	0	chr1	5	60	4M	*	0	0	TTTT	IIII
    EOF
    apptainer exec ~/Work2/Container/MONSDA/umicollapse-1.5.0.sif \
        samtools view -b -o tiny.bam tiny.sam

Verified this dedups to 2 reads via:

    apptainer exec --bind $PWD:$PWD ~/Work2/Container/MONSDA/umicollapse-1.5.0.sif \
        umicollapse bam -i $PWD/tiny.bam -o $PWD/out.bam --paired
