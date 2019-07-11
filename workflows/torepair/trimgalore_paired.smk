# trimgalore_paired.smk ---
#
# Filename: trimgalore_paired.smk
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Apr 30 10:19:14 2019 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Apr 30 10:19:38 2019 (+0200)
#           By: Joerg Fallmann
#     Update #: 1
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
#
#
#
#

# Change Log:
#
#
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
#
#

# Code:
rule trimgalore_trim_paired:
    input:  "FASTQ/{file}.fastq.gz", "QC/{file}_fastqc.zip" if "ON" in config["QC"] else "FASTQ/{file}.fastq.gz"
    output: "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/trimgalore.yaml"
    threads: 1
    params: odir=lambda wildcards,output:os.path.dirname(output[0])
    shell:  "trim_galore --no_report_file --gzip -q 25 --length 8 -e 0.15 -o {params.odir} {input[0]} > {log[0]}"

rule trimgalore_rename_paired:
    input:  "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
    output: "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    shell:  "mv {input} {output}"

#
# trimgalore_paired.smk ends here
