# trimming.smk ---
#
# Filename: trimming.smk
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep  4 14:50:40 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Feb 27 15:35:14 2019 (+0100)
#           By: Joerg Fallmann
#     Update #: 26
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
#q
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
rule cutadapt_trim:
    input:  "FASTQ/{file}.fastq.gz", "QC/{file}_fastqc.zip" if "ON" in config["QC"] else  "FASTQ/{file}.fastq.gz"
    output: "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    log:    "LOGS/{file}_trimmed.log"
    params: ada=ADAPTERS
    conda: "../envs/trimm.yaml"
    threads: 1
    shell:  "cutadapt -a file:{params.ada} -q25 -M 95 -m 8 -e 0.15 -o {output} {input[0]} > {log}"

#
# trimming.smk ends here
