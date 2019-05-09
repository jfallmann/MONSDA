# trimmomatic.smk ---
#
# Filename: trimmomatic.smk
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep  4 15:26:04 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Feb 27 11:08:27 2019 (+0100)
#           By: Joerg Fallmann
#     Update #: 4
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
rule trimmomatic_trim:
        input:  "FASTQ/{file}.fastq.gz", "QC/{file}_fastqc.zip" if config["QC"] == "ON" else  "FASTQ/{file}.fastq.gz"
        output: "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
        log:    "LOGS/{file}_trimmed.log"
        threads: 1
        conda: "../envs/trimm.yaml"
        shell:  "java -jar {BINS}/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 {input[0]} {output} -threads {threads} -trimlog {log} MINLEN:25"
#
# trimmomatic.smk ends here
