# bbduk.smk ---
#
# Filename: bbduk.smk
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep  4 15:26:24 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Thu Oct 11 10:09:29 2018 (+0200)
#           By: Joerg Fallmann
#     Update #: 3
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
rule bbduk_trim:
        input:  "FASTQ/{file}.fastq.gz", "QC/{file}_fastqc.zip" if config["QC"] == "ON" else  "FASTQ/{file}.fastq.gz"
        output: "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
        log:    "LOGS/{file}_trimmed.log"
        params: ada={ADAPTERS}
        conda: "../envs/trimm.yaml"
        threads: 1
        shell:  "bbduk.sh in={input[0]} out={output} ref={params.ada} ktrim=r k=23 mink=11 hdist=1 tpe tbo"
#
# bbduk.smk ends here
