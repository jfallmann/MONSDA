#!/usr/bin/env python
# Code from:https://github.com/vivekbhr/Subread_to_DEXSeq

import sys
import collections
import itertools
import os.path
import optparse
import os
import gzip
import warnings
import traceback as tb
import logging
import inspect

cmd_subfolder = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ),
    "../../lib",
)
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

try:
    scriptname = os.path.basename(inspect.stack()[-1].filename).replace(".py", "")
    log = logging.getLogger(scriptname)

    lvl = log.level if log.level else "INFO"
    for handler in log.handlers[:]:
        handler.close()
        log.removeHandler(handler)

    handler = logging.FileHandler("LOGS/MONSDA.log", mode="a")
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(levelname)-8s %(name)-12s %(message)s",
            datefmt="%m-%d %H:%M",
        )
    )
    log.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(levelname)-8s %(name)-12s %(message)s",
            datefmt="%m-%d %H:%M",
        )
    )
    log.addHandler(handler)
    log.setLevel(lvl)

except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type,
        exc_value,
        exc_tb,
    )
    print("".join(tbe.format()), file=sys.stderr)


optParser = optparse.OptionParser(
    usage="python %prog [options] <in.gtf> <out.gff>",
    description="Script to prepare annotation for DEXSeq."
    + "This script takes an annotation file in Ensembl GTF format"
    + "and outputs a 'flattened' annotation file suitable for use "
    + "with the count_in_exons.py script ",
    epilog="Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology "
    + "Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General "
    + "Public License v3. Part of the 'DEXSeq' package. "
    + "Modified by Vivek Bhardwaj (just a bit!) to write featurecounts gtf as an option.",
)

optParser.add_option(
    "-r",
    "--aggregate",
    type="choice",
    dest="aggregate",
    choices=("no", "yes"),
    default="yes",
    help="'yes' or 'no'. Indicates whether two or more genes sharing an exon should be merged into an 'aggregate gene'. If 'no', the exons that can not be assiged to a single gene are ignored.",
)

# add option for featurecounts output
optParser.add_option(
    "-f",
    "--featurecountsgtf",
    type="string",
    dest="fcgtf",
    action="store",
    help="gtf file to write for featurecounts.",
)
optParser.add_option(
    "-s",
    "--stranded",
    dest="fcstrand",
    action="store_true",
    default=False,
    help="Use strand in gffparse if stranded sequencing protocol was used.",
)
##

(opts, args) = optParser.parse_args()
print(opts, args)

if len(args) != 2:
    sys.stderr.write("Script to prepare annotation for DEXSeq.\n\n")
    sys.stderr.write(
        "Usage: python %s <in.gtf> <out.gff>\n\n" % os.path.basename(sys.argv[0])
    )
    sys.stderr.write("This script takes an annotation file in Ensembl GTF format\n")
    sys.stderr.write("and outputs a 'flattened' annotation file suitable for use\n")
    sys.stderr.write("with the count_in_exons.py script.\n")
    sys.exit(1)

try:
    import HTSeq
except ImportError:
    sys.stderr.write(
        "Could not import HTSeq. Please install the HTSeq Python framework\n"
    )
    sys.stderr.write("available from http://www-huber.embl.de/users/anders/HTSeq\n")
    sys.exit(1)


gtf_file = args[0]
out_file = args[1]

aggregateGenes = opts.aggregate == "yes"

# Step 1: Store all exons with their gene and transcript ID
# in a GenomicArrayOfSets
try:
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=opts.fcstrand)
except:
    warnings.warn(
        "There was an issue setting stranded to True when creating the HTSeq.GenomicArrayOfSets. This is probably due to an error in the annotation allocating one gene ID to both strands. Will retry without strand information, please carefully check your results."
    )
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=False)

for f in HTSeq.GFF_Reader(gtf_file):
    if f.type != "exon":
        continue
    f.attr["gene_id"] = f.attr["gene_id"].replace(":", "_")
    exons[f.iv] += (f.attr["gene_id"], f.attr["transcript_id"])


# Step 2: Form sets of overlapping genes

# We produce the dict 'gene_sets', whose values are sets of gene IDs. Each set
# contains IDs of genes that overlap, i.e., share bases (on the same strand).
# The keys of 'gene_sets' are the IDs of all genes, and each key refers to
# the set that contains the gene.
# Each gene set forms an 'aggregate gene'.

if aggregateGenes == True:
    gene_sets = collections.defaultdict(lambda: set())
    for iv, s in exons.steps():
        # For each step, make a set, 'full_set' of all the gene IDs occuring
        # in the present step, and also add all those gene IDs, whch have been
        # seen earlier to co-occur with each of the currently present gene IDs.
        full_set = set()
        for gene_id, transcript_id in s:
            full_set.add(gene_id)
            full_set |= gene_sets[gene_id]
        # Make sure that all genes that are now in full_set get associated
        # with full_set, i.e., get to know about their new partners
        for gene_id in full_set:
            assert gene_sets[gene_id] <= full_set
            gene_sets[gene_id] = full_set


# Step 3: Go through the steps again to get the exonic sections. Each step
# becomes an 'exonic part'. The exonic part is associated with an
# aggregate gene, i.e., a gene set as determined in the previous step,
# and a transcript set, containing all transcripts that occur in the step.
# The results are stored in the dict 'aggregates', which contains, for each
# aggregate ID, a list of all its exonic_part features.

aggregates = collections.defaultdict(lambda: list())
for iv, s in exons.steps():
    # Skip empty steps
    if len(s) == 0:
        continue
    gene_id = list(s)[0][0]
    ## if aggregateGenes=FALSE, ignore the exons associated to more than one gene ID
    if aggregateGenes == False:
        check_set = set()
        for geneID, transcript_id in s:
            check_set.add(geneID)
        if len(check_set) > 1:
            continue
        else:
            aggregate_id = gene_id
    # Take one of the gene IDs, find the others via gene sets, and
    # form the aggregate ID from all of them
    else:
        assert set(gene_id for gene_id, transcript_id in s) <= gene_sets[gene_id]
        aggregate_id = "+".join(gene_sets[gene_id])
    # Make the feature and store it in 'aggregates'
    f = HTSeq.GenomicFeature(aggregate_id, "exonic_part", iv)
    f.source = os.path.basename(sys.argv[0])
    #   f.source = "camara"
    f.attr = {}
    f.attr["gene_id"] = aggregate_id
    transcript_set = set((transcript_id for gene_id, transcript_id in s))
    f.attr["transcripts"] = "+".join(transcript_set)
    aggregates[aggregate_id].append(f)


# Step 4: For each aggregate, number the exonic parts

aggregate_features = []
for l in aggregates.values():
    for i in range(len(l) - 1):
        if l[i].iv.strand != l[i + 1].iv.strand:
            raise ValueError(
                "Same name found on two strands: %s, %s" % (str(l[i]), str(l[i + 1]))
            )
        assert l[i].name == l[i + 1].name, str(l[i + 1]) + " has wrong name"
        assert l[i].iv.end <= l[i + 1].iv.start, str(l[i + 1]) + " starts too early"
        if l[i].iv.chrom != l[i + 1].iv.chrom:
            raise ValueError(
                "Same name found on two chromosomes: %s, %s"
                % (str(l[i]), str(l[i + 1]))
            )
    aggr_feat = HTSeq.GenomicFeature(
        l[0].name,
        "aggregate_gene",
        HTSeq.GenomicInterval(
            l[0].iv.chrom, l[0].iv.start, l[-1].iv.end, l[0].iv.strand
        ),
    )
    aggr_feat.source = os.path.basename(sys.argv[0])
    aggr_feat.attr = {"gene_id": aggr_feat.name}
    for i in range(len(l)):
        l[i].attr["exonic_part_number"] = "%03d" % (i + 1)
    aggregate_features.append(aggr_feat)


# Step 5: Sort the aggregates, then write everything out

aggregate_features.sort(key=lambda f: (f.iv.chrom, f.iv.start))

if ".gz" in out_file:
    with gzip.open(out_file, "wb") as fout:
        for aggr_feat in aggregate_features:
            fout.write(bytes(aggr_feat.get_gff_line(), encoding="UTF8"))
            for f in aggregates[aggr_feat.name]:
                fout.write(bytes(f.get_gff_line(), encoding="UTF8"))

else:
    with open(out_file, "w") as fout:
        for aggr_feat in aggregate_features:
            fout.write(aggr_feat.get_gff_line())
            for f in aggregates[aggr_feat.name]:
                fout.write(f.get_gff_line())

## modify file to print gtf if featurecounts gtf requested
fcountgtf = opts.fcgtf

if fcountgtf:
    if ".gz" in fcountgtf:
        os.system(
            "zcat "
            + out_file
            + "|sed s/aggregate_gene/gene/g|sed s/exonic_part/exon/g |cut -d';' -f1|gzip > "
            + fcountgtf
        )
    else:
        os.system(
            "sed s/aggregate_gene/gene/g "
            + out_file
            + "|sed s/exonic_part/exon/g |cut -d';' -f1 > "
            + fcountgtf
        )
