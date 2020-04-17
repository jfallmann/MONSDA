#!/usr/bin/env python3

"""DIEGO
"""

# Gero Doose gero@bioinf.uni-leipzig.de
import sys
import argparse
import operator
import math
from datetime import datetime
import itertools
import numpy as np
from scipy.stats import ranksums
from scipy.stats import nbinom
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='DIEGO (DIfferential altErnative splicinG detectiOn) \n')
parser.add_argument("-a", "--table", dest='my_a', required=True, type=argparse.FileType('r'),
                    help="table of splice junction supports per sample with gene information "
                         "(created with pre_std.py, pre_star.py or pre_segemehl.py)")
parser.add_argument("-b", "--list", dest='my_b', required=True, type=argparse.FileType('r'),
                    help="condition to sample relation in the format: condition tab-delimiter sampleName "
                         "(you can specify a subset of the samples included in the table, "
                         "however the sampleName has to be the same as used in table) "
                         "(you can specify a label in the third column that is used (if present) "
                         "when the clustering mode is executed)")
parser.add_argument("-x", "--base_condition", dest='my_base_condition', required=True,
                    help="specify base condition (direction of change between the two conditions)")
parser.add_argument("-c", "--minsupp", dest='my_ms', type=int,
                    help="min support per splice site (at least -d samples have to show this min support)(default:10)")
parser.add_argument("-d", "--minsamples", dest='my_mt', type=int,
                    help="min amount of samples showing  in at least one of the junctions the min supp (default:3)")
parser.add_argument("-q", "--significanceThreshold", dest='my_qvalue', type=float,
                    help="significance level alpha (default 0.01)")
parser.add_argument("-z", "--foldchangeThreshold", dest='my_fc', type=float,
                    help="abundance change threshold (default 1.0)")
parser.add_argument("-e", "--cluster", action="store_true", dest='cluster',
                    help="provide -e if clustering should be performed")
parser.add_argument("-f", "--dendogram", dest='my_f', help="prefix specifying the dendogram plot")
parser.add_argument("-r", "--random", action="store_true", dest='my_r',
                    help="provide -r if random seed should be used (not deterministic mode)")


args = parser.parse_args()
if args.my_ms:
    minsupp = args.my_ms
else:
    minsupp = 10
if args.my_mt:
    min_samples = args.my_mt
else:
    min_samples = 3
if args.my_qvalue:
    myqvalue = args.my_qvalue
else:
    myqvalue = 0.01
if args.my_fc:
    myfc = args.my_fc
else:
    myfc = 1.0
perform_cluster = args.cluster
if args.my_r:
    my_seed = datetime.now().microsecond
else:
    my_seed = 102261126

program_version = "0.1.2"

sys.stderr.write("==============================================\n")
sys.stderr.write("====================DIEGO=====================\n\n")
sys.stderr.write("  _______     _____   _________     ______      ____                 \n")
sys.stderr.write(" |_      `.  |_   _| |_   ___  |  .' ___  |   .'    `.                \n")
sys.stderr.write("   |   __  \   | |     | |_  \_| /  /   \_|  |   __   |                \n")
sys.stderr.write("   |  /  \ |   | |     |  _|  _  | |   ____  |  /  \  |                 \n")
sys.stderr.write("  _|  \__/ /  _| |_   _| |___/ | \  \___]  | \  \__/  /                \n")
sys.stderr.write(" |_______.'  |_____| |_________|  `._____.'   `.____.'               \n\n")
sys.stderr.write("====================DIEGO=====================\n")
sys.stderr.write("==============================================\n")
if perform_cluster:
    if args.my_f:
        dendro_name = args.my_f
    else:
        dendro_name = "cluster_dendrogram"
    sys.stderr.write("You are running the version {} of DIEGO "
                     "(DIfferential altErnative splicinG detectiOn) in cluster mode\n".format(program_version))

else:
    sys.stderr.write("You are running the version {} of DIEGO "
                     "(DIfferential altErnative splicinG detectiOn) in normal mode\n".format(program_version))
sys.stderr.write("Parameters: minimum read support: {}, min amount of samples per condition: {}, "
                 "minimum abundance change cutoff: {}, significance threshold: {}, base condition: {}, "
                 "random seed: {} \n".format(minsupp, min_samples, myfc, myqvalue, args.my_base_condition, my_seed))

geneH = {}
scH = {}
design = []
baseC = ""
compC = ""
useClusterNames = False
sampleList = []
sampleListCluster = []
p_valueL = []
correspondingL = []


class Junction(object):
    """Subclass of Module

    Junction
    """
    left = -1
    right = -1
    significant = "none"
    p_value = -1
    q_value = -1
    ab_change = 0
    jtype = "none"
    zero_replaced = "none"
    bio_sig = -1

    def __init__(self, start, end, jtype):
        self.left = int(start)
        self.right = int(end)
        self.significant = "no"
        self.p_value = 1
        self.q_value = 1
        self.ab_change = 0
        self.jtype = jtype
        self.zero_replaced = "False"
        self.bio_sig = -1


def make_Junction(left, right, jtype):
    """initialize Junction object

    """
    myJunction = Junction(left, right, jtype)
    return myJunction


class Gene(object):
    """Subclass of Module

    Gene
    """
    ensg = ""
    name = ""
    chrom = ""
    conditions = {}
    junctionL = []
    centreH = {}
    centreH_clr = {}
    centreDist = -1

    def __init__(self, ensg, name, chrom):
        self.chrom = chrom
        self.ensg = ensg
        self.name = name
        self.conditions = {}
        self.junctionL = []
        self.centreH = {}
        self.centreH_clr = {}
        self.centreDist = -1


def make_Gene(ensg, name, chrom):
    """initialize Gene object

    """
    myGene = Gene(ensg, name, chrom)
    return myGene


class Gene_cluster(object):
    """Subclass of Module

    Genecluster
    """
    ensg = ""
    name = ""
    chrom = ""
    sampleH = {}
    sampleL = {}
    junctionL = []
    sampleDistL = []

    # the class constructor
    def __init__(self, ensg, name, chrom):
        self.chrom = chrom
        self.ensg = ensg
        self.sampleH = {}
        self.junctionL = []
        self.sampleDistL = []
        self.name = name
        self.sampleL = []


def make_Gene_cluster(ensg, name, chrom):
    """initialize Gene cluster object

    """
    myGene = Gene_cluster(ensg, name, chrom)
    return myGene


def read_sample_condition_info():
    """creates the design.

    """
    global baseC, compC
    condiH = {}
    for line in args.my_b:
        cols = line.strip().split('\t')
        if cols[1] in scH:
            sys.stderr.write("sample is multiple times present: "+str(cols[1])+" \n")
            sys.exit(1)
        else:
            scH[cols[1]] = cols[0]
            condiH[cols[0]] = 1
    if len(condiH) > 2:
        sys.stderr.write("more than two conditions specified with -b\n")
        sys.exit(1)
    else:
        for key in condiH:
            if key == args.my_base_condition:
                baseC = str(key)
                firstB = True
            else:
                compC = str(key)
        if not firstB:
            sys.stderr.write("specified base condition is not present in file provided by -b\n")
            sys.exit(1)


def read_sample_info_cluster():
    """reads cluster info name.

    """
    global useClusterNames
    for line in args.my_b:
        cols = line.strip().split('\t')
        sampleList.append(cols[1])
        scH[cols[1]] = 1
        if len(cols) > 2:
            useClusterNames = True
            sampleListCluster.append(cols[2])
        else:
            useClusterNames = False


def create_design(line):
    """creates the design.

    """
    cols = line.strip().split('\t')
    for idx, ele in enumerate(cols):
        if 1 < idx < len(cols) - 2:
            if ele in scH:
                design.append(ele)
            else:
                design.append("skip_me")


def read_table():
    """creates the design.

    """
    read_sample_condition_info()
    firstLine = True
    for line in args.my_a:
        if firstLine:
            create_design(line)
        else:
            cols = line.strip().split('\t')
            chrom = cols[0].strip().split(':')
            leftright = chrom[1].strip().split('-')
            jtype = cols[1]
            geneID = cols[-2]
            geneName = cols[-1]
            if geneID not in geneH:
                geneH[geneID] = make_Gene(geneID, geneName, chrom[0])
            geneH[geneID].junctionL.append(make_Junction(leftright[0], leftright[1], jtype))
            for idxx, ele in enumerate(cols):
                if 1 < idxx < len(cols) - 2:
                    idx = idxx - 2
                    if design[idx] != "skip_me":
                        if scH[design[idx]] not in geneH[geneID].conditions:
                            geneH[geneID].conditions[scH[design[idx]]] = {}
                        if design[idx] not in geneH[geneID].conditions[scH[design[idx]]]:
                            geneH[geneID].conditions[scH[design[idx]]][design[idx]] = []
                        geneH[geneID].conditions[scH[design[idx]]][design[idx]].append(int(ele))
        firstLine = False


def read_table_cluster():
    """creates the design.

    """
    read_sample_info_cluster()
    firstLine = True
    for line in args.my_a:
        if firstLine:
            create_design(line)
        else:
            cols = line.strip().split('\t')
            chrom = cols[0].strip().split(':')
            leftright = chrom[1].strip().split('-')
            jtype = cols[1]
            geneID = cols[-2]
            geneName = cols[-1]
            if geneID not in geneH:
                geneH[geneID] = make_Gene_cluster(geneID, geneName, chrom[0])
            junction_obj = make_Junction(leftright[0], leftright[1], jtype)
            geneH[geneID].junctionL.append(junction_obj)
            for idxx, ele in enumerate(cols):
                if 1 < idxx < len(cols) - 2:
                    idx = idxx - 2
                    if design[idx] != "skip_me":
                        if design[idx] not in geneH[geneID].sampleH:
                            geneH[geneID].sampleH[design[idx]] = []
                        geneH[geneID].sampleH[design[idx]].append(int(ele))
        firstLine = False


def zeroReplacement_multiplicative(gene):
    """replaces the zeros in count matrix.

    :type gene: Gene
    """
    mdelta = 0.5/len(gene.junctionL)
    for condi in gene.conditions:
        for sample in gene.conditions[condi]:
            comp_sum = 0
            zero_counter = 0
            for ele in gene.conditions[condi][sample]:
                comp_sum += ele
                if ele == 0:
                    zero_counter += 1
            for idx, ele in enumerate(gene.conditions[condi][sample]):
                if ele == 0:
                    gene.junctionL[idx].zero_replaced = "True"
                    gene.conditions[condi][sample][idx] = mdelta
                else:
                    gene.conditions[condi][sample][idx] = (1-zero_counter*mdelta)*(float(ele)/comp_sum)


def zeroReplacement_nbinom(gene, nbinnormalL):
    """replaces the zeros in count matrix.

    :type gene: Gene
    :type nbinnormalL: list
    """
    idx_nbinom = 0
    for condi in gene.conditions:
        sample_amount = len(gene.conditions[condi])
        count_replacementL = [0]*len(gene.junctionL)
        for sample in gene.conditions[condi]:
            comp_sum = 0
            for idx, ele in enumerate(gene.conditions[condi][sample]):
                if ele == 0:
                    idx_nbinom += 1
                    count_replacementL[idx] += 1
                    gene.conditions[condi][sample][idx] = nbinnormalL[idx_nbinom]
                comp_sum += gene.conditions[condi][sample][idx]
            for idx, ele in enumerate(gene.conditions[condi][sample]):
                gene.conditions[condi][sample][idx] = (float(ele)/comp_sum)
        for _idx in range(len(gene.junctionL)):
            if count_replacementL[_idx] >= sample_amount/2:
                gene.junctionL[_idx].zero_replaced = "True"


def zeroReplacement_nbinom_cluster(gene, nbinnormalL):
    """replaces the zeros in count matrix.

    :type gene: Gene
    :type nbinnormalL: list
    """
    sample_amount = len(gene.sampleH)
    count_replacementL = [0]*len(gene.junctionL)
    idx_nbinom = 0
    for sample in gene.sampleH:
        comp_sum = 0
        for idx, ele in enumerate(gene.sampleH[sample]):
            if ele == 0:
                idx_nbinom += 1
                count_replacementL[idx] += 1
                gene.sampleH[sample][idx] = nbinnormalL[idx_nbinom]
            comp_sum += gene.sampleH[sample][idx]
        for idx, ele in enumerate(gene.sampleH[sample]):
            gene.sampleH[sample][idx] = (float(ele)/comp_sum)
    for _idx in range(len(gene.junctionL)):
        if count_replacementL[_idx] >= sample_amount/2:
            gene.junctionL[_idx].zero_replaced = "True"


def reduceMatrix(gene, nbinnormalL):
    """reduces the count matrix.

    :type gene: Gene
    :type nbinnormalL: list
    :rtype reference_file: Bool
    """
    minjunctions = 2
    maxjunctions = 60
    keep_colH = {}
    keep_colHpre = {}
    for condi in gene.conditions:
        toRemoveC = []
        for sample in gene.conditions[condi]:
            keep_composition = False
            for idx, ele in enumerate(gene.conditions[condi][sample]):
                if ele >= minsupp:
                    keep_composition = True
                    if idx not in keep_colHpre:
                        keep_colHpre[idx] = 0
                    keep_colHpre[idx] += 1
            if not keep_composition:
                toRemoveC.append(sample)
        for sample in toRemoveC:
            del gene.conditions[condi][sample]
    for key, value in keep_colHpre.items():
        if value >= min_samples:
            keep_colH[key] = 1
    myIndices = sorted(keep_colH.keys())
    for condi in gene.conditions:
        for sample in gene.conditions[condi]:
            gene.conditions[condi][sample] = [gene.conditions[condi][sample][i] for i in myIndices]
    gene.junctionL = [gene.junctionL[i] for i in myIndices]
    if len(gene.junctionL) < minjunctions or len(gene.junctionL) > maxjunctions:
        return True
    else:
        for condi in gene.conditions:
            if len(gene.conditions[condi]) < min_samples:
                return True
    zeroReplacement_nbinom(gene, nbinnormalL)
    return False


def reduceMatrix_cluster(gene, nbinnormalL):
    """reduces the count matrix.

    :type gene: Gene
    :type nbinnormalL: list
    :rtype reference_file: Bool
    """
    minjunctions = 2
    maxjunctions = 60
    keep_colH = {}
    keep_colHpre = {}
    toRemoveC = []
    for sample in gene.sampleH:
        keep_composition = False
        for idx, ele in enumerate(gene.sampleH[sample]):
            if ele >= minsupp:
                keep_composition = True
                if idx not in keep_colHpre:
                    keep_colHpre[idx] = 0
                keep_colHpre[idx] += 1
        if not keep_composition:
            toRemoveC.append(sample)
    for sample in toRemoveC:
        del gene.sampleH[sample]
    for key, value in keep_colHpre.items():
        if value >= min_samples:
            keep_colH[key] = 1
    myIndices = sorted(keep_colH.keys())
    for sample in gene.sampleH:
        gene.sampleH[sample] = [gene.sampleH[sample][i] for i in myIndices]
    gene.junctionL = [gene.junctionL[i] for i in myIndices]
    if len(gene.junctionL) < minjunctions or len(gene.junctionL) > maxjunctions:
        return True
    else:
        if len(gene.sampleH) < min_samples:
            return True
    zeroReplacement_nbinom_cluster(gene, nbinnormalL)
    return False


def dist_eukledian(firstL, secondL):
    """calc eukledian
    :type firstL: list
    :type secondL: list
    :rtype : float
    """
    diste = 0
    for i in range(len(firstL)):
        diste += math.pow(firstL[i] - secondL[i], 2)
    return math.pow(diste, 0.5)


def build_centre_clr_transformed_cluster(gene):
    """build cluster

    :type gene: Gene
    :rtype center_list: list
    """
    center_list = []
    for i in range(len(gene.junctionL)):
        column_values = []
        for sample in gene.sampleH:
            column_values.append(gene.sampleH[sample][i])
        center_list.append(arithmean(column_values))
    return center_list


def build_centre_clr_transformed(gene, condi):
    """Info
    :type gene: Gene
    :type condi: str
    """
    center_list = []
    for i in range(len(gene.junctionL)):
        column_values = []
        for sample in gene.conditions[condi]:
            column_values.append(gene.conditions[condi][sample][i])
        center_list.append(arithmean(column_values))
    return center_list


def geomean(my_list):
    """calc geomean
    :type my_list: list
    :rtype : float
    """
    gs = 1.0
    for ele in my_list:
        gs = gs * ele
    if gs <= 0.0:
        sys.stderr.write(" list product is zero or negative! (in geomean()) \n")
        sys.stderr.write(str(my_list)+"\n")
        sys.exit(1)
    elif len(my_list) == 0:
        sys.stderr.write(" list is empty (in geomean()) \n")
        sys.stderr.write(str(my_list)+"\n")
    else:
        return math.pow(gs, 1.0/len(my_list))


def arithmean(my_list):
    """calc geomean
    :type my_list: list
    :rtype : float
    """
    if len(my_list) > 0:
        return float(sum(my_list))/len(my_list)
    else:
        sys.stderr.write(" ERROR empty List to arithmetic mean")
        sys.exit(1)


def transform(gene):
    """transform
    :type gene: Gene
    """
    for condi in gene.conditions:
        for sample in gene.conditions[condi]:
            geom_mean = geomean(gene.conditions[condi][sample])
            for idx, ele in enumerate(gene.conditions[condi][sample]):
                gene.conditions[condi][sample][idx] = (math.log((ele/geom_mean)))


def transform_cluster(gene):
    """transform
    :type gene: Gene
    """
    for sample in gene.sampleH:
        geom_mean = geomean(gene.sampleH[sample])
        for idx, ele in enumerate(gene.sampleH[sample]):
            gene.sampleH[sample][idx] = (math.log((ele/geom_mean)))


def dist_simplex(firstL, secondL):
    """calc dist
    :type firstL: list
    :type secondL: list
    :rtype : float
    """
    geom_first = geomean(firstL)
    geom_second = geomean(secondL)
    diste = 0
    for i in range(len(firstL)):
        diste += math.pow((math.log((firstL[i]/geom_first))-math.log((secondL[i]/geom_second))), 2)
    return pow(diste, 0.5)


def build_centre(gene, condi):
    """build centre
    :type gene: Gene
    :type condi: str
    :rtype : list
    """
    center_list = []
    for i in range(len(gene.junctionL)):
        column_values = []
        for sample in gene.conditions[condi]:
            column_values.append(gene.conditions[condi][sample][i])
        if len(column_values) == 0:
            sys.stderr.write(" ERROR empty list to geom mean \n")
            sys.exit(1)
        center_list.append(geomean(column_values))
    geom_sum = sum(center_list)
    center_list_result = [x / geom_sum for x in center_list]
    return center_list_result


def dist_centre_2_centre(gene):
    """build centre
    :type gene: Gene
    """
    for condis in gene.conditions:
        gene.centreH[condis] = build_centre(gene, condis)
    myDist = dist_simplex(gene.centreH[baseC], gene.centreH[compC])
    gene.centreDist = myDist


def dist_all_sample_2_centre_cluster(gene):
    """build centre
    :type gene: Gene
    :rtype : float
    """
    distance = 0
    centre1 = build_centre_clr_transformed_cluster(gene)
    for sample in gene.sampleH:
        if sample in gene.sampleH:
            distance += dist_eukledian(centre1, gene.sampleH[sample])
    return distance/len(gene.sampleH)


def printDists(genes):
    """print distances
    :type gene: Gene
    """
    sys.stdout.write(genes.ensg+"\t"+str(genes.centreDist)+"\n")
    sys.stdout.write("centre_"+str(baseC)+":\t"+str(genes.centreH[baseC])+"\n")
    sys.stdout.write("centre_"+str(compC)+":\t"+str(genes.centreH[compC])+"\n")


def dist_centre_2_centre_clr_transformed(gene):
    """calc
    :type gene: Gene
    """
    for condis in gene.conditions:
        gene.centreH_clr[condis] = build_centre_clr_transformed(gene, condis)
    myDist = dist_eukledian(gene.centreH_clr[baseC], gene.centreH_clr[compC])
    gene.centreDist = myDist


def benjaminihochberg_correction(pvalues):
    """calc q values
    :type pvalues: list
    :rtype new_pvalues: list
    """
    n = float(len(pvalues))
    new_pvalues = list(range(int(n)))
    values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
    values.sort(reverse=True)
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


def test_juctions(gene):
    """calc
    :type gene: Gene
    """
    for indx, junction in enumerate(gene.junctionL):
        junction.ab_change = gene.centreH_clr[baseC][indx]-gene.centreH_clr[compC][indx]
        if gene.centreH_clr[baseC][indx] >= gene.centreH_clr[compC][indx]:
            junction.bio_sig = gene.centreH_clr[baseC][indx]
        else:
            junction.bio_sig = gene.centreH_clr[compC][indx]
        if abs(junction.ab_change) >= myfc:
            condition_base = []
            condition_compare = []
            for sample in gene.conditions[baseC]:
                condition_base.append(gene.conditions[baseC][sample][indx])
            for sample in gene.conditions[compC]:
                condition_compare.append(gene.conditions[compC][sample][indx])
            p_value = ranksums(condition_base, condition_compare)[1]
            p_valueL.append(p_value)
            gene.junctionL[indx].p_value = p_value
            myTuple = (gene.ensg, indx)
            correspondingL.append(myTuple)


def printJunctions(genes):
    """calc
    :type gene: Gene
    """
    sig_juncs = []
    for junction in genes.junctionL:
        if junction.significant == "yes":
            sig_juncs.append(junction)
    for junction in genes.junctionL:
        sys.stdout.write(str(genes.chrom)+":"+str(junction.left)+"-"+str(junction.right)+"\t")
        sys.stdout.write(junction.jtype+"\t")
        sys.stdout.write(str(junction.ab_change)+"\t")
        sys.stdout.write(str(junction.p_value)+"\t"+str(junction.q_value)+"\t")
        sys.stdout.write(genes.ensg+"\t"+genes.name+"\t"+str(len(genes.junctionL))+"\t"+str(len(sig_juncs))+"\t"+str(genes.centreDist)+"\t")
        sys.stdout.write(junction.significant+"\t"+junction.zero_replaced+"\n")


def printMatrix(gene):
    """calc
    :type gene: Gene
    """
    sys.stderr.write(gene.ensg+"\n")
    sys.stderr.write("gene _junction_laenge:: "+str(str(len(gene.junctionL)))+"\n")
    sys.stderr.write("gene_condition_laenge:: "+str(str(len(gene.conditions)))+"\n")
    for condi in gene.conditions:
        sys.stderr.write("condi: "+str(condi)+"\n")
        for sample in gene.conditions[condi]:
            sys.stderr.write(sample+"\t"+str(gene.conditions[condi][sample])+"\n")


def calculate_pairwise_dist_matrix():
    """calc  pairwise_dist_matrix
    """
    amount_of_genes = len(geneH)
    samplePairL = list(itertools.combinations(sampleList, 2))
    cond_dist_matrix = np.zeros((len(samplePairL),), dtype=np.float)
    for idx, samplePair in enumerate(samplePairL):
        substractCounter = 0
        samplePairDist = 0
        for gene in geneH.values():
            if samplePair[0] in gene.sampleH and samplePair[1] in gene.sampleH:
                samplePairDist += dist_eukledian(gene.sampleH[samplePair[0]], gene.sampleH[samplePair[1]])
            else:
                substractCounter += 1
        cond_dist_matrix[idx] = samplePairDist/(amount_of_genes-substractCounter)
    linkageResult = linkage(cond_dist_matrix, method='average')
    if useClusterNames:
        #dendrogram(linkageResult, color_threshold=3.46, leaf_font_size=4,labels = sampleListCluster)
        dendrogram(linkageResult, leaf_font_size=4, labels=sampleListCluster)
    else:
        dendrogram(linkageResult, leaf_font_size=4, labels=sampleList)
    plt.savefig(dendro_name+".pdf")


def find_interesting_genes():
    """find_interesting_genes

    """
    genNameVarL = []
    for gene in geneH.values():
        genNameVarL.append((gene.ensg, dist_all_sample_2_centre_cluster(gene)))
    genNameVarL.sort(key=operator.itemgetter(1))
    for ele in genNameVarL[:int(len(genNameVarL) * 0.75)]:
        del geneH[ele[0]]


def perform_clustering():
    """perform clustering

    """
    #import matplotlib
    #matplotlib.use('Agg')
    #import matplotlib.pyplot as plt
    
    nbinnormalL = nbinom.rvs(1, 0.4, loc=1, size=1000000, random_state=my_seed)
    toRemoveG = []
    read_table_cluster()
    sys.stderr.write("Number of genes before matrix reduction: "+str(len(geneH))+"\n")
    for gene in geneH.values():
        if reduceMatrix_cluster(gene, nbinnormalL):
            toRemoveG.append(gene.ensg)
    for ele in toRemoveG:
        del geneH[ele]
    sys.stderr.write("Number of genes remaining: "+str(len(geneH))+"\n")
    if len(geneH) < 1:
        sys.exit(1)
    for gene in geneH.values():
        transform_cluster(gene)
    find_interesting_genes()
    calculate_pairwise_dist_matrix()


def perform_DAS_detection():
    """run algo

    """
    sig_counter = 0
    nbinnormalL = nbinom.rvs(1, 0.4, loc=1, size=1000000, random_state=my_seed)
    toRemoveG = []
    read_table()
    sys.stderr.write("Number of genes before matrix reduction: "+str(len(geneH))+"\n")
    for gene in geneH.values():
        if reduceMatrix(gene, nbinnormalL):
            toRemoveG.append(gene.ensg)
    for ele in toRemoveG:
        del geneH[ele]
    sys.stderr.write("Number of genes remaining: "+str(len(geneH))+"\n")
    if len(geneH) < 1:
        sys.exit(1)
    for gene in geneH.values():
        transform(gene)
        dist_centre_2_centre_clr_transformed(gene)
        test_juctions(gene)
    q_valueL = benjaminihochberg_correction(p_valueL)
    for idx, ele in enumerate(correspondingL):
        geneH[ele[0]].junctionL[ele[1]].q_value = q_valueL[idx]
        if q_valueL[idx] < myqvalue:
            geneH[ele[0]].junctionL[ele[1]].significant = "yes"
            sig_counter += 1
    sys.stdout.write("junction\tjunction_type\tabundance_change\t")
    sys.stdout.write("p_val\tq_val\tgeneID\tgeneName\tjunctions\tsig_junctions\tcentre_dist\tsignificant\tzero_replacement\n")
    for gene in geneH.values():
        printJunctions(gene)
    sys.stderr.write('Number of significant differential alternative splicing events: '+str(sig_counter)+'\n')


def run_tool():
    """run tool

    """
    if perform_cluster:
        perform_clustering()
    else:
        perform_DAS_detection()
    sys.stderr.write('Have fun with your results. \n')
    sys.stdout.flush()
    sys.stdout.close()


run_tool()
