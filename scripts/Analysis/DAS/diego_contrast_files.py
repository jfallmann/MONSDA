#! /usr/bin/env python3

import sys, argparse, os, inspect, gzip, glob, re, logging
import traceback as tb

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='diego_contrast_files.py')
    parser.add_argument('-g', '--groupfile', type=str, default="", help='')
    parser.add_argument('-c', '--comparisons', type=str, default="", help='')
    parser.add_argument('-o', '--outdir', type=str, default="", help='')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def create_tables(groupfile,comparisons,outdir):
    comps = comparisons.split(",")

    sample_dict = {}
    with open(groupfile, "r") as gf:
        for line in gf:
            l = line.replace("\n","")
            if l == "":
                continue
            sample_dict[l.split("\t")[0]] = [i for i in l.split("\t")[1].split("|")]


    for c in comps:
        contrast_name = c.split(":")[0]
        contrast_group1 = [i for i in c.split(":")[1].split("-vs-")[0].split("+")]
        contrast_group2 = [i for i in c.split(":")[1].split("-vs-")[1].split("+")]

        samples_1 = []
        for condition in contrast_group1:
            samples_1 += sample_dict[condition]
        samples_2 = []
        for condition in contrast_group2:
            samples_2 += sample_dict[condition]

        s1 = "\t".join(samples_1)
        s2 = "\t".join(samples_2)

        line1 = f"{contrast_name}_1\t{s1}"
        line2 = f"{contrast_name}_2\t{s2}"

        with open(f"{outdir}{contrast_name}_contrast.txt","w") as outfile:
            outfile.write(line1+"\n"+line2)


####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        makelogdir('LOGS')
        try:
            log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
            log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog
        except:
            log = logging.getLogger(os.path.basename(inspect.stack()[-1].filename))

        create_tables(args.groupfile, args.comparisons, args.outdir)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


# diego_contrast_files.py ends here
