#!/usr/bin/env python3

import sys, argparse, os, inspect, gzip, glob, re, logging
import traceback as tb

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='Creates contrast list for DIEGO')
    parser.add_argument('-g', '--groupfile', type=str, default="", help='')
    parser.add_argument('-c', '--comparisons', type=str, default="", help='')
    parser.add_argument('-o', '--outdir', type=str, default="", help='')
    parser.add_argument('--loglevel', default='INFO', help="Log verbosity" )

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def create_tables(groupfile,comparisons,outdir):
    logid = scriptname+'.create_tables: '

    comps = comparisons.split(",")

    sample_dict = {}
    with open(groupfile, "r") as gf:
        for line in gf:
            l = line.replace("\n","")
            if l == "":
                continue
            sample_dict[l.split("\t")[0]] = [i for i in l.split("\t")[1].split("|")]

    log.debug(logid+'COMPS: '+str(comps))

    for c in comps:
        contrast_name = c.split(":")[0]
        log.debug('C: '+str(c) + '\tcn: '+str(contrast_name))
        contrast_group1 = [i for i in c.split(":")[1].split("-vs-")[0].split("+")]
        contrast_group2 = [i for i in c.split(":")[1].split("-vs-")[1].split("+")]


        outstring = ""
        for condition in contrast_group1:
            for sample in sample_dict[condition]:
                outstring += f"{contrast_name}_1\t{sample}\n"
        for condition in contrast_group2:
            for sample in sample_dict[condition]:
                outstring += f"{contrast_name}_2\t{sample}\n"

        with open(f"{outdir}{contrast_name}_contrast.txt","w") as outfile:
            outfile.write(outstring)


####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        print(scriptname)
        print(os.path.basename(inspect.stack()[-1].filename))
        try:
            makelogdir('LOGS')
            log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
            log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog
        except:
            log = logging.getLogger(os.path.basename(inspect.stack()[-1].filename))

        log.info('RUNNING!')
        create_tables(args.groupfile, args.comparisons, args.outdir)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


# diego_contrast_files.py ends here