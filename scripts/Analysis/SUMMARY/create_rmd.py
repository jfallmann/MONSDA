import os, sys, re
import argparse
import subprocess
from lib.Logger import *

log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M')
log = setup_logger(name=scriptname, log_file='stderr', logformat='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M')

from lib.Collection import *

def parseargs():
    parser = argparse.ArgumentParser(description='Wrapper around snakemake to run config based jobs automatically')
    parser.add_argument("-w", "--workflows", type=str, help='which workflows can be summarized')
    parser.add_argument("-h", "--head", type=str, help='header Rmd file')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_known_args()

def create_rmd(workflows, headfile, optionalargs=None):

    head = os.path.abspath(os.path.join('nextsnakes','scripts','SUMMARY','SUM_HEAD.Rmd'))







smko = os.path.abspath(os.path.join(subdir,'summary_subsnake.smk'))
if os.path.exists(smko):
    os.rename(smko,smko+'.bak')
with open(rmd, 'a') as smkout:
    with open(smkf,'r') as smk:
        for line in smk.readlines():
            line = re.sub(logfix,'loglevel="'+loglevel+'"',line)
            line = re.sub(condapath,'conda:  "../',line)
            smkout.write(line)
    smkout.write('\n\n')




####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.create_Rmd: '
    try:
        args=parseargs()
        knownargs=args[0]
        optionalargs=args[1:]

        log.setLevel(knownargs.loglevel)

        create_rmd(knownargs.workflows, knownargs.head, optionalargs[0])

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
