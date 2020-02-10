#!/usr/bin/env python3
# Configurator.py ---
#
# Filename: Configurator.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Mon Feb 10 08:09:48 2020 (+0100)
# Version:
# Package-Requires: ()
# Last-Updated: Mon Feb 10 08:55:23 2020 (+0100)
#           By: Joerg Fallmann
#     Update #: 12
# URL:
# Doc URL:
# Keywords:
# Compatibility:
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

import glob, os, sys, inspect, json, shutil
from collections import defaultdict
import traceback as tb
from snakemake import load_configfile
from snakemake.utils import validate, min_version
import argparse
import subprocess
import re
min_version("5.8.2")

#cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib")
#if cmd_subfolder not in sys.path:
#    sys.path.insert(0, cmd_subfolder)

from lib.Collection import *
from lib.Logger import *
scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='Helper to create initial config file used for workflow processing')
    parser.add_argument("-c", "--configfile", type=str, help='Configuration json to write to')
    parser.add_argument("--pre", "--preprocess", action="store_true", choices=['SRA', 'BASECALL'], help='Which preprocessing steps to conduct')
    parser.add_argument("-w", "--workflows", action="store_true", choices=['MAPPING', 'TRIMMING', 'QC'], help='Which workflow steps to conduct')
    parser.add_argument("--post", "--postprocess", action="store_true", choices=['COUNTING','UCSC','PEAKS','ANNOTATE','DE','DEX'], help='Which workflow steps to conduct')
    parser.add_argument("-u", "--use-conda", action="store_true", default=True, help='Should conda be used')
    parser.add_argument("-b", "--binaries", type="str", help='Path to binary directory')
    parser.add_argument("-s", "--scripts", type="str", help='Path to script for execution')
    parser.add_argument("-j", "--procs", type=int, default=1, help='Maximum number of parallel processes to start snakemake with, represented by MAXTHREADS in config')
    parser.add_argument("-v", "--loglevel", type=str, default='INFO', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_known_args()

####################
#### FUNCTIONS  ####
####################



####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        knownargs=args[0]
        optionalargs=args[1:]
        makelogdir('LOGS')
        log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=knownargs.loglevel)
        log.addHandler(logging.StreamHandler(sys.stdout))  # streamlog

        MIN_PYTHON = (3,7)
        if sys.version_info < MIN_PYTHON:
            log.error("This script requires Python version >= 3.7")
            sys.exit("This script requires Python version >= 3.7")
        log.info(logid+'Running '+scriptname+' on '+str(knownargs.procs)+' cores')

        run_snakemake(knownargs.configfile, knownargs.debug_dag, knownargs.filegraph, knownargs.directory, knownargs.use_conda, knownargs.procs, knownargs.skeleton, knownargs.unlock, optionalargs[0])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))



# Configurator.py ends here
