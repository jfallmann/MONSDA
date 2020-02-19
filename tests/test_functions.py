#!/usr/bin/env python3
# test_functions.py ---
#
# Filename: test_functions.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Wed Feb 19 12:31:13 2020 (+0100)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Feb 19 13:40:24 2020 (+0100)
#           By: Joerg Fallmann
#     Update #: 42
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
import glob, os, sys, inspect, json, shutil
from collections import defaultdict
import traceback as tb
from snakemake import load_configfile
from snakemake.utils import validate, min_version
import argparse
import subprocess
import re

cmd_subfolder = [os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../snakes/lib"),os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"snakes/lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"lib")]
for x in cmd_subfolder:
    if x not in sys.path:
        sys.path.insert(0, x)

from Collection import *
from Logger import *

log = setup_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='DEBUG')
logid = 'TESTER: '

config = load_configfile('/home/fall/Work/Tests/SnakemakeTest/SubSnakes/Dm6_unpaired_std_subconfig.json')
REFERENCE=config['REFERENCE']

file=r'Dm6/unpaired/std/GSM461177_untreat_paired_subset_r1'
dir = 'Dm6'

SAMPLES=list()
SAMPLES.append(file)

print(logid+'SAMPLES: '+str(SAMPLES))

MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')
print(logid+'GENOME: '+genome(file, config))

print(logid+'OPTIONS: '+' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items()))

#print(logid+'ANNO: '+str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'MAPPING')['ANNOTATION']]))

#print(logid+'GENPATH: '+"{ref}/{dir}/{map}/{extension}/".format(ref=REFERENCE, dir=dir, map=MAPPERENV, extension=check_tool_params(SAMPLES[0], None ,config, 'MAPPING',2)))

#print(logid+)

#print(logid+)

#
# test_functions.py ends here
