# Logger.py ---
#
# Filename: Logger.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Mon Aug 12 10:26:55 2019 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Apr 29 16:42:40 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 91
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
import logging
import os
import sys
import traceback as tb


def makelogdir(logdir):
    if not os.path.isabs(logdir):
        logdir = os.path.abspath(logdir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    return logdir


def setup_logger(name, log_file, filemode='w', logformat=None, datefmt=None, level='WARNING', delay=False):
    """Function setup as many loggers as you want"""

    logger = logging.getLogger(name)
    if log_file != 'stdout' and log_file != 'stderr':
        makelogdir(os.path.dirname(log_file))
        if not os.path.isfile(os.path.abspath(log_file)):
            open(os.path.abspath(log_file), 'a').close()
        handler = logging.FileHandler(os.path.abspath(log_file), mode=filemode, delay=delay)
    else:
        handler = logging.StreamHandler()

    handler.setFormatter(logging.Formatter(fmt=logformat, datefmt=datefmt))

    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


if __name__ == '__main__':
    try:
        # set up logging to file
        log = setup_logger(
            name='',
            log_file='stderr',
            logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
            datefmt='%m-%d %H:%M',
            level='WARNING',
        )

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        logging.error(''.join(tbe.format()))


# Logger.py ends here
