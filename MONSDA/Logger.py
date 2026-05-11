# Logger.py ---


# Code:
import logging
import os
import sys
import traceback as tb


def makelogdir(logdir: str) -> str:
    """create log directory if it does not exist

    Parameters
    ----------
    logdir : str
        directory to create

    Returns
    -------
    str
        directory path
    """
    if not os.path.isabs(logdir):
        logdir = os.path.abspath(logdir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    return logdir


def setup_logger(
    name: str,
    log_file: str,
    filemode: str = "w",
    logformat: str = None,
    datefmt: str = None,
    level: str = "WARNING",
    delay: bool = False,
) -> logging.Logger:  # pragma: no cover
    """set up a logger

    Parameters
    ----------
    name : str
        name of the logger
    log_file : str
        file to log to
    filemode : str, optional
        mode, by default "w"
    logformat : str, optional
        format, by default None
    datefmt : str, optional
        dateformat, by default None
    level : str, optional
        loglevel, by default "WARNING"
    delay : bool, optional
        delay logging, by default False

    Returns
    -------
    logging.Logger
        logger instance
    """

    logger = logging.getLogger(name)
    if log_file != "stdout" and log_file != "stderr":
        makelogdir(os.path.dirname(log_file))
        if not os.path.isfile(os.path.abspath(log_file)):
            open(os.path.abspath(log_file), "a").close()
        handler = logging.FileHandler(
            os.path.abspath(log_file), mode=filemode, delay=delay
        )
    else:
        handler = logging.StreamHandler()

    handler.setFormatter(logging.Formatter(fmt=logformat, datefmt=datefmt))

    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


if __name__ == "__main__":
    try:
        # set up logging to file
        log = setup_logger(
            name="",
            log_file="stderr",
            logformat="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt="%m-%d %H:%M",
            level="WARNING",
        )

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        logging.error("".join(tbe.format()))


# Logger.py ends here
