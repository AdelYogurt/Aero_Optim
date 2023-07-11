import os
from shutil import copy2
import inspect
import time
import logging


class MyLogger():
    '''
    create a logger object
    '''
    def __init__(self, log_filename=time.asctime()+".log", loglevel=logging.INFO):
        # create logger
        self.logger = logging.getLogger(log_filename)
        self.logger.setLevel(loglevel)

        # default create a file log
        file_handler = logging.FileHandler(log_filename)
        file_handler.setLevel(loglevel)
        formatter = logging.Formatter(
            '[%(levelname)s]%(asctime)s %(filename)s:%(lineno)d: %(message)s')
        if not self.logger.handlers:
            file_handler.setFormatter(formatter)

            # add handler
            self.logger.addHandler(file_handler)
            self.logger.info("initialize logger")

    def getlog(self):
        self.logger.info("get logger")
        return self.logger

    def addStreamOut(self, loglevel=logging.INFO):
        # create a stream handler
        formatter = logging.Formatter(
            '[%(levelname)s]%(asctime)s %(filename)s:%(lineno)d: %(message)s')
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(loglevel)
        stream_handler.setFormatter(formatter)

        # add handler
        self.logger.addHandler(stream_handler)


def safeMakeDirs(dir, my_logger=None):
    '''
    safety make dir, error will be recorf if cannot create dir
    '''
    try:
        os.makedirs(dir)
    except OSError:
        if (not os.path.isdir(dir)) and (my_logger is not None):
            my_logger.logger.error(
                "could not create directory at {}".format(dir))
        return 1
    return 0


def safeCopyDirs(files, old_dir, new_dir, my_logger):
    """safe copy file from old dir to new dir

    Args:
        files (_type_): _description_
        old_dir (_type_): _description_
        new_dir (_type_): _description_
        my_logger (_type_): _description_

    Returns:
        _type_: _description_
    """
    safeMakeDirs(new_dir, my_logger)
    if type(files) is str:
        files = (files, )
    for f in files:
        if os.path.normcase(old_dir) != os.path.normcase(new_dir):
            try:
                copy2("{0}/{1}".format(old_dir, f),
                      "{0}/{1}".format(new_dir, f))
            except FileNotFoundError as e:
                my_logger.logger.error(e)
                return 1
    return 0


def safeSplitFileStr(filestr):
    """safety split str file to file name and file dir

    Args:
        filestr (_type_): _description_
        my_logger (_type_): _description_

    Raises:
        RuntimeError: _description_

    Returns:
        _type_: _description_
    """
    filename = os.path.basename(filestr)
    filedir = os.path.dirname(filestr)
    if not os.path.exists(filedir):
        # if filedir is not exist, default file sort in work_dir
        cur_dir = os.getcwd()
        filedir = os.path.join(cur_dir, filedir)
        filestr = os.path.join(cur_dir, filestr)
    if not os.path.isabs(filedir):
        cur_dir = os.getcwd()
        filedir = os.path.join(cur_dir, filedir)
        filestr = os.path.join(cur_dir, filestr)

    # # check if file exist
    # if not (os.path.exists(filestr)):
    #     # obtain call function name
    #     call_fun_name = inspect.getframeinfo(inspect.currentframe().f_back)[2]
    #     str_error = "model."+call_fun_name+": file do not exist!"
    #     my_logger.logger.error(str_error)
    #     raise RuntimeError(str_error)

    return filename, filedir, filestr
