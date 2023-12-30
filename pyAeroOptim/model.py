from subprocess import Popen, PIPE
from shutil import rmtree, copy2
from multiprocessing import Pool, current_process
import os
import sys
import time
import logging


if os.name == 'nt':
    mpi_command = 'mpiexec'
else:
    mpi_command = 'mpirun'


class MyLogger():
    '''Create a logger object
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


def safeMakeDirs(dir, out_logger=None):
    '''safety make dir, error will be recorf if cannot create dir

    Args:
        dir (str): File dirtory
        out_logger (logger): Logger to record error
    '''
    try:
        os.makedirs(dir)
    except OSError:
        if (not os.path.isdir(dir)) and (out_logger is not None):
            out_logger.logger.error(
                "could not create directory at {}".format(dir))
        return 1
    return 0


def safeCopyDirs(filename, old_dir, new_dir, out_logger):
    """Safe copy file from old dir to new dir

    Args:
        filename (str): File name
        old_dir (str): old file dirtory
        new_dir (str): New file dirtory
        out_logger (logger): Logger to record error
    """
    safeMakeDirs(new_dir, out_logger)
    if type(filename) is str:
        filename = (filename, )
    for f in filename:
        if os.path.normcase(old_dir) != os.path.normcase(new_dir):
            try:
                copy2(os.path.join(old_dir, f),
                      os.path.join(new_dir, f))
            except FileNotFoundError as e:
                if out_logger is not None:
                    out_logger.logger.error(e)
                return 1
    return 0


def safeSplitFileStr(filestr):
    """Safety split str file to file name and file dir

    Args:
        filestr (str): Filestr

    Returns:
        filename: File name
        filedir: File dirtory
        filestr: File string
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

    return filename, filedir, filestr


def readSU2CSV(filestr):
    """Read analysis data of SU2

    Args:
        filestr (str): File string of data csv

    Returns:
        data (dict): Dict for data of 'filestr'.csv file.
    """

    type_list = []
    value_list = []

    if not os.path.isfile(filestr):
        return None

    # read out file data
    with open(filestr, "r") as file_data:
        file_data.seek(0, 2)
        eof = file_data.tell()  # end of file
        file_data.seek(0, 0)

        # read head
        str_type = file_data.readline().replace(
            " ", "").replace('"', "").replace('[', "").replace(']', "").replace('\n', "")
        type_list = str_type.split(",")

        while file_data.tell() < eof:
            str_data = file_data.readline().replace(" ", "")
            data = str_data.split(",")
            value_list.append([float(data_unit) for data_unit in data])

    value_list = [[value_list[j][i] for j in range(
        len(value_list))] for i in range(len(value_list[0]))]

    # make dictionary
    data = dict(zip(type_list, value_list))

    return data


def runSU2CFD(cfg_param, partitions: int, dir_temp=os.path.join(os.getcwd(), 'SU2_temp'),
              run_desc: str = None, REMOVE_TEMP=False, out_logger=None):
    '''
    Interface of SU2_CFD
    Base on input cfg_param and partitions to run SU2_CFD

    Args:
        cfg_param: Cfg_filestr or class of 'SU2.io.Config'.
        partitions: Processes number of parallel run SU2_CFD

    Returns:
        SU2_data: Dict for ouput data of 'CONV_FILENAME'.csv file(if exist).
        SU2_history: Dict for history data of 'CONV_FILENAME'.csv file(if exist).
        SU2_surface: Dict for surface data of 'SURFACE_FILENAME'.csv file(if exist).
        SU2_CFD_info: Shell output of SU2_CFD.
    '''

    sys.path.append(os.environ['SU2_RUN'])
    import SU2

    dir_cur = os.getcwd()

    # Config
    if isinstance(cfg_param, dict):
        config = cfg_param
        cfg_filestr = cfg_param._filename
    elif isinstance(cfg_param, str):
        cfg_filestr = cfg_param
        cfg_filename, cfg_filedir, cfg_filestr = safeSplitFileStr(cfg_filestr)
        config = SU2.io.Config(cfg_param)
    else:
        raise RuntimeError(
            'runSU2CFD: error input, config is not dict or cfg_filestr')
    config['NUMBER_PART'] = partitions

    if config['SOLVER'] == 'MULTIPHYSICS':
        print("Parallel computation script not compatible with MULTIPHYSICS solver.")
        exit(1)

    # create dir work
    procid = int(current_process().pid)
    if run_desc is not None:
        dir_work = os.path.join(dir_temp, run_desc)
    else:
        dir_work = os.path.join(dir_temp, "PID_{}".format(procid))
    if out_logger is not None:
        out_logger.logger.info("SU2 calculating ... procid={}".format(procid))

    # create dir work
    safeMakeDirs(dir_work, out_logger)

    # process mesh file
    if 'MESH_FILENAME' not in config.keys():
        raise RuntimeError('runSU2CFD: config lack MESH_FILENAME define')
    mesh_filestr = config['MESH_FILENAME']
    mesh_filename, mesh_filedir, mesh_filestr = safeSplitFileStr(mesh_filestr)
    if not os.path.isfile(mesh_filestr):
        raise RuntimeError('runSU2CFD: mesh file do not exist')
    config['MESH_FILENAME'] = mesh_filename
    if 'cgns' in mesh_filename:
        config['MESH_FORMAT'] = 'CGNS'
    elif 'su2' in mesh_filename:
        config['MESH_FORMAT'] = 'SU2'
    else:
        raise RuntimeError('runSU2CFD: unsupport mesh type')
    safeCopyDirs((mesh_filename), mesh_filedir, dir_work, out_logger)

    # process restart file
    if 'RESTART_FILENAME' in config:
        restart_filestr = config['RESTART_FILENAME']
        restart_filename, restart_filedir, restart_filestr = safeSplitFileStr(
            restart_filestr)
        if os.path.isfile(restart_filestr):
            config['RESTART_SOL'] = 'YES'
            config['RESTART_FILENAME'] = restart_filename
            safeCopyDirs((restart_filename), restart_filedir,
                         dir_work, out_logger)

    # State
    config.dump(os.path.join(dir_work, 'config_CFD.cfg'))
    # state = SU2.io.State()
    # out_logger.logger.info(state)
    if out_logger is not None:
        out_logger.logger.info('mesh input file: {0}'.format(mesh_filestr))
        out_logger.logger.info('cfg file: {0}'.format(cfg_filestr))
        out_logger.logger.info('AOA: {0}, SIDESLIP_ANGLE: {1}, Ma: {2}, T: {3}, P: {4}'.format(
            config['AOA'], config['SIDESLIP_ANGLE'], config['MACH_NUMBER'], config['FREESTREAM_TEMPERATURE'], config['FREESTREAM_PRESSURE']))

    # run SU2 CFD
    # SU2.run.CFD(config)
    if out_logger is not None:
        out_logger.logger.info("begin run SU2 CFD")

    info_file = open(os.path.join(dir_work, 'SU2_CFD_info.log'), 'w')
    run_command = mpi_command+" -np " + \
        str(config.NUMBER_PART)+" SU2_CFD config_CFD.cfg"
    sys.stdout.flush()
    process = Popen(run_command, shell=True, stdin=PIPE,
                    stdout=info_file, stderr=PIPE, cwd=dir_work)
    stdout, stderr = process.communicate()
    return_code = process.returncode
    error_message = stderr.decode()
    info_file.close()

    retry_time = 0
    while ((return_code != 0) and (retry_time < 2)):
        if 'retrying' in error_message:
            # retry again
            if out_logger is not None:
                out_logger.logger.warning("node busy, retry run SU2 CFD")

            info_file = open(os.path.join(dir_work, 'SU2_CFD_info.log'), 'w')
            run_command = mpi_command+" -np " + \
                str(config.NUMBER_PART)+" SU2_CFD config_CFD.cfg"
            sys.stdout.flush()
            process = Popen(run_command, shell=True, stdin=PIPE,
                            stdout=info_file, stderr=PIPE, cwd=dir_work)
            stdout, stderr = process.communicate()
            return_code = process.returncode
            error_message = stderr.decode()
            info_file.close()

            retry_time = retry_time+1
        else:
            err_file = open(os.path.join(dir_work, 'SU2_CFD_err.log'), 'w')
            err_file.write(error_message)
            err_file.close()
            if out_logger is not None:
                out_logger.logger.error(dir_work+error_message+'\n')
            raise RuntimeError(
                "EvalSU2.runSU2: fatal error with SU2 CFD, proid={}".format(procid))
    info_file = open(os.path.join(dir_work, 'SU2_CFD_info.log'), 'r')
    SU2_CFD_info = info_file.read()
    info_file.close()

    if out_logger is not None:
        out_logger.logger.info("end run SU2 CFD")

    # obtain result
    history_filestr = os.path.join(dir_work,config['CONV_FILENAME']+".csv")
    if os.path.isfile(history_filestr):
        SU2_history = readSU2CSV(os.path.join(dir_work, history_filestr))
        type_list = list(SU2_history.keys())
        value_list = list(SU2_history.values())
        value_list = [value_list[idx][-1] for idx in range(len(value_list))]
        SU2_data = dict(zip(type_list, value_list))
    else:
        SU2_history = None
        SU2_data = None

    surface_filestr = os.path.join(dir_work,config['SURFACE_FILENAME']+".csv")
    if os.path.isfile(surface_filestr):
        SU2_surface = readSU2CSV(os.path.join(dir_work, surface_filestr))
    else:
        SU2_surface = None

    # delete temp file
    if REMOVE_TEMP:
        if out_logger is not None:
            out_logger.logger.info(
                'cleaning temp directory and files: {}'.format(dir_work))
        rmtree(dir_work)

    return SU2_data, SU2_history, SU2_surface, SU2_CFD_info


def runSU2DEF(cfg_param, partitions: int, dir_temp=os.path.join(os.getcwd(), 'SU2_temp'),
              run_desc: str = None, REMOVE_TEMP=False, out_logger=None):
    '''
    Interface of SU2_DEF
    Base on input cfg_param and partitions to run SU2_DEF

    Args:
        cfg_param: Cfg_filestr or class of 'SU2.io.Config'.
        partitions: Processes number of parallel run SU2_DEF

    Returns:
        mesh_out_filestr: Filestr of output mesh
        SU2_DEF_info: Shell output of SU2_DEF.
    '''

    sys.path.append(os.environ['SU2_RUN'])
    import SU2

    dir_cur = os.getcwd()

    # Config
    if isinstance(cfg_param, dict):
        config = cfg_param
        cfg_filestr = cfg_param._filename
    elif isinstance(cfg_param, str):
        cfg_filestr = cfg_param
        cfg_filename, cfg_filedir, cfg_filestr = safeSplitFileStr(cfg_filestr)
        config = SU2.io.Config(cfg_param)
    else:
        raise RuntimeError(
            'runSU2CFD: error input, config is not dict or cfg_filestr')
    config['NUMBER_PART'] = partitions

    # create dir work
    procid = int(current_process().pid)
    if run_desc is not None:
        dir_work = os.path.join(dir_temp, run_desc)
    else:
        dir_work = os.path.join(dir_temp, "PID_{}".format(procid))
    if out_logger is not None:
        out_logger.logger.info("SU2 calculating ... procid={}".format(procid))

    # create dir work
    safeMakeDirs(dir_work, out_logger)

    # process mesh file
    if 'MESH_FILENAME' not in config.keys():
        raise RuntimeError('runSU2CFD: config lack MESH_FILENAME define')
    mesh_filestr = config['MESH_FILENAME']
    mesh_filename, mesh_filedir, mesh_filestr = safeSplitFileStr(mesh_filestr)
    if not os.path.isfile(mesh_filestr):
        raise RuntimeError('runSU2CFD: mesh file do not exist')
    config['MESH_FILENAME'] = mesh_filename
    if 'cgns' in mesh_filename:
        config['MESH_FORMAT'] = 'CGNS'
    elif 'su2' in mesh_filename:
        config['MESH_FORMAT'] = 'SU2'
    else:
        raise RuntimeError('runSU2DEF: unsupport mesh type')
    safeCopyDirs((mesh_filename), mesh_filedir, dir_work, out_logger)

    # process dat file
    if 'DV_FILENAME' in config:
        dat_filestr = config['DV_FILENAME']
        dat_filename, dat_filedir, dat_filestr = safeSplitFileStr(dat_filestr)
        if os.path.isfile(dat_filestr):
            config['DV_FILENAME'] = dat_filename
            safeCopyDirs((dat_filename), dat_filedir, dir_work, out_logger)

    # process mesh out file
    if 'MESH_OUT_FILENAME' in config:
        mesh_out_filestr = os.path.join(dir_cur, config['MESH_OUT_FILENAME'])
    else:
        mesh_out_filestr = os.path.join(dir_cur, 'mesh_out.su2')
    mesh_out_filename, mesh_out_filedir, mesh_out_filestr = safeSplitFileStr(
        mesh_out_filestr)
    config['MESH_OUT_FILENAME'] = mesh_out_filename

    # State
    config.dump(os.path.join(dir_work, 'config_DEF.cfg'))
    if out_logger is not None:
        out_logger.logger.info('mesh input file: {0}'.format(mesh_filestr))
        out_logger.logger.info('cfg file: {0}'.format(cfg_filestr))
        if dat_filestr is not None:
            out_logger.logger.info('dat file: {0}'.format(dat_filestr))
        out_logger.logger.info('mesh out file: {0}'.format(mesh_out_filestr))
    # state = SU2.io.State()
    # out_logger.logger.info(state)

    # run SU2 DEF
    # SU2.run.DEF(config)
    if out_logger is not None:
        out_logger.logger.info("begin run SU2 DEF")

    info_file = open(os.path.join(dir_work, 'SU2_DEF_info.log'), 'w')
    run_command = mpi_command+" -np " + \
        str(config.NUMBER_PART)+" SU2_DEF config_DEF.cfg"
    sys.stdout.flush()
    process = Popen(run_command, shell=True, stdin=PIPE,
                    stdout=info_file, stderr=PIPE, cwd=dir_work)
    stdout, stderr = process.communicate()
    return_code = process.returncode
    error_message = stderr.decode()
    info_file.close()

    retry_time = 0
    while ((return_code != 0) and (retry_time < 2)):
        if 'retrying' in error_message:
            # retry again
            if out_logger is not None:
                out_logger.logger.warning("node busy, retry run SU2 DEF")

            info_file = open(os.path.join(dir_work, 'SU2_DEF_info.log'), 'w')
            run_command = mpi_command+" -np " + \
                str(config.NUMBER_PART)+" SU2_DEF config_DEF.cfg"
            sys.stdout.flush()
            process = Popen(run_command, shell=True, stdin=PIPE,
                            stdout=info_file, stderr=PIPE, cwd=dir_work)
            stdout, stderr = process.communicate()
            return_code = process.returncode
            error_message = stderr.decode()
            info_file.close()

            retry_time = retry_time+1
        else:
            err_file = open(os.path.join(dir_work, 'SU2_DEF_err.log'), 'w')
            err_file.write(error_message)
            err_file.close()
            if out_logger is not None:
                out_logger.logger.error(dir_work+error_message+"\n")
            raise RuntimeError(
                "EvalSU2.runSU2DEF: fatal error with SU2, proid={}".format(procid))
    info_file = open(os.path.join(dir_work, 'SU2_DEF_info.log'), 'r')
    SU2_DEF_info = info_file.read()
    info_file.close()

    if out_logger is not None:
        out_logger.logger.info("end run SU2 DEF")

    # move out
    safeCopyDirs(mesh_out_filename, dir_work, mesh_out_filedir, out_logger)

    # delete temp file
    if REMOVE_TEMP:
        if out_logger is not None:
            out_logger.logger.info(
                'cleaning temp directory and files: {}'.format(dir_work))
        rmtree(dir_work)

    return mesh_out_filestr, SU2_DEF_info


def getAtmosEnv(Z):
    '''
    base on altitude calculate atmosphere parameter
    Z is altitude m
    return temperature pressure density speed_of_sound acceleration_of_gravity
    '''
    from math import sqrt, exp
    Z = Z/1e3

    R_earth = 6.356766e3  # earth ratio km
    T_SL = 2.8815e2       # sea level temperature
    P_SL = 1.01325e5      # sea level pressure
    rou_SL = 1.2250       # sea level density
    g_SL = 9.80665        # sea level acceleration of gravity
    a_SL = 3.40294e2      # sea level speed of sound

    # geopotential height H
    H = Z/(1+Z/R_earth)

    # calculate parameter
    if Z <= 11.0191:
        W = 1-H/44.3308
        T = W*T_SL
        P = W ** 5.2553*P_SL
        rho = W ** 4.2559*rou_SL

    elif Z <= 20.0631:
        W = exp((14.9647-H)/6.3416)
        T = 216.650
        P = 1.1953*1e-1*W*P_SL
        rho = 1.5898*1e-1*W*rou_SL

    elif Z <= 32.1619:
        W = 1+(H-24.9021)/221.552
        T = 221.552*W
        P = 2.5158*1e-2*W ** -34.1629*P_SL
        rho = 3.2722*1e-2*W ** -35.1829*rou_SL

    elif Z <= 47.3501:
        W = 1+(H-39.4799)/89.4107
        T = 250.350*W
        P = 2.8338*1e-3*W ** -12.2011*P_SL
        rho = 3.2617*1e-3*W ** -13.2011*rou_SL

    elif Z <= 51.4125:
        W = exp((48.6252-H)/7.9223)
        T = 270.650
        P = 8.9155*1e-4*W*P_SL
        rho = 9.4920*1e-4*W*rou_SL

    elif Z <= 71.8020:
        W = 1-(H-59.4390)/88.2218
        T = 247.021*W
        P = 2.1671*1e-4*W ** 12.2011*P_SL
        rho = 2.5280*1e-4*W ** 11.2011*rou_SL

    elif Z <= 86.0000:
        W = 1-(H-78.0303)/100.2950
        T = 200.590*W
        P = 1.2274*1e-5*W ** 17.0816*P_SL
        rho = 1.7632*1e-5*W ** 16.0816*rou_SL

    elif Z <= 91.0000:
        W = exp((87.2848-H)/5.4700)
        T = 186.870
        P = (2.2730+1.042*1e-3*H)*1e-6*W*P_SL
        rho = 3.6411*1e-6*W*rou_SL

    else:
        T = -0.0040*Z ** 3+1.5054*Z ** 2-177.5620*Z+6.8929*1e3
        P = 2.532*1e6*exp(-0.1829*Z)+0.1403*exp(-0.03698*Z)
        rho = 70.22*1e6*exp(-0.1874*Z)+1.734*1e-5*exp(-0.05828*Z)

    a = 20.0468*sqrt(T)
    g = g_SL/(1+Z/R_earth) ** 2

    mu = 1.716*1e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4)

    return T, P, rho, a, mu, g


if __name__ == "__main__":
    # flight_cond=[5e3,7e3,10e3,15e3,20e3,25e3,30e3,40e3,50e3,60e3,70e3,80e3]
    # for cond in flight_cond:
    #     print(getAtmosEnv(cond))

    sys.path.append(os.environ['SU2_RUN'])
    import SU2
    
    config = SU2.io.Config('Airfoil/SU2/NACA0012.cfg')
    config['MESH_FILENAME'] = 'Airfoil/SU2/NACA0012.cgns'
    SU2_data, SU2_history, SU2_surface, SU2_CFD_info = runSU2CFD(config, 1)
    
    config = SU2.io.Config('Airfoil/SU2/NACA0012_deform.cfg')
    config['MESH_FILENAME'] = 'Airfoil/SU2/NACA0012.cgns'
    mesh_out_filestr, SU2_DEF_info = runSU2DEF(config, 1)
