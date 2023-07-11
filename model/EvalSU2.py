from subprocess import Popen, PIPE
from shutil import rmtree
from multiprocessing import Pool, current_process
import os
import sys

sys.path.append(os.environ['SU2_RUN'])
str_filedir__ = os.path.dirname(os.path.abspath(__file__))
sys.path.append(str_filedir__)

import SU2
from util import MyLogger, safeCopyDirs, safeMakeDirs, safeSplitFileStr

if os.name == 'nt':
    mpi_command = 'mpiexec'
else:
    mpi_command = 'mpirun'


def runSU2CFD(mesh_filestr, cfg_filestr, partitions: int, CFD_config: dict = None, restart_filestr=None,
              dir_temp=os.path.dirname(os.path.abspath(__file__))+"/SU2_temp",
              description: str = None, remove_dump=False, out_logger=None):
    if out_logger is None:
        out_logger = MyLogger(os.path.join(str_filedir__, "SU2.log"))
    dir_cur = os.getcwd()

    mesh_filename, mesh_filedir, mesh_filestr = safeSplitFileStr(
        mesh_filestr)
    cfg_filename, cfg_filedir, cfg_filestr = safeSplitFileStr(
        cfg_filestr)
    if restart_filestr is not None:
        restart_filename, restart_filedir, restart_filestr = safeSplitFileStr(
            restart_filestr)

    # Config
    config = SU2.io.Config(cfg_filestr)
    config['NUMBER_PART'] = partitions

    if config['SOLVER'] == "MULTIPHYSICS":
        print("Parallel computation script not compatible with MULTIPHYSICS solver.")
        exit(1)

    # process input
    if restart_filestr is not None:
        config['RESTART_SOL'] = 'YES'
        config['RESTART_FILENAME'] = restart_filename

    config['MESH_FILENAME'] = mesh_filename
    if 'cgns' in mesh_filename:
        config['MESH_FORMAT'] = 'CGNS'
    else:
        config['MESH_FORMAT'] = 'SU2'

    # add input config into config
    if CFD_config is not None:
        config.updata(CFD_config)

    procid = int(current_process().pid)
    out_logger.logger.info("SU2 calculating ... procid={}".format(procid))
    if description is not None:
        dir_work = os.path.join(dir_temp, description+"_{}".format(procid))
    else:
        dir_work = os.path.join(dir_temp, "_{}".format(procid))

    safeMakeDirs(dir_work, out_logger)
    # copy mesh file to dir work
    safeCopyDirs((mesh_filename), mesh_filedir, dir_work, out_logger)
    if restart_filestr is not None:
        # copy restart data to dir work
        safeCopyDirs((restart_filename), restart_filedir, dir_work, out_logger)

    # State
    config.dump(os.path.join(dir_work, 'config_CFD.cfg'))
    # state = SU2.io.State()
    # out_logger.logger.info(state)
    out_logger.logger.info('mesh input file: {0}'.format(mesh_filestr))
    out_logger.logger.info('cfg file: {0}'.format(cfg_filestr))
    out_logger.logger.info('AOA: {0}, SIDESLIP_ANGLE: {1}, Ma: {2}, T: {3}, P: {4}'.format(
        config['AOA'], config['SIDESLIP_ANGLE'], config['MACH_NUMBER'], config['FREESTREAM_TEMPERATURE'], config['FREESTREAM_PRESSURE']))

    # run SU2 CFD
    # SU2.run.CFD(config)
    out_logger.logger.info("begin run SU2 CFD")
    info_file = open(os.path.join(dir_work, 'SU2_CFD_info.log'), 'w')
    run_command = mpi_command+" -np " + \
        str(config.NUMBER_PART)+" SU2_CFD config_CFD.cfg"
    sys.stdout.flush()
    process = Popen(run_command, shell=True, stdin=PIPE,
                    stdout=info_file, stderr=PIPE, cwd=dir_work)
    stdout, stderr = process.communicate()
    error_message = stderr.decode()
    info_file.close()

    retry_time = 0
    while ((error_message != '') and (retry_time < 2)):
        if error_message.find('retrying'):
            # retry again
            out_logger.logger.warning("node busy, retry run SU2 CFD")

            info_file = open(os.path.join(dir_work, 'SU2_CFD_info.log'), 'w')
            run_command = mpi_command+" -np " + \
                str(config.NUMBER_PART)+" SU2_CFD config_CFD.cfg"
            sys.stdout.flush()
            process = Popen(run_command, shell=True, stdin=PIPE,
                            stdout=info_file, stderr=PIPE, cwd=dir_work)
            stdout, stderr = process.communicate()
            error_message = stderr.decode()
            info_file.close()

            retry_time = retry_time+1
        else:
            err_file = open(os.path.join(dir_work, 'SU2_CFD_err.log'), 'w')
            err_file.write(error_message)
            err_file.close()
            out_logger.logger.error(dir_work+error_message+"\n")
            raise RuntimeError(
                "EvalSU2.runSU2: fatal error with SU2 CFD, proid={}".format(procid))

    out_logger.logger.info("end run SU2 CFD")

    # obtain result
    history_filestr = config.CONV_FILENAME+".csv"
    SU2_history = readSU2CSV(os.path.join(dir_work, history_filestr))
    surface_filestr = config.SURFACE_FILENAME+".csv"
    SU2_surface = readSU2CSV(os.path.join(dir_work, surface_filestr))

    # delete temp file
    if remove_dump:
        out_logger.logger.info("cleaning temp files...")
        rmtree(dir_work)

    return SU2_history, SU2_surface


def runSU2DEF(mesh_filestr, cfg_filestr, dat_filestr, mesh_out_filestr, partitions: int, DEF_config: dict = None,
              dir_temp=os.path.dirname(os.path.abspath(__file__))+"/SU2_temp",
              description: str = None, remove_dump=False, out_logger=None):
    if out_logger is None:
        out_logger = MyLogger(os.path.join(str_filedir__, "SU2.log"))
    dir_cur = os.getcwd()

    mesh_filename, mesh_filedir, mesh_filestr = safeSplitFileStr(
        mesh_filestr)
    cfg_filename, cfg_filedir, cfg_filestr = safeSplitFileStr(
        cfg_filestr)
    dat_filename, dat_filedir, dat_filestr = safeSplitFileStr(
        dat_filestr)
    mesh_out_filename, mesh_out_filedir, mesh_out_filestr = safeSplitFileStr(
        mesh_out_filestr)

    # Config
    config = SU2.io.Config(cfg_filestr)
    config['NUMBER_PART'] = partitions

    # process input
    config['DV_FILENAME'] = dat_filename
    config['MESH_FILENAME'] = mesh_filename
    config['MESH_OUT_FILENAME'] = mesh_out_filename

    # add input config into config
    if DEF_config is not None:
        config.updata(DEF_config)

    procid = int(current_process().pid)
    out_logger.logger.info("SU2 calculating ... procid={}".format(procid))
    if description is not None:
        dir_work = os.path.join(dir_temp, description+"_{}".format(procid))
    else:
        dir_work = os.path.join(dir_temp, "_{}".format(procid))

    safeMakeDirs(dir_work, out_logger)
    # copy mesh file to dir work
    safeCopyDirs((mesh_filename), mesh_filedir, dir_work, out_logger)
    # copy dat data to dir work
    safeCopyDirs((dat_filename), dat_filedir, dir_work, out_logger)

    # State
    config.dump(os.path.join(dir_work, 'config_DEF.cfg'))
    out_logger.logger.info('mesh input file: {0}'.format(mesh_filestr))
    out_logger.logger.info('cfg file: {0}'.format(cfg_filestr))
    out_logger.logger.info('dat file: {0}'.format(dat_filestr))
    out_logger.logger.info('mesh output file: {0}'.format(mesh_out_filestr))
    # state = SU2.io.State()
    # out_logger.logger.info(state)

    # run SU2 deform mesh
    out_logger.logger.info("begin run SU2 DEF")
    info_file = open(os.path.join(dir_work, 'SU2_DEF_info.log'), 'w')
    run_command = mpi_command+" -np " + \
        str(config.NUMBER_PART)+" SU2_DEF config_DEF.cfg"
    sys.stdout.flush()
    process = Popen(run_command, shell=True, stdin=PIPE,
                    stdout=info_file, stderr=PIPE, cwd=dir_work)
    stdout, stderr = process.communicate()
    error_message = stderr.decode()
    info_file.close()

    retry_time = 0
    while ((error_message != '') and (retry_time < 2)):
        if error_message.find('retrying'):
            # retry again
            out_logger.logger.warning("node busy, retry run SU2 DEF")
            info_file = open(os.path.join(dir_work, 'SU2_DEF_info.log'), 'w')
            run_command = mpi_command+" -np " + \
                str(config.NUMBER_PART)+" SU2_DEF config_DEF.cfg"
            sys.stdout.flush()
            process = Popen(run_command, shell=True, stdin=PIPE,
                            stdout=info_file, stderr=PIPE, cwd=dir_work)
            stdout, stderr = process.communicate()
            error_message = stderr.decode()
            info_file.close()
            retry_time = retry_time+1
        else:
            err_file = open(os.path.join(dir_work, 'SU2_DEF_err.log'), 'w')
            err_file.write(error_message)
            err_file.close()
            out_logger.logger.error(dir_work+error_message+"\n")
            raise RuntimeError(
                "EvalSU2.runSU2DEF: fatal error with SU2, proid={}".format(procid))

    out_logger.logger.info("end run SU2 DEF")

    # move out
    safeCopyDirs((config['MESH_OUT_FILENAME']),
                 dir_work, mesh_out_filedir, out_logger)

    # delete temp file
    if remove_dump:
        out_logger.logger.info("cleaning temp files...")
        rmtree(dir_work)

    return


def runSU2CFDParallel(model_num: int, parallel_num: int, mesh_filestr_list, cfg_filestr_list, partitions: int, CFD_config_list: dict = None,
                      dir_temp=os.path.dirname(
                          os.path.abspath(__file__))+"/SU2_temp",
                      description_list: list = "", remove_dump=False, out_logger=None):
    if out_logger is None:
        out_logger = MyLogger(os.path.join(str_filedir__, "SU2.log"))
    if parallel_num > model_num:
        parallel_num = model_num

    cur_flag = 0  # current complete task number

    # process input parameter
    mesh_filestr_list = prePar(mesh_filestr_list, model_num, out_logger)
    cfg_filestr_list = prePar(cfg_filestr_list, model_num, out_logger)
    CFD_config_list = prePar(CFD_config_list, model_num, out_logger)
    description_list = prePar(description_list, model_num, out_logger)

    # parallel analysis
    out_logger.logger.info("parallel run SU2 begin")
    SU2_out_list = []
    SU2_out_histroy = []

    pool = Pool(processes=parallel_num)
    out_list = []
    for model_index in range(model_num):
        # generate input
        args = (mesh_filestr_list[model_index],
                cfg_filestr_list[model_index], partitions,
                CFD_config_list[model_index],
                out_logger, dir_temp, description_list[model_index], remove_dump)

        out = pool.apply_async(runSU2CFD, args)
        out_list.append(out)

    pool.close()
    pool.join()

    # get result
    for out in out_list:
        SU2_out, SU2_history = out.get()
        SU2_out_list.append(SU2_out)
        SU2_out_histroy.append(SU2_history)

    out_logger.logger.info("parallel run SU2 done")

    return SU2_out_list, SU2_out_histroy


def runSU2CFDSerial(model_num: int, mesh_filestr_list, cfg_filestr_list, partitions: int, CFD_config_list: dict = None,
                    dir_temp=os.path.dirname(
                        os.path.abspath(__file__))+"/SU2_temp",
                    description_list: list = "", remove_dump=False, out_logger=None):
    if out_logger is None:
        out_logger = MyLogger(os.path.join(str_filedir__, "SU2.log"))
    # process input parameter
    mesh_filestr_list = prePar(mesh_filestr_list, model_num, out_logger)
    cfg_filestr_list = prePar(cfg_filestr_list, model_num, out_logger)
    CFD_config_list = prePar(CFD_config_list, model_num, out_logger)
    description_list = prePar(description_list, model_num, out_logger)

    # serial analysis
    out_logger.logger.info("serial run SU2 begin")
    SU2_out_list = []
    SU2_out_histroy = []

    for model_index in range(model_num):
        SU2_out, SU2_history = runSU2CFD(mesh_filestr_list[model_index],
                                         cfg_filestr_list[model_index], partitions,
                                         CFD_config_list[model_index],
                                         out_logger, dir_temp, description_list[model_index], remove_dump)

        SU2_out_list.append(SU2_out)
        SU2_out_histroy.append(SU2_history)

    out_logger.logger.info("serial run SU2 done")

    return SU2_out_list, SU2_out_histroy


def prePar(par, num, out_logger):
    """extend parameter to number

    Args:
        par (_type_): _description_
        num (_type_): _description_
        out_logger (_type_): _description_

    Raises:
        RuntimeError: _description_

    Returns:
        _type_: _description_
    """
    if not isinstance(par, list):
        par = [par]
    if len(par) == 1:
        par = par*num
    if len(par) != num:
        str_error = "FvalSHABP.analysisSHABPParallel: number of input parameter no equal to model_number"
        out_logger.logger.error(str_error)
        raise RuntimeError(str_error)
    return par


def readSU2CSV(filestr):
    """read analysis data of SU2

    Args:
        filestr (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    type_list=[]
    value_list=[]
    
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

    value_list=[[value_list[j][i] for j in range(len(value_list))] for i in range(len(value_list[0]))]
    
    # make dictionary
    data = dict(zip(type_list, value_list))

    return data


if __name__ == "__main__":
    import pickle
    import scipy.io as scio

    # out_logger=MyLogger("analysis_SU2.log")

    # mesh_filestr = "test_model/mesh_NACA0012_inv.su2"
    # cfg_filestr = "test_model/inv_NACA0012.cfg"

    # AOA = 5
    # SU2_out, SU2_history = runSU2(mesh_filestr, cfg_filestr, 1, AOA=AOA,description='test')
    # print(SU2_out)

    # mesh_filestr = "test_model/mesh_NACA0012_inv.su2"
    # cfg_filestr = "test_model/inv_NACA0012.cfg"

    # AOA_list = [0, 2, 4]

    # SU2_out, SU2_history = runSU2CFDParallel(
    #     3, 2, mesh_filestr, cfg_filestr, 2, AOA_list=AOA_list)

    # SU2_out, SU2_history = runSU2Serial(
    #     3, mesh_filestr, cfg_filestr, 2, AOA_list=AOA_list)

    # with open("result.pkl",'wb') as output:
    #     pickle.dump(SU2_out,output)

    # CL_list=[SU2_out_unit['CL'] for SU2_out_unit in SU2_out]
    # CD_list=[SU2_out_unit['CD'] for SU2_out_unit in SU2_out]

    # scio.savemat('data.mat',{'CL_list':CL_list,'CD_list':CD_list})

    # SU2_out, SU2_history = readHistory('history.csv')
    # print(list(SU2_out.keys()))

    SU2_surface = readSU2CSV('model/history.csv')
    scio.savemat('data.mat', {'SU2_surface': SU2_surface})
