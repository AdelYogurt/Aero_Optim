import os,sys
from multiprocessing import Pool, current_process
import pandas as pd
from scipy import io
from shutil import rmtree
from subprocess import Popen, PIPE
from pathlib import Path
from parea.main import _munge_triangle
from shapely.geometry import Polygon
from shapely.ops import unary_union
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from util import MyLogger, safeCopyDirs, safeMakeDirs, safeSplitStrFile

SHABP_logger = MyLogger("SHABP.log")

def runSHABP(str_inp_file, str_wgs_file=None,
                  Ma=10.0, AOA=0.0, cmethods=[], emethods=[],
                  dir_temp=os.path.dirname(os.path.abspath(__file__))+'/SHABP_temp',
                  remove_dump=True,parallel_flag:bool=False):
    """simple function to run SHABP. if do not exist inp file, funciton will generate inp file from input Ma, AOA, cmethods, emethods

    Args:
        str_inp_file (_type_): _description_
        str_wgs_file (_type_, optional): _description_. Defaults to None.
        Ma (float, optional): _description_. Defaults to 10.0.
        AOA (float, optional): _description_. Defaults to 0.0.
        cmethods (list, optional): _description_. Defaults to [].
        emethods (list, optional): _description_. Defaults to [].
        SHABP_logger (_type_, optional): _description_. Defaults to None.
        dir_temp (_type_, optional): _description_. Defaults to os.path.dirname(os.path.abspath(__file__))+'/SHABP_temp'.
        remove_dump (bool, optional): _description_. Defaults to True.
        parallel_flag (bool, optional): _description_. Defaults to False.

    Raises:
        RuntimeError: _description_

    Returns:
        SHABP_out: dict, include{'alpha':, 'beta':, 'cx':, 'cy':, 'cz':, 'cmx':, 'cmy':, 'cmz':, 'clift':, 'cdrag':, 'cpitch':}
    """
    wgs_file_name, wgs_file_dir, str_wgs_file = safeSplitStrFile(
        str_wgs_file, SHABP_logger)

    # process inp file
    if (not (isinstance(str_inp_file, str)) or not (os.path.isfile(str_inp_file)) or not (os.path.isfile(os.path.join(os.getcwd(), str_inp_file)))):
        # if inp file do not exist, create inp file in temp dir
        inp_file_name = wgs_file_name[0:-4]+'.inp'

        # support parallel mode
        if parallel_flag:
            procid = int(current_process().pid)
            SHABP_logger.logger.info("parallel make inp")
            dir_work = os.path.join(dir_temp, "{}".format(procid))

            safeMakeDirs(dir_work, SHABP_logger.logger)
            inp_file_dir = dir_work
            inp_file_name = inp_file_name
        else:
            inp_file_dir = dir_temp

        str_inp_file = os.path.join(inp_file_dir, inp_file_name)

        genINP(wgs_file_name=wgs_file_name,str_wgs_file=str_wgs_file, str_inp_file=str_inp_file, Ma=Ma, AOA=AOA,
               cmethods=cmethods, emethods=emethods,SHABP_logger=SHABP_logger)
    else:
        inp_file_name = os.path.basename(str_inp_file)
        inp_file_dir = os.path.dirname(str_inp_file)
        if not os.path.exists(inp_file_dir):
            # if inp_file_dir is empty, default inp_file sort in work_dir
            inp_file_dir = os.getcwd()
            str_inp_file = os.path.join(inp_file_dir, inp_file_name)
        # check if inp file exist
        if not (os.path.exists(str_inp_file)):
            str_error = "EvalSHABP.analysisSHABP: inp file do not exist!"
            SHABP_logger.logger.error(str_error)
            raise RuntimeError(str_error)

    # run SHABP
    procid = int(current_process().pid)
    SHABP_logger.logger.info("S/HABP calculating ... procid={}".format(procid))
    dir_work = os.path.join(dir_temp, "{}".format(procid))

    safeMakeDirs(dir_work, SHABP_logger)
    # copy wgs file to dir work
    safeCopyDirs((wgs_file_name), wgs_file_dir, dir_work, SHABP_logger)
    # copy inp file to dir work
    safeCopyDirs((inp_file_name), inp_file_dir, dir_work, SHABP_logger)

    process = Popen("hyper", stdin=PIPE, stdout=PIPE,
                    stderr=PIPE, cwd=dir_work)
    stdout, stderr = process.communicate(
        inp_file_name.encode("UTF-8"))  # notice: donot use wait()!!!
    message = stderr.decode()
    if message != '':
        SHABP_logger.logger.error("\n".join((dir_work, str(message))))
        raise RuntimeError(
            "EvalSHABP.runSHABP: fatal error with S/HABP, proid={}".format(procid))

    if os.path.isfile(os.path.join(dir_work, "hyper.err")):
        with open(os.path.join(dir_work, "hyper.err")) as f:
            if "Error termination." in f.read():
                SHABP_logger.logger.critical(
                    "\n".join((dir_work, "fatal error with S/HABP")))
                raise RuntimeError(
                    'EvalSHABP.runSHABP: fatal error with S/HABP')

    str_out_file = dir_work + "/hyper.out"
    SHABP_out = readOut(str_out_file)

    # delete temp file
    if remove_dump:
        SHABP_logger.logger.info("cleaning temp files...")
        rmtree(dir_work)

    return SHABP_out


def genINP(wgs_file_name, str_wgs_file, str_inp_file, Ma, AOA, cmethods, emethods,SHABP_logger=None,):
    """generate inp file by input

    Args:
        wgs_file_name (_type_): _description_
        str_inp_file (_type_): _description_
        Ma (_type_): _description_
        AOA (_type_): _description_
        cmethods (_type_): _description_
        emethods (_type_): _description_
        SHABP_logger (_type_, optional): _description_. Defaults to None.

    Raises:
        RuntimeError: _description_
    """
    SHABP_logger.logger.info("generating .inp file...")
    from pyPanair.preprocess import wgs_creator
    from stl import mesh  # pip install numpy-stl

    # read wgs file and generate stl file to calculate basical parameter
    wgs_this = wgs_creator.read_wgs(str_wgs_file)
    str_stl_file = str_wgs_file.replace(".wgs", ".stl")
    wgs_this.create_stl(str_stl_file)

    dictRef = {"cbar": 0, "span": 0, "sref": 0, "xref": 0}
    mesh_this = mesh.Mesh.from_file(str_stl_file)
    dictRef["cbar"] = (mesh_this.x.max()-mesh_this.x.min())
    dictRef["span"] = (mesh_this.y.max()-mesh_this.y.min())

    dictRef["xref"] = mesh_this.get_mass_properties()[1][0]
    v0, v1, v2 = mesh_this.v0[:,
                              :2], mesh_this.v1[:, :2], mesh_this.v2[:, :2]

    # calculate S_ref
    listTri = []
    for idx in range(v0.shape[0]):
        verts = _munge_triangle(
            tuple(v0[idx]), tuple(v1[idx]), tuple(v2[idx]))
        if len(verts) == 3:
            listTri.append(Polygon(verts))
    dictRef["sref"] = unary_union(listTri).area

    # method number
    network_num = len(wgs_this._networks)
    if not isinstance(cmethods, list):
        cmethods = [int(cmethods)]
    if len(cmethods) == 1:
        cmethods = cmethods*network_num
    if not isinstance(emethods, list):
        emethods = [int(emethods)]
    if len(emethods) == 1:
        emethods = emethods*network_num
    if ((len(cmethods) != network_num) | (len(emethods) != network_num)):
        str_error = "model.analysisSHABP: number of cmethods or emethods do not equal to wgs network number!"
        SHABP_logger.logger.error(str_error)
        raise RuntimeError(str_error)

    # generate inp file
    str_inp = "&hyp  title='SHABP simulation',\n  wgsFileName='"+wgs_file_name+"',\n"
    str_inp += "  cmethods=" + ", ".join(map("{}".format, cmethods)) + "\n"
    str_inp += "  emethods=" + ", ".join(map("{}".format, emethods)) + "\n"
    str_inp += "  mach={:.1f}, alpha={:.1f},\n".format(
        Ma, AOA)
    str_inp += "  sref={:.6f}, cbar={:.1f},\n".format(
        dictRef["sref"], dictRef["cbar"])
    str_inp += "  xref={:.6f}, span={:.6f}/\n".format(
        dictRef["xref"], dictRef["span"])

    # write data into inp file
    with open(str_inp_file, "w") as f:
        f.write(str_inp)
        
    return 0


def runSHABPParallel(model_num: int, parallel_num: int, str_wgs_file_list, str_inp_file_list,
                          Ma_list=None, AOA_list=None, cmethods_list=None, emethods_list=None,
                          dir_temp=os.path.dirname(os.path.abspath(__file__))+'/SHABP_temp',
                          remove_dump=True):
    """parallel run SHABP to calculate multi input

    Args:
        model_num (int): _description_
        parallel_num (int): _description_
        str_wgs_file_list (_type_): _description_
        str_inp_file_list (_type_): _description_
        Ma_list (_type_, optional): _description_. Defaults to None.
        AOA_list (_type_, optional): _description_. Defaults to None.
        cmethods_list (_type_, optional): _description_. Defaults to None.
        emethods_list (_type_, optional): _description_. Defaults to None.
        SHABP_logger (_type_, optional): _description_. Defaults to None.
        dir_temp (_type_, optional): _description_. Defaults to os.path.dirname(os.path.abspath(__file__))+'/SHABP_temp'.
        remove_dump (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    if parallel_num > model_num:
        parallel_num = model_num

    cur_flag = 0  # current complete task number

    # process input parameter
    str_wgs_file_list = prePar(str_wgs_file_list, model_num, SHABP_logger)
    str_inp_file_list = prePar(str_inp_file_list, model_num, SHABP_logger)
    Ma_list = prePar(Ma_list, model_num, SHABP_logger)
    AOA_list = prePar(AOA_list, model_num, SHABP_logger)
    cmethods_list = prePar(cmethods_list, model_num, SHABP_logger)
    emethods_list = prePar(emethods_list, model_num, SHABP_logger)

    # parallel analysis
    SHABP_out_list = []
    
    while cur_flag < model_num:
        process_number=min(parallel_num,model_num-cur_flag)

        pool=Pool(processes=process_number)
        out_list=[]
        for pro_index in range(process_number):
            # generate input
            args=(str_wgs_file_list[cur_flag+pro_index],
                               str_inp_file_list[cur_flag+pro_index],
                               Ma_list[cur_flag+pro_index],AOA_list[cur_flag+pro_index],cmethods_list[cur_flag+pro_index],emethods_list[cur_flag+pro_index],
                               SHABP_logger,dir_temp,remove_dump,True)
            out = pool.apply_async(runSHABP, args)
            out_list.append(out)

        pool.close()
        pool.join()
        
        # get result
        for out in out_list:
            SHABP_out=out.get()
            SHABP_out_list.append(SHABP_out)
        
        # add analysis nunber
        cur_flag=cur_flag+parallel_num

    return SHABP_out_list


def prePar(par, num, SHABP_logger):
    """extend parameter to number

    Args:
        par (_type_): _description_
        num (_type_): _description_
        SHABP_logger (_type_): _description_

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
        SHABP_logger.logger.error(str_error)
        raise RuntimeError(str_error)
    return par


def readOut(file):
    """read analysis data of SHABP

    Args:
        file (_type_): _description_

    Returns:
        _type_: _description_
    """
    out_body, out_wind, out_norm = list(), list(), list()
    columns = ["alpha", "beta", "cx", "cy", "cz", "cmx", "cmy", "cmz", "clift",
               "cdrag", "cpitch"]

    # read out file data
    with open(file, "r") as f:
        f.seek(0, 2)
        eof = f.tell()  # end of file
        f.seek(0, 0)

        while f.tell() < eof:
            str_data = f.readline()
            data = str_data.split()
            if not data:
                pass
            elif (str_data.find("SOLUTIONS") != -1) and (str_data.find("BODY") != -1):
                # read body direction coefficient
                str_data = f.readline()
                data = f.readline().split()
                out_body = data[1:]
                out_body = [float(out_body_unit) for out_body_unit in out_body]
            elif (str_data.find("SOLUTIONS") != -1) and (str_data.find("WIND") != -1):
                # read wind direction coefficient
                str_data = f.readline()
                data = f.readline().split()
                out_wind = data[3:]
                out_wind = [float(out_wind_unit) for out_wind_unit in out_wind]
            elif data[0] == "COMPRESSION":
                break
    
    SHABP_out = dict(zip(columns, out_body+out_wind))
    return SHABP_out


if __name__ == "__main__":
    Ma = 13.8
    AOA = 8
    str_wgs_file = 'test_model/waverider_wing_dia.wgs'
    str_inp_file = None
    cmethods = 14
    emethods = 3

    # SHABP_logger = MyLogger("analysis_SHABP.log")

    # SHABP_out = runSHABP(str_inp_file, str_wgs_file,
    #                           Ma, AOA, cmethods, emethods, SHABP_logger,parallel_flag=True)
    # print(SHABP_out)
    # SHABP_logger.logger.info("S/HABP analysis finished")
    
    # Ma = 13.8
    # AOA = 8
    # str_wgs_file = 'test_model/waverider_wing_dia.wgs'
    # str_inp_file = None
    # cmethods = 14
    # emethods = 3

    # SHABP_logger = MyLogger("analysis_SHABP.log")

    # SHABP_out = runSHABPParallel(2,2,str_inp_file,str_wgs_file,
    #                           Ma, AOA, cmethods, emethods, SHABP_logger)
    # print(SHABP_out)
    # SHABP_logger.logger.info("S/HABP analysis finished")
