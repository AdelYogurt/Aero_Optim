#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import os

from mpl_toolkits.mplot3d import Axes3D
from typing import Union, List, Tuple, Optional


class CST2DCurve(object):
    # generate curve by CST parameter

    def __init__(self, LX=0.0, LY=0.0, shape_par_Y=None, class_par_Y=None, symmetry_x: bool = False):
        self.tran_fun_Y = None
        self.rotation_matrix = None
        self.translation = None
        # generate 2D CST line by LX, LY, shape_par_Y, class_par_Y

        # xi, x=X
        # psi, y=Y(xi)

        # input:
        # LX,LY,shape_par_Y,class_par_Y,symmetry_x

        # notice:
        # shape_fcn_Y(XI), class_fun_Y[N1, N2]
        self.LX = LX
        self.LY = LY
        self.symmetry_x = symmetry_x

        # preprocess
        if (type(shape_par_Y) is np.ndarray) or (type(shape_par_Y) is list) and (len(shape_par_Y) == 0):
            shape_par_X = None
        if (type(class_par_Y) is np.ndarray) or (type(class_par_Y) is list) and (len(class_par_Y) == 0):
            class_par_Y = None

        if (shape_par_Y is not None):
            if (type(shape_par_Y) is np.ndarray) or (type(shape_par_Y) is list):
                # shape_par_Y only one cross-section parameter
                self.shape_par_Y = np.array(shape_par_Y)
                self.shape_fcn_Y = self.defunShape
            else:
                # function handle
                self.shape_fcn_Y = shape_par_Y
        else:
            self.shape_fcn_Y = None

        self.N1 = class_par_Y[0]
        self.N2 = class_par_Y[1]
        self.class_fun_Y = self.defunClass
        return

        # common function

    # calculate grid coordinate function
    def calCurve(self, XI=10):
        # calculate 2D CST vector by XI
        # X=CX*XI, if symmetry_x X=CX*(XI-0.5)
        # Y=CY*shape_fcn_Y.*class_fun_Y

        # default
        # xi_list=linspace(0,1,xi_gird_num+1) if symmetry_x xi_list=linspace(0.5,1,xi_gird_num+1)

        # xi, x=X
        # psi, y=Y(xi)

        # input:
        # XI
        # xi_grid_number

        # notice:
        # shape_fcn_Y(XI), class_fun_Y[N1, N2]

        # output:
        # X,Y

        # xi
        if isinstance(XI, int):
            xi_grid_numebr = XI
            XI = []
        else:
            # mean input XI matrix
            xi_grid_numebr = XI.shape[1] - 1

        if len(XI) == 0:
            if self.symmetry_x:
                XI = np.linspace(0.5, 1, xi_grid_numebr +
                                 1).reshape((-1, xi_grid_numebr + 1))
            else:
                XI = np.linspace(0, 1, xi_grid_numebr +
                                 1).reshape((-1, xi_grid_numebr + 1))

        X, Y = self.calPoint(XI)
        return X, Y

    def calPoint(self, XI):
        # calculate point on curve
        XI=np.array(XI).reshape((1,-1))

        # calculate origin surface matrix
        if self.symmetry_x:
            X = self.LX * (XI - 0.5)
        else:
            X = self.LX * XI

        if (not self.shape_fcn_Y) and (not self.class_fun_Y):
            if len(X.shape) == 1:
                Y = np.zeros(1, len(X)[0])
            else:
                Y = np.zeros((1, X.shape[1]))
        elif (not self.shape_fcn_Y):
            Y = self.LY * self.class_fun_Y(XI)
        elif (not self.class_fun_Y):
            Y = self.LY * self.shape_fcn_Y(XI)
        else:
            Y = self.LY * \
                self.shape_fcn_Y(XI) * self.class_fun_Y(XI)
        return X, Y

    def calCoordinate(self, X):
        # base on X, Y calculate local coordinate in surface

        XI = X / self.LX
        return XI

    def defunShape(self, XI):
        # default shape function

        NP = self.calNormPar(self.shape_par_Y[0], self.shape_par_Y[1])
        S = (XI ** self.shape_par_Y[0] * (1 - XI) ** self.shape_par_Y[1]) / NP
        return S

    def defunClass(self, XI):
        # default class function
        # default N1, N2 is relation with XI

        N1 = self.N1
        N2 = self.N2
        NP = self.calNormPar(N1, N2)
        XI[XI < 0] = 0
        XI[XI > 1] = 1
        C = (XI ** N1 * (1 - XI) ** N2) / NP
        return C

    def calNormPar(self, N1, N2):
        # calculate normailize class function parameter by N1, N2
        N1 = np.array(N1)
        N2 = np.array(N2)
        N1[np.logical_and(N1 == 0, N2 == 0)] = 1
        N2[np.logical_and(N1 == 0, N2 == 0)] = 1
        nomlz_par = np.array((N1 / (N1 + N2)) ** N1 * (N2 / (N1 + N2)) ** N2)
        nomlz_par[np.logical_and(N1 == 0, N2 == 0)] = 1
        return nomlz_par


class CST3DSurface(object):
    # generate surface by CST parameter
    # process order:
    # generate origin surface, deform, rotation(global y z x), translation

    def __init__(self, LX=0.0, LY=0.0, LZ=0.0, shape_par_X=None, shape_par_Y=None, shape_par_Z=None, class_par_Z=None, symmetry_y: bool = False):
        self.tran_fun_X = None
        self.tran_fun_Y = None
        self.tran_fun_Z = None
        self.rotation_matrix = None
        self.translation = None
        # notice do not support np.array input

        # generate 3D CST surface by LX,LY,LZ,shape_par_Y,shape_par_Z,class_par_Z

        # xi, x=X
        # psi, y=Y(xi)
        # zeta, z=Z(xi,psi)

        # input:
        # LX,LY,LZ,shape_par_Y,shape_par_Z,class_par_Z,symmetry_y

        # notice:
        # shape_fcn_X(PSI), shape_fcn_Y(XI), shape_fcn_Z(XI,PSI), class_fcn_Z{N1_fcn_Z(XI), N2_fcn_Z(XI)}
        # if input N1, N2 == 0, shape_fun will equal to 1
        # one of shape_fcn_X and shape_par_Y should be empty
        self.LX = LX
        self.LY = LY
        self.LZ = LZ
        self.symmetry_y = symmetry_y

        # preprocess
        if (type(shape_par_X) is np.ndarray) or (type(shape_par_X) is list) and (len(shape_par_X) == 0):
            shape_par_X = None
        if (type(shape_par_Y) is np.ndarray) or (type(shape_par_Y) is list) and (len(shape_par_Y) == 0):
            shape_par_Y = None
        if (type(shape_par_Z) is np.ndarray) or (type(shape_par_Z) is list) and (len(shape_par_Z) == 0):
            shape_par_Z = None
        if (type(class_par_Z) is np.ndarray) or (type(class_par_Z) is list) and (len(class_par_Z) == 0):
            class_par_Z = None

        if (shape_par_X is not None) and (shape_par_Y is not None):
            raise Exception(
                'CST3DSurface: donot support both X and Y deform')
        elif shape_par_X is not None:
            # deform direction y
            self.deform_ID = 1
        elif shape_par_Y is not None:
            # deform direction x
            self.deform_ID = 2
        else:
            self.deform_ID = 0

        if (shape_par_X is not None):
            if (type(shape_par_X) is np.ndarray) or (type(shape_par_X) is list):
                # shape_par_Y only one cross-section parameter
                self.shape_par_X = np.array(shape_par_X)
                self.shape_fcn_X = self.defunShapeX
            else:
                # function handle
                self.shape_fcn_X = shape_par_X
        else:
            self.shape_fcn_X = None

        if (shape_par_Y is not None):
            if (type(shape_par_Y) is np.ndarray) or (type(shape_par_Y) is list):
                # shape_par_Y only one cross-section parameter
                self.shape_par_Y = np.array(shape_par_Y)
                self.shape_fcn_Y = self.defunShapeY
            else:
                # function handle
                self.shape_fcn_Y = shape_par_Y
        else:
            self.shape_fcn_Y = None

        if (shape_par_Z is not None) or (class_par_Z is not None):
            if shape_par_Z is None:
                self.shape_par_Z = None
                self.shape_fcn_Z = None
            else:
                if (type(shape_par_Z) is np.ndarray) or (type(shape_par_Z) is list):
                    # shape_par_Z can have two parameter cross-section parameter
                    # default shape function donot
                    shape_par_Z = np.array(shape_par_Z)
                    if shape_par_Z.ndim == 1:
                        shape_par_Z = np.array([shape_par_Z, shape_par_Z])

                    if shape_par_Z.shape[0] == 1:
                        shape_par_Z = np.tile(shape_par_Z, [2, 1])

                    self.shape_par_Z = shape_par_Z
                    self.shape_fcn_Z = self.defunShapeZ
                else:
                    # function handle
                    self.shape_par_Z = None
                    self.shape_fcn_Z = shape_par_Z

            if (hasattr(class_par_Z, '__call__')):
                self.class_fcn_Z = class_par_Z
            else:
                # class_par_Z is {N1,N2}, N1/N2 hace two parameter(means two cross-section)
                if len(class_par_Z) != 2:
                    raise Exception(
                        'calCST: z class function parameter error input')

                if (type(class_par_Z[0]) is np.ndarray) or (type(class_par_Z[0]) is list):
                    self.N1_fcn_Z = lambda XI: class_par_Z[0][0] * \
                        XI + class_par_Z[0][1] * (1 - XI)
                else:
                    # function handle
                    self.N1_fcn_Z = lambda XI: class_par_Z[0]

                if (type(class_par_Z[1]) is np.ndarray) or (type(class_par_Z[1]) is list):
                    self.N2_fcn_Z = lambda XI: class_par_Z[1][0] * \
                        XI + class_par_Z[1][1] * (1 - XI)
                else:
                    # function handle
                    self.N2_fcn_Z = lambda XI: class_par_Z[1]

                # default class function
                self.class_fcn_Z = self.defunClassZ
        else:
            self.shape_fcn_Z = None
            self.N1_fcn_Z = None
            self.N2_fcn_Z = None
            self.class_fcn_Z = None

        return

    def addDeform(self, tran_fun_X=None, tran_fun_Y=None, tran_fun_Z=None):
        # base on local coordinate deform or translation surface

        # input:
        # tran_fun_X(PSI), tran_fun_Y(XI), tran_fun_Z(XI,PSI)
        # preprocess
        if (type(tran_fun_X) is np.ndarray) or (type(tran_fun_X) is list) and (len(tran_fun_X) == 0):
            tran_fun_X = None
        if (type(tran_fun_Y) is np.ndarray) or (type(tran_fun_Y) is list) and (len(tran_fun_Y) == 0):
            tran_fun_Y = None
        if (type(tran_fun_Z) is np.ndarray) or (type(tran_fun_Z) is list) and (len(tran_fun_Z) == 0):
            tran_fun_Z = None

        if self.deform_ID == 1:
            if tran_fun_Y:
                raise Exception(
                    'CST3DSurface: donot support both X and Y deform')

        elif self.deform_ID == 2:
            if tran_fun_X:
                raise Exception(
                    'CST3DSurface: donot support both X and Y deform')

        self.tran_fun_X = tran_fun_X
        self.tran_fun_Y = tran_fun_Y
        self.tran_fun_Z = tran_fun_Z

        return

    def addRotation(self, ang_x=0, ang_y=0, ang_z=0):
        # base on angle to rotation surface
        # rotation order:
        # y, z, x

        # input:
        # ang_x(deg), ang_y(deg), ang_z(deg)

        matrix = np.eye(3)
        # process rotation
        if ang_y != 0:
            cRY = np.cos(ang_y / 180 * np.pi)
            sRY = np.sin(ang_y / 180 * np.pi)
            matrix = np.matmul(np.array(
                [[cRY, 0, sRY], [0, 1, 0], [- sRY, 0, cRY]]), matrix)

        if ang_z != 0:
            cRZ = np.cos(ang_z / 180 * np.pi)
            sRZ = np.sin(ang_z / 180 * np.pi)
            matrix = np.matmul(
                np.array([[cRZ, - sRZ, 0], [sRZ, cRZ, 0], [0, 0, 1]]), matrix)

        if ang_x != 0:
            cRX = np.cos(ang_x / 180 * np.pi)
            sRX = np.sin(ang_x / 180 * np.pi)
            matrix = np.matmul(np.array(
                [[1, 0, 0], [0, cRX, - sRX], [0, sRX, cRX]]), matrix)

        self.rotation_matrix = matrix
        return

    def addTranslation(self, tran_x=None, tran_y=None, tran_z=None):
        # base on angle to rotation surface

        self.translation = np.array([tran_x, tran_y, tran_z])
        return

    def calSurface(self, XI=10, PSI=10):
        # generate 3D CST matrix by XI,PSI
        # X=CX*XI
        # Y=CY*shape_fcn_Y.*PSI, if symmetry_y Y=CY*shape_fcn_Y.*(PSI-0.5)
        # Z=CZ*shape_fcn_Z.*class_fcn_Z

        # default
        # xi_list=linspace(0,1,xi_gird_num+1)
        # psi_list=linspace(0,1,psi_gird_num+1), if symmetry_y psi_list=linspace(0.5,1,psi_gird_num+1);
        # [XI,PSI]=meshgrid(xi_list,psi_list), colume is LaWGS format line

        # xi, x=X
        # psi, y=Y(xi)
        # zeta, z=Z(xi,psi)

        # input:
        # XI,PSI
        # xi_grid_number,psi_grid_number

        # notice:
        # shape_fcn_Y(XI), shape_fcn_Z(XI,PSI), class_fcn_Z{N1_fcn_Z(XI), N2_fcn_Z(XI)}

        # output:
        # X,Y,Z (colume is LaWGS format line)

        # xi
        if isinstance(XI, int):
            xi_grid_numebr = XI
            XI = []
        else:
            # mean input XI matrix
            xi_grid_numebr = XI.shape[1] - 1

        # psi
        if isinstance(PSI, int):
            psi_grid_numebr = PSI
            PSI = []
        else:
            # mean input PSI matrix
            psi_grid_numebr = PSI.shape[0] - 1

        # calculate local coordinate matrix
        if len(XI) == 0:
            XI = np.tile([np.linspace(
                0, 1, xi_grid_numebr + 1)], [psi_grid_numebr + 1, 1])

        if len(PSI) == 0:
            if self.symmetry_y:
                PSI = np.tile(np.transpose([np.linspace(
                    0.5, 1, psi_grid_numebr + 1)]), [1, xi_grid_numebr + 1])
            else:
                PSI = np.tile(np.transpose([np.linspace(
                    0, 1, psi_grid_numebr + 1)]), [1, xi_grid_numebr + 1])

        X, Y, Z = self.calPoint(XI, PSI)
        return X, Y, Z

    def calPoint(self, XI, PSI):
        # calculate point on surface

        # calculate origin surface matrix
        X = self.LX * XI
        if self.shape_fcn_X:
            X = X * self.shape_fcn_X(PSI)

        if self.symmetry_y:
            Y = self.LY * (PSI - 0.5)
        else:
            Y = self.LY * PSI

        if self.shape_fcn_Y:
            Y = Y * self.shape_fcn_Y(XI)

        if (not self.shape_fcn_Z) and (not self.class_fcn_Z):
            if len(X.shape) == 1:
                Z = np.zeros((len(X)))
            else:
                Z = np.zeros((X.shape[0], X.shape[1]))
        elif (not self.shape_fcn_Z):
            Z = self.LZ * self.class_fcn_Z(XI, PSI)
        elif (not self.class_fcn_Z):
            Z = self.LZ * self.shape_fcn_Z(XI, PSI)
        else:
            Z = self.LZ * \
                self.shape_fcn_Z(XI, PSI) * self.class_fcn_Z(XI, PSI)

        # deform surface
        if self.tran_fun_X:
            X = X + self.tran_fun_X(PSI)

        if self.tran_fun_Y:
            Y = Y + self.tran_fun_Y(XI)

        if self.tran_fun_Z:
            Z = Z + self.tran_fun_Z(XI, PSI)

        # rotation surface
        if self.rotation_matrix is not None:
            matrix = self.rotation_matrix
            X_old = X
            Y_old = Y
            Z_old = Z
            X = matrix[0][0] * X_old + matrix[0][1] * \
                Y_old + matrix[0][2] * Z_old
            Y = matrix[1][0] * X_old + matrix[1][1] * \
                Y_old + matrix[1][2] * Z_old
            Z = matrix[2][0] * X_old + matrix[2][1] * \
                Y_old + matrix[2][2] * Z_old

        # translation surface
        if self.translation is not None:
            X = X + self.translation[0]
            Y = Y + self.translation[1]
            Z = Z + self.translation[2]

        return X, Y, Z

    def calCoordinate(self, X, Y, Z):
        # base on X, Y calculate local coordinate in surface
        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)

        # undone translation surface
        if not (self.translation is None):
            X = X - self.translation[0]
            Y = Y - self.translation[1]
            Z = Z - self.translation[2]

        # undone rotation surface
        if not (self.rotation_matrix is None):
            matrix = np.transpose(self.rotation_matrix)
            X_old = X
            Y_old = Y
            Z_old = Z
            X = matrix[0][0] * X_old + matrix[0][1] * \
                Y_old + matrix[0][2] * Z_old
            Y = matrix[1][0] * X_old + matrix[1][1] * \
                Y_old + matrix[1][2] * Z_old

        # identify base local coordinate
        if self.deform_ID == 1:
            PSI = Y / self.LY
            if self.symmetry_y:
                PSI = PSI + 0.5
            # re deform surface
            if not (self.tran_fun_X is None):
                X = X - self.tran_fun_X(PSI)
            XI = X / self.LX
            if not (self.shape_fcn_X is None):
                XI = XI / self.shape_fcn_X(PSI)
        else:
            if self.deform_ID == 2:
                XI = X / self.LX
                # re deform surface
                if not (self.tran_fun_Y is None):
                    Y = Y - self.tran_fun_Y(XI)
                PSI = Y / self.LY
                if not (self.shape_fcn_Y is None):
                    PSI = PSI / self.shape_fcn_Y(XI)
                if self.symmetry_y:
                    PSI = PSI + 0.5

        return XI, PSI

    def defunShapeX(self, PSI):
        # default shape function of X

        NP = self.calNormPar(self.shape_par_X[0], self.shape_par_X[1])
        S = (PSI ** self.shape_par_X[0] *
             (1 - PSI) ** self.shape_par_X[0]) / NP
        return S

    def defunShapeY(self, XI):
        # default shape function of Y

        NP = self.calNormPar(self.shape_par_Y[0], self.shape_par_Y[1])
        S = (XI ** self.shape_par_Y[0] * (1 - XI) ** self.shape_par_Y[1]) / NP
        return S

    def defunShapeZ(self, XI, PSI):
        # default shape function of Z
        # default N1, N2 is relation with PSI

        N1 = self.shape_par_Z[0][0] * (1 - PSI) + self.shape_par_Z[1][0] * PSI
        N2 = self.shape_par_Z[0][1] * (1 - PSI) + self.shape_par_Z[1][1] * PSI
        NP = self.calNormPar(N1, N2)
        S = XI ** N1 * (1 - XI) ** N2 / NP
        return S

    def defunClassZ(self, XI, PSI):
        # default class function of Z
        # default N1, N2 is relation with XI

        N1 = self.N1_fcn_Z(XI)
        N2 = self.N2_fcn_Z(XI)
        NP = self.calNormPar(N1, N2)
        PSI[PSI < 0] = 0
        PSI[PSI > 1] = 1
        C = (PSI ** N1 * (1 - PSI) ** N2) / NP
        return C

    def calNormPar(self, N1, N2):
        # calculate normailize class function parameter by N1, N2
        N1 = np.array(N1)
        N2 = np.array(N2)
        N1[np.logical_and(N1 == 0, N2 == 0)] = 1
        N2[np.logical_and(N1 == 0, N2 == 0)] = 1
        nomlz_par = np.array((N1 / (N1 + N2)) ** N1 * (N2 / (N1 + N2)) ** N2)
        nomlz_par[np.logical_and(N1 == 0, N2 == 0)] = 1
        return nomlz_par


def write_2D(
    vecPsi: np.ndarray,
    vecZetaUp: np.ndarray,
    vecZetaLow: np.ndarray,
    varX: float = 0.,
    varLen: float = 1.,
    strFileName: str = None,
):
    '''
    Create .csv file using airfoil data points.
    '''
    numPnt = len(vecPsi)
    vecX = vecPsi * varLen + varX
    vecZUp = vecZetaUp * varLen
    vecZLow = vecZetaLow * varLen

    if strFileName is None:
        strFileName = "airfoil.csv"
    hdlFile = open(strFileName, "w")
    hdlWriter = csv.writer(hdlFile)
    vecHeader = ["xup", "zup", "xlow", "zlow"]
    hdlWriter.writerow(vecHeader)
    for idx in range(numPnt):
        strContent = ["{}".format(vecX[idx]), "{}".format(
            vecZUp[idx]), "{}".format(vecX[idx]), "{}".format(vecZLow[idx])]
        hdlWriter.writerow(strContent)
    hdlFile.close()
    return 1


def plot_3D(
    matX: Union[List[float], np.ndarray],
    matY: Union[List[float], np.ndarray],
    matZ1: Union[List[float], np.ndarray],
    matZ2: Union[List[float], np.ndarray] = None,
):
    '''
    Plot 3D plot.

    Inputs
    ------
    matX, matY, matZ : matrices of X, Y, and Z coordinates

    Output
    ------
    plot
    '''
    hdlFig = plt.figure()
    hdlAxe = Axes3D(hdlFig, auto_add_to_figure=False)
    hdlFig.add_axes(hdlAxe)
    # hdlAxe.plot_surface(matX, matY, matZ1, rstride = 1, cstride = 1, cmap='Greys', linewidth=0, antialiased=False)
    # hdlAxe.plot_surface(matX, matY, matZ2, rstride = 1, cstride = 1, cmap='Greys', linewidth=0, antialiased=False)
    hdlAxe.plot_surface(matX, matY, matZ1, rstride=1, cstride=1, cmap='Greys')
    if not matZ2 is None:
        hdlAxe.plot_surface(matX, matY, matZ2, rstride=1,
                            cstride=1, cmap='Greys')
    hdlAxe.set_xlabel(r'$x$', fontsize=18)
    hdlAxe.set_ylabel(r'$y$', fontsize=18)
    hdlAxe.set_zlabel(r'$z$', fontsize=18)
    # hdlAxe.set_axis_off()
    # hdlAxe.view_init(30,-120)
    hdlAxe.view_init(0, 0)
    _set_axes_equal(hdlAxe)
    plt.show()
    return 1


def plot_3D_multiple(
    listMat: List[np.ndarray]
):
    '''
    Draw multiple 3D plots defined by list.

    Inputs
    ------
    listMat : list of matrices, every three items define one part. in the sequence of matX, matY, matZ.

    Output
    ------
    plot.
    '''
    if np.mod(len(listMat), 3) != 0:
        raise RuntimeError('data matrices list define error')
    else:
        numPlt = len(listMat)/3

    hdlFig = plt.figure()
    hdlAxe = Axes3D(hdlFig)
    for idx in range(numPlt):
        matX, matY, matZ = listMat[3*idx], listMat[3*idx+1], listMat[3*idx+2]
        hdlAxe.plot_surface(matX, matY, matZ, rstride=1,
                            cstride=1, cmap='Greys')
    hdlAxe.set_xlabel(r'$x$', fontsize=18)
    hdlAxe.set_ylabel(r'$y$', fontsize=18)
    hdlAxe.set_zlabel(r'$z$', fontsize=18)
    hdlAxe.set_axis_off()
    # hdlAxe.view_init(30,-120)
    hdlAxe.view_init(0, 0)
    _set_axes_equal(hdlAxe)
    plt.show()

    return 1


def _set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])
    return 1


def importdata(str_file):
    with open(str_file) as file_data:
        str_data = file_data.readlines()
    data = [data_unit.split() for data_unit in str_data]
    return data


if __name__ == "__main__":
    # def shape_par_Y(XI): return 1-XI*(1-0.4)
    # shape_par_Z = [0, 0]
    # class_par_Z = [[1, 1], [0.001, 0.001]]
    # test_surf = CST3DSurface(2, 1, 0.3, shape_par_Y,
    #                          [], shape_par_Z, class_par_Z)
    # test_surf.addDeform(shape_par_Y, [], [])
    # test_surf.addRotation(90, 0, 0)
    # test_surf.addTranslation(5, 6, 9)
    # [X, Y, Z] = test_surf.calSurface()
    # plot_3D(X, Y, Z)

    curve = CST2DCurve(1,1,[1,0.5],[0.5,1])
    print(curve.calCurve())