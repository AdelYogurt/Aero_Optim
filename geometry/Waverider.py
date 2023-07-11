import os
import sys
import cmath
import numpy as np
from numpy import fliplr, sqrt
from pyPanair.preprocess import wgs_creator
from scipy import io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from CST import CST3DSurface, plot_3D, plot_3D_multiple

# def genWaveriderWingFei(
#     varCoefW, varClsT, varCoefHu, varCoefHl, varClsNu, varClsNl, varClsM1u, varClsM1l, # coefficients in CST
#     varRtoHead, # ratio of head length
#     varRtoChrd, # ratio of chord length of 2nd section respect to root chord length
#     varSSpan1, varSSpan2, # semi-span of 1st and 2nd sections
#     varArflThck=0.05, # thickness of airfoils in 2nd section
#     varLen=4., # total length
#     varCoefG = 0., varClsF = 0., # coefficients in CST controling rare-end
#     varRad = 0.005, # radius of leading to alleviate heating
#     varArflThckMin = 0.02, # minimum airfoil thickness
#     varSweep1 = 0., # sweep angle variation of wing1
#     varSweep2Norm = None, # normalized sweep angle variation of wing2
#     varExp2Rate = 1.0, # expansion rate of wing2
#     flgHalf = True, # output half model or not
#     strName  = "waverider"
#     ):
#     """
#     Generate parameterized geometry model of waverider with wing.

#     Parameters boundaries
#     ---------------------
#     varCoefW[0.8, 2.0], varClsT[0.5, 1.5], varCoefHu[0.1, 0.5], varCoefHl[0.1,0.5], varClsNu[0.5, 2.0], varClsNl[0.5, 2.0], varClsM1u[0.1, 0.99], varClsM1l[0.1, 0.99], varRtoHead[0.55, 0.85], varRtoChrd[0.4, 1.0], varArflThck[0.02, 0.08], varSSpan1[0.6, 0.85], varSSpan2[0.4, 0.8]


#     Packages used:
#     * https://github.com/yufei96/pyCST
#     * https://github.com/SaTa999/pyPanair
#     """
#     wgs = wgs_creator.LaWGS(strName)

#     # Create waverider head
#     matX, matY, matZUp, matZLow = CST.fitting_waverider(varLen, varCoefW, varClsT, varCoefHu, varCoefHl, varClsNu, varClsNl, varClsM1u, varClsM1l, varCoefG, varClsF, varRad)
#     # matZUp += 5E-4  matZLow -= 5E-4
#     matZUp += varRad  matZLow -= varRad

#     # cutting head to varRtoHead
#     numGrdY, numGrdX = matX.shape
#     numPntY = numGrdY*2
#     numLineHead = int((numGrdX-1) * varRtoHead)
#     listSection = list()
#     listSection.append(zeros((numPntY,3)))
#     for idx1 in range(1, numLineHead+1):
#         line = ones((numPntY, 3))
#         line[:,0] *= matX[0,idx1]
#         line[:,1]  = hstack((matY[:,idx1], matY[::-1,idx1]))
#         line[:,2]  = hstack((matZUp[:,idx1], matZLow[::-1,idx1]))
#         if flgHalf:
#             line = delete(line, where(line[:,1]<0)[0], axis=0)
#         listSection.append(line)
#     if flgHalf:
#         listSection[0] = zeros((numGrdY+1,3))
#     netHead = wgs_creator.Network(listSection)
#     wgs.append_network("head", netHead, 1)
#     # wgs.create_stl(filename="head.stl")
#     # netHead.plot_wireframe()


#     # Create aft-body (longitudinal-ward), upper and lower sides
#     vecY = matY[:, numLineHead]
#     matAftX    = tile(matX[0,numLineHead:], (numGrdY,1))
#     numGrdXAft = matAftX.shape[1]
#     matAftY    = zeros((numGrdY, numGrdXAft))
#     matAftZUp  = zeros((numGrdY, numGrdXAft))
#     matAftZLow = zeros((numGrdY, numGrdXAft))
#     matAftY[:,0] = vecY
#     matAftZUp[:,0], matAftZLow[:,0] = matZUp[:,numLineHead], matZLow[:,numLineHead]
#     varYLim = max(matY[:,numLineHead]) # maximum width defined by head
#     for idx2 in range(numLineHead+1,numGrdX):
#         _vecY1, _vecZUp1, _vecZLow1 = matY[:,idx2], matZUp[:,idx2], matZLow[:,idx2]
#         _vecY2 = linspace(-varYLim, varYLim, num=numGrdY)
#         _ratio = (numGrdX-idx2-1) / (numGrdXAft-1)
#         _vecY3 = _ratio*vecY + (1-_ratio)*_vecY2 # avoid slope changing rapidly
#         # _vecY3 = vecY
#         matAftY[:,idx2-numLineHead] = _vecY3
#         fcnIntpUp  = InterpolatedUnivariateSpline(_vecY1, _vecZUp1)
#         fcnIntpLow = InterpolatedUnivariateSpline(_vecY1, _vecZLow1)
#         matAftZUp[:,idx2-numLineHead]  = fcnIntpUp(_vecY3) # coresponding z coordinates
#         matAftZLow[:,idx2-numLineHead] = fcnIntpLow(_vecY3)
#     matAftXUp, matAftYUp, matAftZUpNew = fliplr(matAftX), fliplr(matAftY), fliplr(matAftZUp)

#     listSectionUp  = list()
#     listSectionLow = list()
#     for idx3 in range(numGrdY):
#         lineUp, lineLow = zeros((numGrdXAft, 3)), ones((numGrdXAft, 3))
#         lineUp[:,0], lineLow[:,0] = matAftXUp[idx3,:], matAftX[idx3,:]
#         lineUp[:,1], lineLow[:,1] = matAftYUp[idx3,:], matAftY[idx3,:]
#         lineUp[:,2], lineLow[:,2] = matAftZUpNew[idx3,:], matAftZLow[idx3,:]
#         if flgHalf:
#             if lineUp[0,1]<-1E-9:
#                 continue
#         listSectionUp.append(lineUp)
#         listSectionLow.append(lineLow)
#     netAftUp  = wgs_creator.Network(listSectionUp)
#     wgs.append_network("aft_up", netAftUp, 1)
#     # netAftUp.plot_wireframe()
#     netAftLow = wgs_creator.Network(listSectionLow)
#     wgs.append_network("aft_low", netAftLow, 1)
#     # wgsAft = wgs_creator.LaWGS("aft")
#     # wgsAft.append_network("aft_up", netAftUp, 1)
#     # wgsAft.append_network("aft_low", netAftLow, 1)
#     # wgsAft.create_stl()
#     # netAftLow.plot_wireframe()


#     # Create wings
#     # actually the root chord length does not equal to the aft-body length, but they are assumed to be identical for geometry smoothness
#     try:
#         idxLeading = where(matAftZUp[-1,:]-matAftZLow[-1,:] >= varArflThckMin)[0][0]
#     except IndexError:
#         idxLeading = 0
#     if idxLeading == matAftZUp.shape[1]-1:
#         idxLeading = 0
#     # idxLeading = 0
#     # idxLeading = where(matAftZUp[-1,:]-matAftZLow[-1,:] >= varArflThckMin)[0][0]
#     numGrdWing = numGrdXAft
#     numGrdEdge = int(round((numGrdWing-1)/2))
#     numGrid0Trail = 5
#     lineSctn0  = ones((2*numGrdWing+2+numGrid0Trail, 3))
#     lineSctn1  = ones((2*numGrdWing+2+numGrid0Trail, 3))
#     lineSctn2  = ones((2*numGrdWing+2+numGrid0Trail, 3))
#     x0, y0, z0 = matAftX[0,idxLeading], matAftY[-1,idxLeading], (matAftZUp[-1,idxLeading]+matAftZLow[-1,idxLeading])/2 #*基准外形优化时用了这个
#     # x0, y0, z0 = matAftX[0,0], matAftY[-1,0], (matAftZUp[-1,0]+matAftZLow[-1,0])/2 #*写论文介绍参数时应该用这个，几何意义更明确
#     varLenChrd1 = varLen - x0
#     z2 = [i*(matAftZUp[-1,-1]-matAftZLow[-1,-1])/(numGrid0Trail+1)+matAftZLow[-1,-1] for i in range(numGrid0Trail,0,-1)]
#     x0trail, y0trail = varLen, y0
#     lineSctn0[:,0] = hstack((matAftX[-1,0], matAftX[-1,:], matAftX[-1,-1]*ones(numGrid0Trail), matAftX[-1,::-1], matAftX[-1,0]))
#     lineSctn0[:,1] *= matAftY[-1,0]
#     lineSctn0[:,2] = hstack((matAftZUp[-1,0], matAftZUp[-1,:], z2, matAftZLow[-1,::-1], matAftZUp[-1,0]))
#     varLenChrd2 = varRtoChrd * varLenChrd1
#     x1, y1, z1 = x0+(varLenChrd1-varLenChrd2), y0+varSSpan1, z0
#     varAngl1 = atan2(x1-x0, y1-y0)
#     varLenLead1 = sqrt((x1-x0)**2 + (y1-y0)**2)
#     varAngl1New = varAngl1 + varSweep1
#     x1new, y1new, z1new = x0 + varLenLead1 * sin(varAngl1New), y0 + varLenLead1 * cos(varAngl1New), z1 # wing1 sweep morph
#     x1trailnew, y1trailnew = x1new + varLenChrd2, y1new
#     varAngl1TrailNew = atan2(x1trailnew-x0trail, y1trailnew-y0trail)
#     lineSctn1[:,0] = hstack((x1new, linspace(x1new,x1new+varLenChrd2,numGrdWing), (x1new+varLenChrd2)*ones(numGrid0Trail), linspace(x1new+varLenChrd2,x1new,numGrdWing), x1new))
#     lineSctn1[:,1] *= y1new
#     lineSctn1[:,2] = hstack((z1new, linspace(z1new+varArflThckMin/2,z1new+varArflThck-varArflThckMin/2,numGrdEdge+1), linspace(z1new+varArflThck-varArflThckMin/2,z1new+varArflThckMin/2,numGrdWing-numGrdEdge)[1:], z1new*ones(numGrid0Trail), (z1new-varArflThckMin/2)*ones((1,numGrdWing))[0], z1new))
#     varLenLead2 = varSSpan2 * varExp2Rate
#     if varSweep2Norm is None:
#         varAngl2New = 0.0
#     else:
#         varAngl2New = varAngl1TrailNew + (varAngl1New - varAngl1TrailNew) * varSweep2Norm
#     lineSctn2[:,0] = lineSctn1[:,0] + varLenLead2 * sin(varAngl2New)
#     lineSctn2[:,1] = lineSctn1[:,1] + varLenLead2 * cos(varAngl2New)
#     lineSctn2[:,2] = lineSctn1[:,2]
#     lenGrid = (matAftX[0,-1]-matAftX[0,0])/numGrdXAft
#     netWing1 = wgs_creator.Line(lineSctn0).linspace(lineSctn1, num=int(varSSpan1/lenGrid))
#     netWing2 = wgs_creator.Line(lineSctn1).linspace(lineSctn2, num=int(varSSpan2/lenGrid))
#     wgs.append_network('wingin_right', flipud(netWing1), 1)
#     wgs.append_network('wingout_right', flipud(netWing2), 1)
#     # wgsWingIn = wgs_creator.LaWGS('wingin')
#     # wgsWingIn.append_network('wingin', netWing1, 1)
#     # wgsWingIn.create_stl('BWB-wingin.stl')
#     lineTipRghtUp  = lineSctn2[1:numGrdWing+1,:]
#     lineTipRghtLow = lineSctn2[-2:-numGrdWing-2:-1,:]
#     netWingTipRght = wgs_creator.Line(lineTipRghtUp).linspace(lineTipRghtLow, 3)
#     wgs.append_network("wing_tip_right", flipud(netWingTipRght), 1)
#     # wgsWingOut = wgs_creator.LaWGS('wingout')
#     # wgsWingOut.append_network('wingout', netWing2, 1)
#     # wgsWingOut.append_network('tip', netWingTipRght, 1)
#     # wgsWingOut.create_stl('BWB-wingout.stl')
#     if not flgHalf:
#         netWing3 = netWing1.copy()  netWing4 = netWing2.copy()
#         netWing3[:,:,1] = -netWing1[:,:,1]  netWing4[:,:,1] = -netWing2[:,:,1]
#         wgs.append_network("wingin_left", netWing3, 1)
#         wgs.append_network('wingout_left', netWing4, 1)
#         netWingTipLeft = netWingTipRght.copy()
#         netWingTipLeft[:,:,1] = -netWingTipLeft[:,:,1]
#         wgs.append_network("wing_tip_left", netWingTipLeft, 1)


#     # Create body base
#     lineBase1 = ones((2*numGrdY+1+2*numGrid0Trail,3))
#     lineBase1[:,0] *= matAftX[0,-1]
#     lineBase1[:,1] = hstack((matAftY[0,-1]*ones(numGrid0Trail), matAftY[:,-1], matAftY[-1,-1]*ones(numGrid0Trail), matAftY[::-1,-1], matAftY[0,-1]))
#     lineBase1[:,2] = hstack((z2[::-1], matAftZUp[:,-1], z2, matAftZLow[:,-1], z2[-1]))
#     lineBase1 = wgs_creator.Line(lineBase1)
#     if flgHalf:
#         lineBase1 = delete(lineBase1, where(lineBase1[:,1]<-1E-9)[0], axis=0)
#     lineBase2 = lineBase1.replace(x=matAftX[0,-1], y=matAftY[int((numGrdY-1)/2),-1], z=(matAftZUp[int((numGrdY-1)/2),-1]+matAftZLow[int((numGrdY-1)/2),-1])/2)
#     netBodyBase = lineBase1.linspace(lineBase2, num=7)
#     # netBodyBase.plot_wireframe()
#     wgs.append_network("body_base", netBodyBase, 5)
#     # wgsBody = wgs_creator.LaWGS('body')
#     # wgsBody.append_network('head', netHead, 1)
#     # wgsBody.append_network('aftup', netAftUp, 1)
#     # wgsBody.append_network('aftlow', netAftLow, 1)
#     # wgsBody.append_network('base', netBodyBase, 1)
#     # wgsBody.create_stl('BWB-body.stl')

#     # wgs.create_stl(include_wake=True)
#     return wgs


# def genWaverider(
#     varCoefW, varClsT, varCoefHu, varCoefHl, varClsNu, varClsNl, varClsMu, varClsMl, # coefficients affecting surfaces shape
#     varLen=4., # total length
#     varCoefG=0., varClsF=0., # coefficients in CST controling rare-end
#     varRad=0.005, # radius of leading to alleviate heating
#     flgHalf=True, # output half model or not
#     strName='waverider'
#     ):
#     '''
#     Generate parameterized geometry model for a pure waverider, without wings.

#     Parameters boundaries
#     ---------------------
#     varCoefHl[0.1,0.5], varClsNl[0.5, 2.0], varClsMl[0.1, 0.99]
#     '''
#     wgs = wgs_creator.LaWGS(strName)

#     matX, matY, matZUp, matZLow = CST.fitting_waverider(varLen, varCoefW, varClsT, varCoefHu, varCoefHl, varClsNu, varClsNl, varClsMu, varClsMl, varCoefG, varClsF, varRad)
#     matZUp += 1E-3  matZLow -= 1E-3

#     numGrdY, numGrdX = matX.shape
#     listSection = list()
#     listSection.append(zeros((numGrdY*2,3)))
#     for idx in range(1, numGrdX):
#         line = ones((numGrdY*2, 3))
#         line[:,0] *= matX[0,idx]
#         line[:,1] = hstack((matY[:,idx], matY[::-1,idx]))
#         line[:,2] = hstack((matZUp[:,idx], matZLow[::-1,idx]))
#         if flgHalf:
#             line = delete(line, where(line[:,1]<0)[0], axis=0)
#         listSection.append(line)
#     if flgHalf:
#         listSection[0] = zeros((numGrdY+1, 3))
#     netBody = wgs_creator.Network(listSection)
#     wgs.append_network('body', netBody, 1)

#     lineBase1 = wgs_creator.Line(listSection[-1])
#     lineBase2 = lineBase1.replace(x=matX[0,-1], y=matY[int((numGrdY-1)/2),-1], z=(matZUp[int((numGrdY-1)/2),-1]+matZLow[int((numGrdY-1)/2),-1])/2)
#     netBodyBase = lineBase1.linspace(lineBase2, num=9)
#     wgs.append_network('body_base', netBodyBase, 5)

#     wgs.create_stl(filename='waverider.stl')
#     wgs.create_wgs(filename='waverider.wgs')

#     return wgs


class WaveriderWingDia(object):
    def __init__(self, total_length=None, par_width=None, par_hight_up=None, par_hight_low=None,
                 par_T=None, par_M_up=None, par_M_low=None, par_N_up=None, par_N_low=None, par_R=None,
                 par_rho1=None, par_rho12=None, par_WS1=None, par_WS2=None, par_WT_up=None, par_WT_low=None, half_flag: bool = False):
        # generate waverider wing wgs data
        # xi, x=X(xi)
        # psi, y=Y(xi)
        # zeta, z=Z(xi,psi)

        # notice:
        # colume vector in matrix is wgs format
        # par_rho1 should between 0 and 1

        if par_rho1 >= 1 or par_rho1 <= 0:
            raise Exception('genWaveriderWing: do not support pure waverider')

        global_symmetry_y = True
        self.total_length = total_length
        self.par_width = par_width
        self.par_hight_up = par_hight_up
        self.par_hight_low = par_hight_low
        self.par_T = par_T
        self.par_M_up = par_M_up
        self.par_M_low = par_M_low
        self.par_N_up = par_N_up
        self.par_N_low = par_N_low
        self.par_R = par_R
        self.par_rho1 = par_rho1
        self.par_rho12 = par_rho12
        self.par_WS1 = par_WS1
        self.par_WS2 = par_WS2
        self.par_WT_up = par_WT_up
        self.par_WT_low = par_WT_low
        self.half_flag = half_flag

        self.first_xi = 0.002  # is local length parameter
        self.stag_xi = 0.01  # is local length parameter

        stag_length = (self.stag_xi*(1-par_rho1))**par_T*par_width/2
        self.stag_length = stag_length

        head_length = (1-par_rho1)*total_length
        self.head_length = head_length

        head_side_length = head_length*(1-self.stag_xi)
        self.head_side_length = head_side_length

        body_length = par_rho1 * total_length
        self.body_length = body_length

        y_cut = (1 - par_rho1) ** par_T * par_width / 2
        self.y_cut = y_cut
        def body_psi_max_fun(XI): return y_cut/(XI**par_T*par_width)-2.2204e-16
        self.body_psi_max_fun = body_psi_max_fun

        # wing base parameter
        psi_end_cut = self.body_psi_max_fun(1)
        tri_wing_length = par_rho1 * total_length
        tri_wing_width = par_WS1
        tri_wing_height_up = (psi_end_cut + 0.5) ** par_N_up * (0.5 -
                                                                psi_end_cut) ** par_N_up / (0.5) ** (2 * par_N_up) * par_hight_up
        tri_wing_height_low = (psi_end_cut + 0.5) ** par_N_low * (
            0.5 - psi_end_cut) ** par_N_low / (0.5) ** (2 * par_N_low) * par_hight_low
        wing_length = par_rho1 * par_rho12 * total_length
        wing_width = par_WS2
        wing_height_up = par_WT_up
        wing_height_low = par_WT_low

        # calculate blunt radius parameter
        # local coordinate, local x is global x, local y is global z
        # head radius
        # calculate gradient of up and low discrete surface to local radius center
        if par_M_up > 1:
            radius_center_up = 0
        else:
            dz_dx = (par_hight_up * (self.first_xi) ** par_M_up) / \
                (total_length * self.first_xi)
            radius_center_up = par_R * dz_dx

        if par_M_low > 1:
            radius_center_low = 0
        else:
            dz_dx = (par_hight_low * (self.first_xi) ** par_M_low) / \
                (total_length * self.first_xi)
            radius_center_low = par_R * dz_dx

        radius_center_head = max(radius_center_up , radius_center_low)
        radius_head_sq = (radius_center_head *
                          radius_center_head + par_R * par_R)
        radius_head = np.sqrt(radius_head_sq)

        def shape_head_edge_x(Y): return np.sqrt(
            radius_head_sq-Y**2)-radius_center_head

        # calculate blunt radius of zox plane of tri wing
        dz_dx = (par_M_up*(1-par_rho1)**(par_M_up-1)/total_length)
        radius_center_up = par_R*dz_dx
        dz_dx = (par_M_low*(1-par_rho1)**(par_M_low-1)/total_length)
        radius_center_low = par_R*dz_dx
        radius_center_tri_wing = max(radius_center_up, radius_center_low)
        radius_tri_wing_sq = (radius_center_tri_wing *
                              radius_center_tri_wing+par_R*par_R)

        def shape_tri_wing_edge_x(Y): return sqrt(
            radius_tri_wing_sq-Y**2)-radius_center_tri_wing

        # calculate blunt radius and center of zox plane of wing
        dz_dx = par_WT_up / wing_length * 2
        radius_center_up = par_R * dz_dx
        dz_dx = par_WT_low / wing_length * 2
        radius_center_low = par_R * dz_dx
        radius_center_wing = max(radius_center_up , radius_center_low)
        radius_wing_sq = (radius_center_wing *
                          radius_center_wing + par_R * par_R)
        radius_wing = np.sqrt(radius_wing_sq)
        def shape_wing_edge_x(Y): return np.sqrt(
            radius_wing_sq-Y**2)-radius_center_wing

        # calculate head side blunt radius of yoz plane
        dz_dx = (0.005)**par_N_up*(0.995)**par_N_up/(0.5)**(2*par_N_up) * \
            par_hight_up*(1-par_rho1)**par_M_up/(0.01*y_cut)
        radius_center_up = par_R*dz_dx
        dz_dx = (0.005)**par_N_low*(0.995)**par_N_low/(0.5)**(2*par_N_low) * \
            par_hight_low*(1-par_rho1)**par_M_low/(0.01*y_cut)
        radius_center_low = par_R*dz_dx
        radius_center_head_side = max(radius_center_low, radius_center_up)
        radius_head_side_sq = (radius_center_head_side *
                               radius_center_head_side+par_R*par_R)
        radius_head_side = sqrt(radius_head_side_sq)

        # define head
        def shape_fcn_Y(XI): return (XI*(1-par_rho1))**par_T
        def shape_fcn_Z_up(XI, PSI): return (XI*(1-par_rho1))**par_M_up
        def shape_fcn_Z_low(XI, PSI): return (XI*(1-par_rho1))**par_M_low
        head_up = CST3DSurface(head_length, par_width, par_hight_up, [],
                               shape_fcn_Y, shape_fcn_Z_up, [par_N_up, par_N_up], global_symmetry_y)
        head_low = CST3DSurface(head_length, par_width, - par_hight_low, [],
                                shape_fcn_Y, shape_fcn_Z_low, [par_N_low, par_N_low], global_symmetry_y)

        # blunt translation
        head_up.addTranslation(0, 0, par_R)
        head_low.addTranslation(0, 0, - par_R)

        # # define stag
        # # local coordinate, local x is global y, local y is global z, local z is global x
        # def shape_fcn_X(PSI): return (sqrt(radius_head_side_sq-(((PSI-0.5)*2)
        #                                                         * par_R)**2)-radius_center_head_side)/stag_length+1

        # def shape_fcn_Z(XI, PSI): return (1-XI) * \
        #     shape_head_edge_x((PSI-0.5)*par_R*2)
        # stag = CST3DSurface(stag_length, 2*par_R, -1,
        #                     shape_fcn_X, [], [0.0, 0.0], shape_fcn_Z)

        # def X_edge(PSI): return total_length*(PSI)**(1/par_T)
        # def tran_fun_Z(XI, PSI): return X_edge(XI*stag_length/(par_width/2))
        # # deform surface
        # stag.addDeform([], [], tran_fun_Z)
        # # rotation surface
        # stag.addRotation(90, 90, 0)
        # # translation surface
        # stag.addTranslation(0, 0, -par_R)

        # define blunt head side
        # local coordinate, local x is global -x, local y is global z, local z is global y
        def shape_fcn_X(PSI): return 1 - \
            shape_tri_wing_edge_x((PSI-0.5)*par_R*2)/head_length+shape_head_edge_x((PSI-0.5)*par_R*2)/head_length

        def shape_fcn_Z(XI, PSI): return sqrt(
            radius_head_side_sq-(((PSI-0.5)*2)*par_R)**2)-radius_center_head_side
        head_side = CST3DSurface(
            head_length, 2 * par_R, 1, shape_fcn_X, [], [0.0, 0.0], shape_fcn_Z)

        def tran_fun_X(PSI): return shape_tri_wing_edge_x((PSI-0.5)*par_R*2)
        def Y_edge(XI): return XI**par_T*par_width/2
        def tran_fun_Z(XI, PSI): return Y_edge((1-XI)*(1-par_rho1))

        # deform surface
        head_side.addDeform(tran_fun_X, [], tran_fun_Z)
        # rotation surface
        head_side.addRotation(90, 180, 0)
        # translation surface
        head_side.addTranslation(head_length, 0, - par_R)

        # define body
        def shape_fcn_Z_up(XI, PSI): return (
            XI*par_rho1+(1-par_rho1))**par_M_up

        def class_fcn_Z_up(XI, PSI): return ((PSI-0.5)*2*body_psi_max_fun(XI*par_rho1+(1-par_rho1))+0.5)**par_N_up*(
            1-((PSI-0.5)*2*body_psi_max_fun(XI*par_rho1+(1-par_rho1))+0.5))**par_N_up/(0.5)**(2*par_N_up)

        def shape_fcn_Z_low(XI, PSI): return (
            XI*par_rho1+(1-par_rho1))**par_M_low

        def class_fcn_Z_low(XI, PSI): return ((PSI-0.5)*2*body_psi_max_fun(XI*par_rho1+(1-par_rho1))+0.5)**par_N_low*(
            1-((PSI-0.5)*2*body_psi_max_fun(XI*par_rho1+(1-par_rho1))+0.5))**par_N_low/(0.5)**(2*par_N_low)
        body_up = CST3DSurface(body_length, y_cut * 2, par_hight_up,
                               [], [], shape_fcn_Z_up, class_fcn_Z_up, global_symmetry_y)
        body_low = CST3DSurface(body_length, y_cut * 2, - par_hight_low,
                                [], [], shape_fcn_Z_low, class_fcn_Z_low, global_symmetry_y)

        # blunt translation
        body_up.addTranslation(head_length, 0, par_R)
        body_low.addTranslation(head_length, 0, - par_R)

        # define body back
        def shape_fcn_Y_up(XI): return par_hight_up*(0.5+XI*(psi_end_cut))**par_N_up * \
            (0.5-XI*(psi_end_cut))**par_N_up/(0.5)**(2*par_N_up)+par_R
        def shape_fcn_Y_low(XI): return par_hight_low*(0.5+XI*(psi_end_cut))**par_N_low*(
            0.5-XI*(psi_end_cut))**par_N_low/(0.5)**(2*par_N_low)+par_R

        # local coordinate, local x is global y, local y is global z
        body_back_up = CST3DSurface(y_cut, 1, 0, [], shape_fcn_Y_up)
        body_back_low = CST3DSurface(y_cut, - 1, 0, [], shape_fcn_Y_low)
        # rotation to global coordinatelambda XI:
        body_back_up.addRotation(90, 90, 0)
        body_back_low.addRotation(90, 90, 0)
        # translation to global coordination
        body_back_up.addTranslation(total_length, 0, 0)
        body_back_low.addTranslation(total_length, 0, 0)

        # define transition wing
        # local coordination, local x is global y, local y is global -x, local z is global z
        def shape_fcn_Y(XI): return (1-XI*(1-par_rho12))

        def shape_body_size_up(PSI): return body_up.LZ * \
            body_up.shape_fcn_Z(1-PSI, 1)*body_up.class_fcn_Z(1-PSI, 1)

        def shape_body_size_low(PSI): return -body_low.LZ * \
            body_low.shape_fcn_Z(1-PSI, 1)*body_low.class_fcn_Z(1-PSI, 1)

        def shape_fcn_Z_up(XI, PSI): return (
            1-XI)*shape_body_size_up(PSI)+XI*self.shapeWing(1-PSI)*(wing_height_up)

        def shape_fcn_Z_low(XI, PSI): return (
            1-XI)*shape_body_size_low(PSI)+XI*self.shapeWing(1-PSI)*(wing_height_low)
        tri_wing_up = CST3DSurface(tri_wing_width, tri_wing_length, 1, [
        ], shape_fcn_Y, shape_fcn_Z_up, [0, 0])
        tri_wing_low = CST3DSurface(
            tri_wing_width, tri_wing_length, - 1, [], shape_fcn_Y, shape_fcn_Z_low, [0, 0])

        # add rotation
        tri_wing_up.addRotation(0, 0, 90)
        tri_wing_low.addRotation(0, 0, 90)
        # add translation
        tri_wing_up.addTranslation(total_length, y_cut, par_R)
        tri_wing_low.addTranslation(total_length, y_cut, - par_R)

        # define blunt tir wing front
        # local coordinate, local x is global -y, local y is global z, local z is global -x
        # if do not want blunt intersection, remove self and change blunt head side class N1
        def shape_fcn_X(PSI): return 1-(sqrt(radius_head_side_sq-(((PSI-0.5)*2)
                                                                  * par_R)**2)-radius_center_head_side)/tri_wing_width

        def shape_fcn_Z(XI, PSI): return (1-XI)*shape_wing_edge_x((PSI-0.5)
                                                                  * 2*par_R)+XI*shape_tri_wing_edge_x((PSI-0.5)*2*par_R)
        tri_wing_front = CST3DSurface(
            tri_wing_width, 2 * par_R, 1, shape_fcn_X, [], [], shape_fcn_Z)
        # deform surface

        def tran_fun_Z(XI, PSI): return XI * \
            (par_rho1*(1-par_rho12)*total_length)
        tri_wing_front.addDeform([], [], tran_fun_Z)
        # rotation surface
        tri_wing_front.addRotation(90, - 90, 0)
        # translation surface
        tri_wing_front.addTranslation(
            total_length - wing_length, y_cut + par_WS1, - par_R)

        # define tri wing back
        # local coordination, local x is global y, local y is global z
        def shape_fcn_Y_up(XI): return (
            1-XI)*(tri_wing_height_up+par_R)+XI*(par_R)
        def shape_fcn_Y_low(XI): return (
            1-XI)*(tri_wing_height_low+par_R)+XI*(par_R)
        tri_wing_back_up = CST3DSurface(
            tri_wing_width, 1, 0, [], shape_fcn_Y_up)
        tri_wing_back_low = CST3DSurface(
            tri_wing_width, - 1, 0, [], shape_fcn_Y_low)

        # rotation to global coordinate
        tri_wing_back_up.addRotation(90, 90, 0)
        tri_wing_back_low.addRotation(90, 90, 0)
        # translation to global coordination
        tri_wing_back_up.addTranslation(total_length, y_cut, 0)
        tri_wing_back_low.addTranslation(total_length, y_cut, 0)

        # define wing
        # local coordination, local x is global y, local y is global -x
        def shape_fcn_Z(XI, PSI): return self.shapeWing(1-PSI)
        wing_up = CST3DSurface(wing_width, wing_length, wing_height_up, [
        ], [], shape_fcn_Z, [[0, 0], [0, 0]])
        wing_low = CST3DSurface(
            wing_width, wing_length, - wing_height_low, [], [], shape_fcn_Z, [[0, 0], [0, 0]])
        # add rotation
        wing_up.addRotation(0, 0, 90)
        wing_low.addRotation(0, 0, 90)
        # add translation
        wing_up.addTranslation(total_length, y_cut + par_WS1, par_R)
        wing_low.addTranslation(total_length, y_cut + par_WS1, - par_R)
        # define blunt wing front

        # local coordinate, local x is global -y, local y is global z, local z is global -x
        def shape_fcn_Z(XI, PSI): return shape_wing_edge_x((PSI-0.5) * 2*par_R)
        wing_front = CST3DSurface(
            wing_width, 2 * par_R, 1, [], [], shape_fcn_Z, [0, 0])
        # rotation surface
        wing_front.addRotation(90, - 90, 0)
        # translation surface
        wing_front.addTranslation(
            total_length - wing_length, y_cut + par_WS1 + par_WS2, - par_R)

        # define wing back
        # local coordination, local x is global y, local y is global z
        wing_back_up = CST3DSurface(wing_width, par_R, 0)
        wing_back_low = CST3DSurface(wing_width, - par_R, 0)
        # rotation to global coordinate
        wing_back_up.addRotation(90, 90, 0)
        wing_back_low.addRotation(90, 90, 0)
        # translation to global coordination
        wing_back_up.addTranslation(total_length, y_cut + par_WS1, 0)
        wing_back_low.addTranslation(total_length, y_cut + par_WS1, 0)

        # define wing side
        # local coordinate, local x is global -x, local y is global z
        def shape_fcn_Y_up(PSI): return self.shapeWing(1-PSI)*par_WT_up+par_R
        def shape_fcn_Y_low(PSI): return self.shapeWing(1-PSI)*par_WT_low+par_R
        wing_side_up = CST3DSurface(wing_length, 1, 0, [], shape_fcn_Y_up)
        wing_side_low = CST3DSurface(wing_length, - 1, 0, [], shape_fcn_Y_low)

        # wing side blunt front
        def shape_fcn_X(PSI): return shape_wing_edge_x((PSI-0.5)*2*par_R)
        wing_side_front = CST3DSurface(1, 2 * par_R, 0, shape_fcn_X)
        wing_side_up.addRotation(90, 180, 0)
        wing_side_low.addRotation(90, 180, 0)
        wing_side_front.addRotation(90, 180, 0)
        # translation to global coordinate
        wing_side_up.addTranslation(total_length, y_cut + par_WS1 + par_WS2, 0)
        wing_side_low.addTranslation(
            total_length, y_cut + par_WS1 + par_WS2, 0)
        wing_side_front.addTranslation(
            total_length - wing_length, y_cut + par_WS1 + par_WS2, - par_R)
        # sort data
        self.surface_list = {'head_up': head_up, 'head_low': head_low, 'head_side': head_side,
                             'body_up': body_up, 'body_low': body_low, 'body_back_up': body_back_up, 'body_back_low': body_back_low,
                             'tri_wing_up': tri_wing_up, 'tri_wing_low': tri_wing_low, 'tri_wing_front': tri_wing_front,
                             'tri_wing_back_up': tri_wing_back_up, 'tri_wing_back_low': tri_wing_back_low,
                             'wing_up': wing_up, 'wing_low': wing_low, 'wing_front': wing_front,
                             'wing_back_up': wing_back_up, 'wing_back_low': wing_back_low,
                             'wing_side_up': wing_side_up, 'wing_side_low': wing_side_low, 'wing_side_front': wing_side_front}

    def calSurfaceMatrix(self, xi_grid_num_head: int, psi_grid_num_head: int, xi_grid_num_body: int, psi_grid_num_wing: int, edge_gird_num: int):
        # generate waverider wing wgs data
        # xi, x=X(xi)
        # psi, y=Y(xi)
        # zeta, z=Z(xi,psi)

        # notice colume vector in matrix is wgs format

        # auto allocation xi
        head_xi_list=np.linspace(0, 1, xi_grid_num_head + 1)
        head_xi_list=head_xi_list**(1/self.par_T)
        # head_xi_list = np.array(
        #     [0, self.first_xi]+list(np.linspace(self.stag_xi, 1, xi_grid_num_head - 1)))
        
        # calculate head
        XI, PSI = np.meshgrid(head_xi_list, np.linspace(
            0.5, 1, psi_grid_num_head + 1))
        X_head, Y_head, Z_head_up = self.surface_list['head_up'].calSurface(
            XI, PSI)
        X_head, Y_head, Z_head_low = self.surface_list['head_low'].calSurface(
            XI, PSI)

        # # calculate stag
        # stag_xi_list = [0, (self.first_xi*(1-self.par_rho1))
        #                 ** self.par_T*self.par_width/2/self.stag_length, 1]
        # [XI, ZETA] = np.meshgrid(
        #     stag_xi_list, np.linspace(0, 1, edge_gird_num*2+1))
        # X_stag, Y_stag, Z_stag = self.surface_list['stag'].calSurface(XI, ZETA)

        # calculate blunt head side
        head_side_xi_list=1-head_xi_list
        XI, ZETA = np.meshgrid(
            head_xi_list, np.linspace(0, 1, edge_gird_num * 2 + 1))
        X_head_side, Y_head_side, Z_head_side = self.surface_list['head_side'].calSurface(
            XI, ZETA)

        # calculate body
        X_body, Y_body, Z_body_up = self.surface_list['body_up'].calSurface(
            xi_grid_num_body, psi_grid_num_head)
        X_body, Y_body, Z_body_low = self.surface_list['body_low'].calSurface(
            xi_grid_num_body, psi_grid_num_head)
        # calculate body back
        X_body_back, Y_body_back, Z_body_back_up = self.surface_list['body_back_up'].calSurface(
            psi_grid_num_head, edge_gird_num)
        X_body_back, Y_body_back, Z_body_back_low = self.surface_list['body_back_low'].calSurface(
            psi_grid_num_head, edge_gird_num)
        # calculate transition wing
        X_tri_wing, Y_tri_wing, Z_tri_wing_up = self.surface_list['tri_wing_up'].calSurface(
            psi_grid_num_wing, xi_grid_num_body)
        X_tri_wing, Y_tri_wing, Z_tri_wing_low = self.surface_list['tri_wing_low'].calSurface(
            psi_grid_num_wing, xi_grid_num_body)
        # calculate blunt tir wing front
        X_tri_wing_front, Y_tri_wing_front, Z_tri_wing_front = self.surface_list['tri_wing_front'].calSurface(
            psi_grid_num_wing, 2 * edge_gird_num)
        # calculate tri wing back
        Y_tri_wing_back, Y_tri_wing_back, Z_tri_wing_back_up = self.surface_list['tri_wing_back_up'].calSurface(
            psi_grid_num_wing, edge_gird_num)
        X_tri_wing_back, Y_tri_wing_back, Z_tri_wing_back_low = self.surface_list['tri_wing_back_low'].calSurface(
            psi_grid_num_wing, edge_gird_num)
        # calculate wing
        X_wing, Y_wing, Z_wing_up = self.surface_list['wing_up'].calSurface(
            psi_grid_num_wing, xi_grid_num_body)
        X_wing, Y_wing, Z_wing_low = self.surface_list['wing_low'].calSurface(
            psi_grid_num_wing, xi_grid_num_body)
        # calculate blunt wing front
        X_wing_front, Y_wing_front, Z_wing_front = self.surface_list['wing_front'].calSurface(
            psi_grid_num_wing, 2 * edge_gird_num)
        # calculate wing back
        X_wing_back, Y_wing_back, Z_wing_back_up = self.surface_list['wing_back_up'].calSurface(
            psi_grid_num_wing, edge_gird_num)
        X_wing_back, Y_wing_back, Z_wing_back_low = self.surface_list['wing_back_low'].calSurface(
            psi_grid_num_wing, edge_gird_num)
        # calculate wing side
        X_wing_side, Y_wing_side, Z_wing_side_up = self.surface_list['wing_side_up'].calSurface(
            xi_grid_num_body, edge_gird_num)
        X_wing_side, Y_wing_side, Z_wing_side_low = self.surface_list['wing_side_low'].calSurface(
            xi_grid_num_body, edge_gird_num)
        X_wing_side_front, Y_wing_side_front, Z_wing_side_front = self.surface_list['wing_side_front'].calSurface(
            edge_gird_num, 2 * edge_gird_num)

        # sort data
        X = {}
        Y = {}
        Z = {}
        X, Y, Z = self.sortSurfaceUpLow(
            X, Y, Z, 'head', X_head, Y_head, Z_head_up, Z_head_low)
        X, Y, Z = self.sortSurface(
            X, Y, Z, 'head_side', X_head_side, Y_head_side, Z_head_side)
        X, Y, Z = self.sortSurfaceUpLow(
            X, Y, Z, 'body', X_body, Y_body, Z_body_up, Z_body_low)
        X, Y, Z = self.sortSurfaceUpLow(
            X, Y, Z, 'body_back', X_body_back, Y_body_back, Z_body_back_up, Z_body_back_low)
        X, Y, Z = self.sortSurfaceUpLow(
            X, Y, Z, 'tri_wing', X_tri_wing, Y_tri_wing, Z_tri_wing_up, Z_tri_wing_low)
        X, Y, Z = self.sortSurface(
            X, Y, Z, 'tri_wing_front', X_tri_wing_front, Y_tri_wing_front, Z_tri_wing_front)
        X, Y, Z = self.sortSurfaceUpLow(
            X, Y, Z, 'tri_wing_back', X_tri_wing_back, Y_tri_wing_back, Z_tri_wing_back_up, Z_tri_wing_back_low)
        X, Y, Z = self.sortSurfaceUpLow(
            X, Y, Z, 'wing', X_wing, Y_wing, Z_wing_up, Z_wing_low)
        X, Y, Z = self.sortSurface(
            X, Y, Z, 'wing_front', X_wing_front, Y_wing_front, Z_wing_front)
        X, Y, Z = self.sortSurfaceUpLow(
            X, Y, Z, 'wing_back', X_wing_back, Y_wing_back, Z_wing_back_up, Z_wing_back_low)
        X, Y, Z = self.sortSurfaceUpLow(
            X, Y, Z, 'wing_side', X_wing_side, Y_wing_side, Z_wing_side_up, Z_wing_side_low)
        X, Y, Z = self.sortSurface(
            X, Y, Z, 'wing_side_front', X_wing_side_front, Y_wing_side_front, Z_wing_side_front)
        return X, Y, Z

    def shapeWing(self, XI=None):
        # height normalize to 1

        PSI = 1 - np.abs(1 - 2 * XI)
        return PSI

    def sortSurface(self, X, Y, Z, name, X_surf, Y_surf, Z_surf):
        X[name] = X_surf
        Y[name] = Y_surf
        Z[name] = Z_surf
        return X, Y, Z

    def sortSurfaceUpLow(self, X, Y, Z, name, X_surf, Y_surf, Z_surf_up, Z_surf_low):
        X[name+'_up'] = X_surf
        Y[name+'_up'] = Y_surf
        Z[name+'_up'] = Z_surf_up
        X[name+'_low'] = X_surf
        Y[name+'_low'] = Y_surf
        Z[name+'_low'] = Z_surf_low
        return X, Y, Z

    def getWGSMesh(self, part_name, xi_grid_num_head: int, psi_grid_num_head: int, xi_grid_num_body: int, psi_grid_num_wing: int, edge_gird_num: int):
        surface_name_list = list(self.surface_list.keys())
        surface_number = len(surface_name_list)
        X_total, Y_total, Z_total = self.calSurfaceMatrix(
            xi_grid_num_head, psi_grid_num_head, xi_grid_num_body, psi_grid_num_wing, edge_gird_num)
        mesh_list = []
        for surface_index in range(surface_number):
            surface_name = surface_name_list[surface_index]
            mesh = {}
            if ('up' in surface_name) or ('front' in surface_name) or ('head_side' in surface_name):
                mesh['X'] = X_total[surface_name]
                mesh['Y'] = Y_total[surface_name]
                mesh['Z'] = Z_total[surface_name]
            else:
                mesh['X'] = fliplr(X_total[surface_name])
                mesh['Y'] = fliplr(Y_total[surface_name])
                mesh['Z'] = fliplr(Z_total[surface_name])
            mesh['element_type'] = 'wgs'
            mesh_list.append(mesh)

        if not self.half_flag:
            for surface_index in range(surface_number):
                surface_name = surface_name_list[surface_index]
                mesh = {}
                if ('up' in surface_name) or ('front' in surface_name) or ('head_side' in surface_name):
                    mesh['X'] = fliplr(X_total[surface_name])
                    mesh['Y'] = -fliplr(Y_total[surface_name])
                    mesh['Z'] = fliplr(Z_total[surface_name])
                else:
                    mesh['X'] = X_total[surface_name]
                    mesh['Y'] = -Y_total[surface_name]
                    mesh['Z'] = Z_total[surface_name]
                mesh['element_type'] = 'wgs'
                mesh_list.append(mesh)

        part = {}
        part['name'] = part_name
        part['mesh_list'] = mesh_list
        return part

    def getWGS(self, part_name, xi_grid_num_head: int, psi_grid_num_head: int, xi_grid_num_body: int, psi_grid_num_wing: int, edge_gird_num: int):
        this_wgs = wgs_creator.LaWGS(part_name)

        X_total, Y_total, Z_total = self.calSurfaceMatrix(
            xi_grid_num_head, psi_grid_num_head, xi_grid_num_body, psi_grid_num_wing, edge_gird_num)

        surface_name_list = list(self.surface_list.keys())
        surface_number = len(surface_name_list)

        for surface_index in range(surface_number):
            surface_name = surface_name_list[surface_index]

            if ('up' in surface_name) or ('front' in surface_name):
                X = X_total[surface_name]
                Y = Y_total[surface_name]
                Z = Z_total[surface_name]
            else:
                X = fliplr(X_total[surface_name])
                Y = fliplr(Y_total[surface_name])
                Z = fliplr(Z_total[surface_name])

            line_num = len(X[0])
            point_num = len(X)
            network = list()
            for line_index in range(line_num):
                line = np.zeros((point_num, 3))
                line[:, 0] = X[:, line_index]
                line[:, 1] = Y[:, line_index]
                line[:, 2] = Z[:, line_index]
                network.append(line)

            # convert to network
            this_network = wgs_creator.Network(network)
            this_wgs.append_network(surface_name, this_network, 1)

        if not self.half_flag:
            for surface_index in range(surface_number):
                surface_name = surface_name_list[surface_index]

                if ('up' in surface_name) or ('front' in surface_name):
                    X = fliplr(X_total[surface_name])
                    Y = -fliplr(Y_total[surface_name])
                    Z = fliplr(Z_total[surface_name])
                else:
                    X = X_total[surface_name]
                    Y = -Y_total[surface_name]
                    Z = Z_total[surface_name]

                line_num = len(X[0])
                point_num = len(X)
                network = list()
                for line_index in range(line_num):
                    line = np.zeros((point_num, 3))
                    line[:, 0] = X[:, line_index]
                    line[:, 1] = Y[:, line_index]
                    line[:, 2] = Z[:, line_index]
                    network.append(line)

                # convert to network
                this_network = wgs_creator.Network(network)
                this_wgs.append_network(surface_name+'_sys', this_network, 1)

        return this_wgs

    def writeWGS(self, part_name, xi_grid_num_head: int, psi_grid_num_head: int, xi_grid_num_body: int, psi_grid_num_wing: int, edge_gird_num: int):
        this_wgs = self.getWGS(part_name, xi_grid_num_head, psi_grid_num_head,
                               xi_grid_num_body, psi_grid_num_wing, edge_gird_num)
        this_wgs.create_wgs(filename=part_name+'.wgs')

    def writeSTL(self, part_name, xi_grid_num_head: int, psi_grid_num_head: int, xi_grid_num_body: int, psi_grid_num_wing: int, edge_gird_num: int):
        this_wgs = self.getWGS(part_name, xi_grid_num_head, psi_grid_num_head,
                               xi_grid_num_body, psi_grid_num_wing, edge_gird_num)
        this_wgs.create_stl(filename=part_name+'.stl')

    def drawBody(self, xi_grid_num_head: int, psi_grid_num_head: int, xi_grid_num_body: int, psi_grid_num_wing: int, edge_gird_num: int):
        surface_name_list = list(self.surface_list.keys())
        surface_number = len(surface_name_list)

        X_total, Y_total, Z_total = self.calSurfaceMatrix(
            xi_grid_num_head, psi_grid_num_head, xi_grid_num_body, psi_grid_num_wing, edge_gird_num)

        fig = plt.figure()
        axe = Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(axe)
        for surface_index in range(surface_number):
            surface_name = surface_name_list[surface_index]
            axe.plot_surface(X_total[surface_name], Y_total[surface_name],
                             Z_total[surface_name], rstride=1, cstride=1, cmap='Greys')
        
        axe.set_xlabel(r'$x$', fontsize=18)
        axe.set_ylabel(r'$y$', fontsize=18)
        axe.set_zlabel(r'$z$', fontsize=18)
        axe.view_init(30, -120)
        plt.axis('auto')
        plt.show()

        return


if __name__ == '__main__':
    # print("\nDemo: Waverider Geometry Generation using CST and pyPanair")
    total_length = 4

    variable_num = 15
    low_bou = np.array([0.5, 0.5, 2.4, 0.5, 0.5, 0.5, 0.1,
                       0.1, 0.005, 0.25, 0.6, 0.6, 0.5, 0.01, 0.01])
    up_bou = np.array([1.0, 1.0, 3.2, 0.9, 2.0, 2.0, 0.5,
                      0.5, 0.02, 0.35, 0.8, 0.8, 0.7, 0.05, 0.05])
    x=np.random.random(variable_num)*(up_bou-low_bou)+low_bou
    # x = (up_bou + low_bou) / 2
    # x=up_bou
    # x=low_bou

    # data = io.loadmat('LHD.mat')
    # X = data['X_HF']
    # x = X[0]

    # x=[0.5,0.7,...
    #     2.4,0.5,2,2,...
    #     0.2,0.5,0.005,...
    #     0.3,0.6,0.6,0.6,0.01,0.02]

    # length parameter
    par_M_up = x[0]
    par_M_low = x[1]
    # width parameter
    par_W = x[2]
    par_T = x[3]
    par_NB_up = x[4]
    par_NB_low = x[5]
    # height parameter
    par_hight_up = x[6]
    par_hight_low = x[7]
    par_R = x[8]
    par_R = 0.5
    # wing parameter
    par_rho1 = x[9]
    par_rho12 = x[10]
    par_WS1 = x[11]
    par_WS2 = x[12]
    par_TWU = x[13]
    par_TWL = x[14]

    # create waverider

    # waverider_wing=WaveriderWingFei...
    #     [total_length,par_W,par_hight_up,par_hight_low,...
    #     par_T,par_M_up,par_M_low,par_NB_up,par_NB_low,par_R,...
    #     par_rho1,par_rho12,par_WS1,par_WS2)

    # waverider_wing=WaveriderWingTri...
    #     (total_length,par_W,par_hight_up,par_hight_low,...
    #     par_T,par_M_up,par_M_low,par_NB_up,par_NB_low,par_R,...
    #     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL)

    waverider_wing = WaveriderWingDia(total_length, par_W, par_hight_up, par_hight_low, par_T, par_M_up,
                                      par_M_low, par_NB_up, par_NB_low, par_R, par_rho1, par_rho12, par_WS1, par_WS2, par_TWU, par_TWL, True)
    # draw

    xi_grid_num_head = 30
    eta_grid_num_head = 20
    xi_grid_num_body = 20
    eta_grid_num_wing = 6

    edge_gird_num = int(np.ceil((par_R / 0.003) / 2) * 2)

    # waverider_wing.drawBody(xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num)

    # write to wgs
    # waverider_wing.writeWGS('waverider_wing_dia', xi_grid_num_head,
    # eta_grid_num_head, xi_grid_num_body, eta_grid_num_wing, edge_gird_num)
    waverider_wing.writeSTL('waverider_wing_dia', xi_grid_num_head,
                            eta_grid_num_head, xi_grid_num_body, eta_grid_num_wing, edge_gird_num)
