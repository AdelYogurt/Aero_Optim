import numpy as np
from numpy import array, zeros, argmin, sqrt, pi, sin, cos, arcsin, arctan, log
import sys
import os


def infStagHeat(R_min, R_max=None, Ma_1=None, P_1=None, T_1=None, T_w=300, H=40, method="Fay-Riddell", mode='CST'):
    """
    BWB thermodynamics simulation using engineering prediction methods.

    Inputs:
    ------
    pntCST: np.ndarray [12, ]
        CST parameters for BWB, including:
        varCoefW, varClsT, varCoefHu, varCoefHl, varClsNu, varClsNl, varClsM1u, varClsM1l: coefficients in CST
        varRtoHead: ratio of head length
        varRtoChrd: ratio of chord length of 2nd section respect to root chord length
        varSSpan1, varSSpan2: semi-span of 1st and 2nd sections
    dictCndtn: dict, flight condition
        mach: mach number
        alpha: angle of attack (deg)
        alt: altitude (km)
    method: str, method for stagnation point heat flux prediction, including:
        `Fay-Riddell`, `Kemp-Riddell`, `Scala`, `Lees`, `Detra-Kemp-Riddell`, `Tauber`.

    Output:
    ------
    q: float
        maximum heat flux, among that at head point, wing1/wing2 leading edge (laminar/turbulence heat flux)
    """

    if method not in ["Fay-Riddell", "Kemp-Riddell", "Scala", "Lees", "Detra-Kemp-Riddell", "Tauber"]:
        raise RuntimeError(
            "Currently supported methods include only Fay-Riddell, Kemp-Riddell, Scala, Lees, Detra-Kemp-Riddell, and Tauber.")

    # ----------------------Head curvature radius estimate-----------------------
    if R_max is not None:
        R_S = sqrt((1 + R_min / R_max) / 2) * R_min
    else:
        R_max = R_min
        R_S = R_min

    # ------------------------------air parameter---------------------------------
    R = 287.0955
    rho_1 = P_1/R/T_1
    k = pred_k(T_1)
    a_1 = sqrt(R*k*T_1)
    q_1 = 0.5*rho_1*(a_1*Ma_1)**2

    Ma_1_sq = Ma_1*Ma_1
    Cp_stag = 2/k/Ma_1_sq*(((k+1)**2*Ma_1_sq/(4*k*Ma_1_sq-2*(k-1)))
                           ** (k/(k-1))*((1-k+2*k*Ma_1_sq)/(k+1))-1)
    P_w = Cp_stag*q_1+P_1

    # ------------------------------Head heat flux--------------------------------
    if method == "Fay-Riddell":
        def fcnSphere(R): return Fay_Riddell(
            R, Ma_1, rho_1, P_1, T_1, P_w, T_w)
    elif method == "Kemp-Riddell":
        def fcnSphere(R): return Kemp_Riddell(
            R, Ma_1, rho_1, P_1, T_1, P_w, T_w)
    elif method == "Scala":
        def fcnSphere(R): return Scala(R, Ma_1, rho_1, P_1, T_1, P_w, T_w, H)
    elif method == "Lees":
        def fcnSphere(R): return Lees(R, Ma_1, rho_1, P_1, T_1, P_w, T_w)
    elif method == "Detra-Kemp-Riddell":
        def fcnSphere(R): return Detra_Kemp_Riddell(
            R, Ma_1, rho_1, P_1, T_1, P_w, T_w, R_min, R_max)
    elif method == "Tauber":
        def fcnSphere(R): return Tauber(R, Ma_1, rho_1, P_1, T_1, P_w, T_w)
    HF = fcnSphere(R_S)

    return HF


def Fay_Riddell(RN, Ma_1, rho_1, P_1, T_1, P_w, T_w):
    """
    Stagnation point heat flux estimation in zero angle of attack: Fay-Riddell function.

    Inputs:
    ------
    RN: curvature radius of stagnation point
    Ma_1: Mach number of free stream
    rho_1: density of free stream
    P_1: pressure of free stream
    T_1: temperature of free stream

    Output:
    ------
    qs: stagnation point heat flux

    References:
    ----------
    [1] FAY J A, RIDDELL F R. Theory of stagnation point heat transfer in dissociated air[J]. Journal of the Aerospace Sciences, 1958, 25(2): 73-85.
    [2] ILICH Z, GROSSIR G, CHAZOT O. Evaluation of the stagnation-point velocity gradient in low-enthalpy hypersonic flows[M]. 33rd AIAA Aerodynamic Measurement Technology and Ground Testing Conference. American Institute of Aeronautics and Astronautics. 2017.
    [3] 张志成. 高超声速气动热和热防护[M]. 北京: 国防工业出版社, 2003.
    """

    # constants
    Pr = 0.71  # Prandtl constant
    Le = 1.4  # Lewis number, 1.4 for continuous flow
    AOA = 0.52  # power of Lewis number
    mu0 = 1.716E-5
    T0 = 273.15  # reference viscosity and temperature in Sutherland law
    Smu = 110.4  # Sutherland constant for standard air
    hO = 1.5420E7
    hN = 3.3867E7  # unit mass dissociation enthalpy of O and N (J/kg)
    COs = 0.20946
    CNs = 0.78084  # mole fraction of O and N in atmosphere
    W = 0.02897  # averaged molecular mass of air
    R = 8.31446  # ideal gas constant (J/(K·mol))

    # intermidiate parameters
    V_1 = Ma_1 * sqrt(pred_k(T_1) * 287.05 * T_1)  # freestream velocity
    rho_2, P_2, T_2, rho_es, P_es, T_es = calShockThermo(rho_1, P_1, T_1, Ma_1)
    H_w = pred_Cp(T_w)*T_w
    rho_w = airDensity(H_w, P_w)
    muw = mu0 * (T_w / T0) ** 1.5 * (T0 + Smu) / (T_w + Smu)  # Sutherland law
    mus = mu0 * (T_es / T0) ** 1.5 * (T0 + Smu) / \
        (T_es + Smu)  # Sutherland law
    due_ds = 1/RN * sqrt(2 * (P_es - P_w) / rho_es)  # modified Newton theory
    # averaged dissociation enthalpy of air
    hd = (COs * hO + CNs * hN) / (COs + CNs)
    # sensible enthalpy of air at stagnation point
    H_es = 1.003 * pred_Cp(T_es) * T_es
    # hw = R/W * (3.51715*T_w - 2.5041E-4*T_w**2 + 5.5079E-7*T_w**3 - 1.7197E-10*T_w**4) # enthalpy of air at wall points
    # hw = pred_Cp(T_1) * T_w # enthalpy of air at wall points
    H_w = pred_Cp(T_w) * T_w  # enthalpy of air at wall points

    qs = 0.763 * Pr**(-0.6) * (rho_w*muw / (rho_es*mus))**0.1 * \
        sqrt(rho_es*mus * due_ds) * \
        (1 + (Le**AOA - 1) * hd/H_es) * (H_es - H_w)

    return qs


def Kemp_Riddell(RN, Ma_1, rho_1, P_1, T_1, P_w, T_w):
    """
    Stagnation point heat flux estimation in zero angle of attack: Kemp-Riddell function.

    Inputs:
    ------
    RN: curvature radius of stagnation point
    Ma_1: Mach number of free stream
    rho_1: density of free stream
    P_1: pressure of free stream
    T_1: temperature of free stream

    Output:
    ------
    qs: stagnation point heat flux

    Reference:
    ---------
    [1] KEMP N H, RIDDELL F. Heat transfer to satellite vehicles re-entering the atmosphere[J]. Journal of Jet Propulsion, 1957, 27(2): 132-137.
    [2] 杜涛, 陈闽慷, 李凰立, 等. 高超声速气动加热关联方法的适应性分析[J]. 宇航学报, 2018, 39(09): 1039-1046.
    [3] 张志成. 高超声速气动热和热防护[M]. 北京: 国防工业出版社, 2003.
    """

    # constants
    rhosl = 1.225  # air density at sea level
    vc = 7925  # reference velocity
    W = 0.02897  # averaged molecular mass of air
    R = 8.31446  # ideal gas constant (J/(K·mol))

    # intermediate parameters
    V_1 = Ma_1 * sqrt(pred_k(T_1) * 287.05 * T_1)  # freestream velocity
    rho_2, P_2, T_2, rho_es, P_es, T_es = calShockThermo(rho_1, P_1, T_1, Ma_1)
    # hs = 1.003 * pred_Cp(T_1) * T_1 + V_1**2 / 2 # sensible enthalpy of air at stagnation point
    # hw = R/W * (3.51715*T_w - 2.5041E-4*T_w**2 + 5.5079E-7*T_w**3 - 1.7197E-10*T_w**4) # enthalpy of air at wall points
    # h300K = R/W * (3.51715*300 - 2.5041E-4*300**2 + 5.5079E-7*300**3 - 1.7197E-10*300**4) # enthalpy of air at wall points when T_w=300(K)
    # sensible enthalpy of air at stagnation point
    H_es = pred_Cp(T_es) * T_es
    H_w = pred_Cp(T_w) * T_w  # enthalpy of air at wall points
    # enthalpy of air at wall points when T_w=300(K)
    H_300K = pred_Cp(300) * 300

    qs = 1.318842E8 / sqrt(RN) * sqrt(rho_1/rhosl) * \
        (V_1/7900)**3.25 * (H_es - H_w) / (H_es - H_300K)

    return qs


def Scala(RN, Ma_1, rho_1, P_1, T_1, P_w, T_w, H):
    """
    #! deprecated
    Stagnation point heat flux estimation in zero angle of attack: Scala function.

    Inputs:
    ------
    RN: curvature radius of stagnation point
    Ma_1: Mach number of free stream
    rho_1: density of free stream
    P_1: pressure of free stream
    T_1: temperature of free stream
    H: altitude (m)

    Output:
    ------
    qs: stagnation point heat flux

    Reference:
    ---------
    [1] SCALA S M. A study of hypersonic ablation/untersuchung der hyperschall-ablation/etude de l’ablation en régime hypersonique Proceedings of the Xth International Astronautical Congress London 1959/X Internationaler Astronautischer Kongress/Xe Congrès International d’Astronautique, F, 1960[C]. Springer.
    [2] PARK S-H, NEEB D, PLYUSHCHEV G, et al. A study on heat flux predictions for re-entry flight analysis[J]. Acta Astronautica, 2021, 187: 271-280.
    [3] 张志成. 高超声速气动热和热防护[M]. 北京: 国防工业出版社, 2003.
    """

    # intermediate parameters
    V_1 = Ma_1 * sqrt(pred_k(T_1) * 287.05 * T_1)  # freestream velocity
    rho_2, P_2, T_2, rho_es, P_es, T_es = calShockThermo(rho_1, P_1, T_1, Ma_1)
    # a = - (0.9689 + 6.9984E-5 * T_es) * (5.626 + 3.2285E-5 * H)
    # b = (0.9793 + 4.6715E-5 * T_es) * (2.838 + 9.843E-7 * H)
    a = - (0.9689 + 6.998E-5 * T_w) * (5.626 + 9.84E-6 * 3.28*H)
    b = (0.9793 + 4.672E-5 * T_w) * (2.830 + 3E-7 * 3.28*H)

    # qs = 12488 / sqrt(RN) * 10**a * (3.281E-3 * V_1)**b
    qs = 0.564 / sqrt(RN*2*3.28) * 10**a * (V_1*3.28/10)**b * 5.6783*3600

    return qs


def Lees(RN, Ma_1, rho_1, P_1, T_1, P_w, T_w):
    """
    #! deprecated
    Stagnation point heat flux estimation in zero angle of attack: Lees function.

    Inputs:
    ------
    RN: curvature radius of stagnation point
    Ma_1: Mach number of free stream
    rho_1: density of free stream
    P_1: pressure of free stream
    T_1: temperature of free stream

    Output:
    ------
    qs: stagnation point heat flux

    Reference:
    ---------
    [1] LEES L. Laminar heat transfer over blunt-nosed bodies at hypersonic flight speeds[J]. Journal of Jet Propulsion, 1956, 26(4): 259-269.
    [2] PARK S-H, NEEB D, PLYUSHCHEV G, et al. A study on heat flux predictions for re-entry flight analysis[J]. Acta Astronautica, 2021, 187: 271-280.
    """

    # constants
    Pr = 0.71  # Prandtl constant
    mu0 = 1.716E-5
    T0 = 273.15  # reference viscosity and temperature in Sutherland law
    Smu = 110.4  # Sutherland constant for standard air
    gamma = 1.15  # [1.1, 1.2] at high temperatures

    # intermediate parameters
    k = pred_k(T_1)
    V_1 = Ma_1 * sqrt(k * 287.05 * T_1)  # freestream velocity
    rho_2, P_2, T_2, rho_es, P_es, T_es = calShockThermo(rho_1, P_1, T_1, Ma_1)
    mus = mu0 * (T_es / T0) ** 1.5 * (T0 + Smu) / \
        (T_es + Smu)  # Sutherland law
    # hs = 1.003 * pred_Cp(T_1) * T_1 + V_1**2 / 2 # sensible enthalpy of air at stagnation point
    hs = V_1**2 / 2
    G = ((gamma-1)/gamma)**0.25 * (1 + 2 / ((k-1) * Ma_1**2))**0.25 * \
        (1 - 1/(k*Ma_1**2))**0.25  # intermediate variable in Lees equation

    qs = 1/sqrt(RN) * sqrt(2)/2 / Pr**(2/3) * sqrt(rho_es * mus * V_1) * hs * G

    return qs


def Detra_Kemp_Riddell(RN, Ma_1, rho_1, P_1, T_1, P_w, T_w, Rmin, Rmax):
    """
    Stagnation point heat flux estimation in zero angle of attack: Detra-Kemp-Riddell function.

    Inputs:
    ------
    RN: curvature radius of stagnation point (act as placeholder)
    Ma_1: Mach number of free stream
    rho_1: density of free stream
    P_1: pressure of free stream
    T_1: temperature of free stream
    Rmin: smaller principal curvature radius of stagnation point
    Rmax: larger principal curvature radius of stagnation point

    Output:
    ------
    qs: stagnation point heat flux

    Reference:
    [1] CRABTREE L, DOMMETT R, WOODLEY J. Estimation of heat transfer to flat plates, cones and blunt bodies[R]. Reports and Memoranda No. 3637, 1965.
    [2] PARK S-H, NEEB D, PLYUSHCHEV G, et al. A study on heat flux predictions for re-entry flight analysis[J]. Acta Astronautica, 2021, 187: 271-280.
    ---------
    """

    # constants
    rho_sl = 1.225  # air density at sea level

    # intermediate parameters
    k = pred_k(T_1)
    V_1 = Ma_1 * sqrt(k * 287.05 * T_1)  # freestream velocity
    rho_2, P_2, T_2, rho_es, P_es, T_es = calShockThermo(rho_1, P_1, T_1, Ma_1)
    # sensible enthalpy of air at stagnation point
    # H_es = pred_Cp(T_1) * T_1 + V_1**2 / 2
    H_es = pred_Cp(T_es) * T_es
    H_w = pred_Cp(T_1) * T_w  # enthalpy of air at wall points

    qs = 1.1037E8/sqrt(Rmin) * sqrt(1.1 + 0.9*sqrt(Rmin/Rmax)) * \
        sqrt(rho_1/rho_sl) * (V_1/7925)**3.5 * (H_es - H_w) / (H_es - 3.0145E5)

    return qs


def Tauber(RN, Ma_1, rho_1, P_1, T_1, P_w, T_w):
    """
    Stagnation point heat flux estimation in zero angle of attack: Tauber function.

    Inputs:
    ------
    RN: curvature radius of stagnation point (act as placeholder)
    Ma_1: Mach number of free stream
    rho_1: density of free stream
    P_1: pressure of free stream
    T_1: temperature of free stream

    Output:
    ------
    qs: stagnation point heat flux

    Reference:
    [1] TAUBER M. 1989. A review of high-speed, convective, heat-transfer computation methods[R]. NASA Ames Research Center Moffett Field, CA, United States.
    [2] PARK S-H, NEEB D, PLYUSHCHEV G, et al. A study on heat flux predictions for re-entry flight analysis[J]. Acta Astronautica, 2021, 187: 271-280.
    """

    # intermediate parameters
    k = pred_k(T_1)
    V_1 = Ma_1 * sqrt(k * 287.05 * T_1)  # freestream velocity
    rho_2, P_2, T_2, rho_es, P_es, T_es = calShockThermo(rho_1, P_1, T_1, Ma_1)
    # sensible enthalpy of air at stagnation point
    hs = pred_Cp(T_1) * T_1 + V_1**2 / 2
    hw = pred_Cp(T_1) * T_w  # enthalpy of air at wall points

    qs = 1.83E-4 / sqrt(RN) * sqrt(rho_1) * V_1**3 * (1 - hw/hs)

    return qs


def LE_laminar(qsp, AOA, gamma):
    """
    Laminar heat flux modification for wing leading edge with AoA AOA and sweep angle gamma.

    Inputs:
    ------
    qsp: heat flux of sphere with the same radius as wing leading edge
    AOA: angle of attack (rad)
    gamma: sweep angle (rad, -pi/2 to pi/2)

    Output:
    ------
    qsl: laminar heat flux for wing leading edge
    """

    if abs(gamma) >= 0 and abs(gamma) <= pi/3:
        n = 1.5
    elif abs(gamma) > pi/3 and abs(gamma) < pi/2:
        n = 1.0
    else:
        raise RuntimeError("Sweep anagle is defined in (-pi/2, pi/2)")

    qsl = qsp * sqrt((1 - sin(gamma)**2 * cos(AOA)**2)**n / 2)

    return qsl


def LE_turbulence(qsp, AOA, gamma, RLE, Ma_1, rho_1, P_1, T_1):
    """
    Turbulence heat flux modification for wing leading edge with AoA AOA and sweep angle gamma.

    Inputs:
    ------
    qsp: heat flux of sphere with the same radius as wing leading edge
    AOA: angle of attack (rad)
    gamma: sweep angle (rad, -pi/2 to pi/2)
    RLE: radius of wing leading edge
    Ma_1: Mach number of free stream
    rho_1: density of free stream
    P_1: pressure of free stream
    T_1: temperature of free stream

    Output:
    ------
    qsl: turbulence heat flux for wing leading edge
    """

    # constants
    mu0 = 1.716E-5
    T0 = 273.15  # reference viscosity and temperature in Sutherland law
    Smu = 110.4  # Sutherland constant for standard air

    # intermediate parameters
    V_1 = Ma_1 * sqrt(pred_k(T_1) * 287.05 * T_1)
    mu1 = mu0 * (T_1/T0)**1.5 * (T0+Smu) / (T_1+Smu)  # Sutherland law
    gammae = arcsin(sin(gamma) * cos(AOA))

    qsl = qsp * 1.5/sqrt(2) * (2 * rho_1 * V_1 * RLE / mu1)**0.3 * \
        (0.01714 + 0.01235 * sin(3.53 * (gammae - pi/18)))

    return qsl


def wingLE(RN, Ma_1, rho_1, P_1, T_1, T_w, AOA, gamma):
    """
    #! deprecated
    Calculate wing leading edge heat flux on the analogy of yawing cylinder.

    References:
    [1] RESHOTKO E, BECKWITH I E. 1957. Compressible laminar boundary layer over a yawed infinite cylinder with heat transfer and arbitrary prandtl number[R]. Washington, USA: National Advisory Committee for Aeronautics.
    [2] 孟竹喧. 高超声速飞行器气动热环境仿真与防热结构分析[D]. 湖南长沙: 国防科学技术大学, 2014.
    """
    # constants
    g = 9.8067  # gravitational acceleration
    Pr = 0.71  # Prandtl constant
    Mbar = 28.9  # mass of 1 kmol air
    R0 = 848  # gas constant, kg·m/kmol

    # intermediate parameters
    k = pred_k(T_1)  # specific heat ratio
    gammae = arcsin(sin(gamma) * cos(AOA))  # effective sweeping angle
    Mae = Ma_1 * cos(gammae)
    rho_2, P_2, T_2, rho_es, P_es, T_es = calShockThermo(rho_1, P_1, T_1, Ma_1)
    # T10 = T_1 * (1 + (k-1)/2 * Ma_1**2) # freestream total temperature
    T10 = T_es
    # total temperature w.r.t. kinetic energy perpendicular to wing leading edge
    Tn0 = T_1 * (1 + (k-1)/2 * Mae**2)
    # Pwslg = P_1 * ((k+1)/2 * Mae**2) ** (k/(k-1)) * ((1+k) / (2*k*Mae - (k-1))) ** (k/(k-1))
    Pwslg = P_1 * ((k+1)/2 * Mae**2) ** (k/(k-1)) * \
        ((1+k) / (2*k*Mae**2 - (k-1))) ** (k/(k-1))
    # velocity gradient on stagnation line
    du_ds_sl = 1/RN * sqrt(2 * (Pwslg - P_1) / rho_1)
    # du_ds_sl = 1/RN * sqrt(2 * (P_w - P_1) / rho_1) # velocity gradient on stagnation line
    muwg = 1.73E-6 * (T_w/261)**(3/2) * (375/T_w + 114)  # sutherland law (?)
    rhowg = Pwslg * Mbar / (R0 * T_w * g)
    rhon0 = Pwslg * Mbar / (R0 * Tn0 * g)
    # rhowg = P_w * Mbar / (R0 * T_w  * g)
    # rhon0 = P_w * Mbar / (R0 * Tn0 * g)
    Tr = 1/2 * (T_w + T_es * cos(gammae)**2)  # reference temperature

    thetaw0 = 0.0014 * (T10/Tn0)**2.113 - 0.0109 * \
        (T10/Tn0)*1.113 + 0.516 * (T10/Tn0)**1.113
    thetaw = (1 + 1.5 * thetaw0**3.5 * T_w / T10) * thetaw0
    aw = thetaw * Pr**(-0.54) * g * sqrt(rhowg*muwg) * \
        pred_Cp(T_w) * sqrt(du_ds_sl)
    qsl = aw * (Tr - T_w)

    return qsl


def calShockThermo(rho_1, P_1, T_1, Ma_1):
    """
    Calculate thermodynamics parameters of wall and stagnation points.

    Inputs:
    ------
    rho: density of free stream
    P: pressure of free stream
    T: temperature of free stream
    Ma_1: freestream Mach number

    Outputs:
    -------
    rho_w, P_w, T_w: density, pressure, and temperature of wall points
    rho_es, P_es, T_es: density, pressure, and temperature of stagnation point
    """

    k = pred_k(T_1)  # specific heat ratio

    # thermodynamics parameters of wall and stagnation point
    rho_2 = (k+1) * Ma_1**2 / ((k-1) * Ma_1**2 + 2) * rho_1
    rho_es = (k+1) * Ma_1**2 / ((k-1) * Ma_1**2 + 2) * rho_1
    P_2 = (2*k * Ma_1**2 - (k-1)) / (k+1) * P_1
    P_es = ((k+1)**(k+1) * Ma_1**(2*k) /
            (2**(k+1) * k * Ma_1**2 - (k-1))) ** (1/(k-1)) * P_1
    T_2 = (2*k * Ma_1**2 - (k-1)) * \
        ((k-1) * Ma_1**2 + 2) / ((k+1)**2 * Ma_1**2) * T_1
    T_es = ((k-1) * Ma_1**2 + 2) / 2 * T_1

    return rho_2, P_2, T_2, rho_es, P_es, T_es


def pred_Cp(T):
    """
    Predict constant pressure specific heat capacity Cp (J/(kg·K)) w.r.t. temperature T (K).
    Perfect gas assumption (T < 2500K) adopted.
    Values at standard atmosphere (atm) (101325kPa) used.

    Reference:
    ---------
    [1] HILSENRATH J. Tables of thermal properties of gases: Comprising tables of thermodynamic and transport properties of air, argon, carbon dioxide, carbon monoxide, hydrogen, nitrogen, oxygen, and steam[M]. US Department of Commerce, National Bureau of Standards, 1955.
    """
    return 1.0436E3 - 0.3502*T + 9.0412E-4*T**2 - 5.8095E-7*T**3 + 1.2648E-10*T**4


def pred_k(T):
    """
    Predict specific heat ratio k w.r.t. temperature T (K).
    Perfect gas assumption (T < 2500K) adopted.
    Values at standard atmosphere (atm) (101325kPa) used.

    Reference:
    ---------
    [1] HILSENRATH J. Tables of thermal properties of gases: Comprising tables of thermodynamic and transport properties of air, argon, carbon dioxide, carbon monoxide, hydrogen, nitrogen, oxygen, and steam[M]. US Department of Commerce, National Bureau of Standards, 1955.
    """
    return 1.4028 + 8.5338E-5*T - 3.7143E-7*T**2 + 3.2068E-10*T**3 - 1.1551E-13*T**4 + 1.4685E-17*T**5


def calFricCalibra(Ma_1, T_1, P_1, ref_L):
    """
    Friction drag coefficient calibration for S/HABP.

    Inputs
    ------
    Ma: float, Mach number
    H: float, altitude (m)
    L: float, reference length for Reynolds number

    Output
    ------
    Cf: float, friction drag coefficient

    Reference
    ---------
    [1] MARINI M, PEZZELLA G, SCHETTINO A, et al. Numerical and experimental aerodynamic characterization of the hexafly-int hypersonic glider Proceedings of the 21st AIAA International Space Planes and Hypersonics Technologies Conference, F, 2017[C].
    """
    # constants
    mu0 = 1.716E-5
    T_sl = 273.15  # reference viscosity and temperature in Sutherland law
    Smu = 110.4  # Sutherland constant for standard air

    R = 287.0955
    rho_1 = P_1/R/T_1

    k = pred_k(T_1)
    a_1 = sqrt(R*k*T_1)
    V_1 = Ma_1*a_1

    # intermediate parameters
    rho_2, P_2, T_2, rho_es, P_es, T_es = calShockThermo(rho_1, P_1, T_1, Ma_1)
    H_w = pred_Cp(T_w)*T_w
    rho_w = airDensity(H_w, P_es)

    mu1 = mu0 * (T_1/T_sl)**1.5 * (T_sl+Smu) / (T_1+Smu)  # Sutherland law
    Rew = rho_w * V_1 * ref_L / mu1  # Reynolds number at wall

    Cf = 0.42 / (log(Rew)**2.55 * (1 + 0.25*Ma_1**2)**0.31)
    return Cf


def airDensity(H, P):
    # base pressure and enthalpy calculate entropy
    # H is kJ/kg
    P_1 = 101325
    if H <= 1755.52e3:  # 167.5kJ/kg<h<=1755.52kJ/kg
        a = 0.972
    else:  # 1755.52kJ/kg<h<35026kJ/kg
        a = 0.718+1.38974e-2*log(P/P_1)

    rho = 0.213833*(P/P_1)*(H/1755.52e3)**-a
    return rho


if __name__ == "__main__":
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from EvalAtmos import getAtmosphereEnvironment
    
    Z = 45000
    Ma = 13.8
    [T, P, rho, a, miu, g] = getAtmosphereEnvironment(Z)
    V = Ma*a
    T_w = 294

    print(calFricCalibra(Ma, T, P, 4))

    total_length = 4
    par_M_up = 0.75
    par_M_low = 0.75
    par_W = 2.8
    par_T = 0.7
    par_hight_up = 0.3
    par_hight_low = 0.3
    par_R = 0.013
    if par_M_up > 1:
        radius_center_up = 0
    else:
        dz_dx = (par_hight_up * (0.002) ** par_M_up) / \
            (total_length * 0.002)
        radius_center_up = par_R * dz_dx

    if par_M_low > 1:
        radius_center_low = 0
    else:
        dz_dx = (par_hight_low * (0.002) ** par_M_low) / \
            (total_length * 0.002)
        radius_center_low = par_R * dz_dx

    radius_center_head = (radius_center_up + radius_center_low) / 2
    radius_head_sq = (radius_center_head *
                      radius_center_head + par_R * par_R)
    radius_head_ZX = np.sqrt(radius_head_sq)
    # convert to y=f(x) can calculate x=0 point curvature
    radius_head_XY = 2**(1/par_T-2) * 2*par_W/(1/par_T/2)

    print(infStagHeat(radius_head_ZX, radius_head_ZX,
          Ma, P, T, T_w, None, 'Detra-Kemp-Riddell'))

    # # blunt
    # Ma_1 = 10.6
    # P_1 = 132
    # T_1 = 47.3397
    # radius_head = 0.0095
    # T_w=294
    # print(infStagHeat(radius_head, None, Ma_1, P_1, T_1,T_w, None, 'Detra-Kemp-Riddell'))
