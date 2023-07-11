import numpy as np
import copy
import math
import random
from scipy import io
from scipy.optimize import minimize,Bounds
import matplotlib.pyplot as plt
from scipy import special


class GPC():
    def __init__(self,X:np.array,Y,hyp = {'mean':[0],'cov':np.zeros(2)}):
        # generate gaussian process classifier model
        # version 6,this version is assembly of gpml-3.6 EP method
        # X is x_num x vari_num matirx,Y is x_num x 1 matrix
        # low_bou,up_bou is 1 x vari_num matrix
        # only support binary classification,-1 and 1

        # input:
        # X,Y,hyp(mean,cov(len,eta))

        # abbreviation:
        # pred: predicted,nomlz: normalization,num: num
        # var: variance
        self.X=np.array(X)
        self.Y=np.array(Y).reshape(len(Y),1)
        self.hyp=hyp

        self.x_num,self.vari_num = X.shape

        # normalization data
        self.aver_X = np.mean(X,axis=0).reshape((1,self.vari_num))
        self.stdD_X = np.std(X,ddof =1,axis=0).reshape((1,self.vari_num))
        self.stdD_X[self.stdD_X==0]=1
        
        self.inf_fcn=self.infEP
        self.mean_fcn=self.meanConst
        self.cov_fcn=self.calCov
        self.lik_fcn=self.likErf
        
        self.last_ttau=None
        self.last_tnu=None
        
        self.post=None
        
        pass
    
    
    def train(self): 
        self.X_nomlz = (self.X - self.aver_X) / self.stdD_X
        
        def objNLLGPC(x,inf_fcn): 
            hyp_iter={'mean':[x[0]],'cov':x[1:],'lik':[]}
            _,nlZ,_ = inf_fcn(hyp_iter)
            return nlZ
    
        def objNLLGPCGrad(x,inf_fcn): 
            hyp_iter={'mean':[x[0]],'cov':x[1:],'lik':[]}
            __,_,dnlZ = inf_fcn(hyp_iter)
            return np.append(dnlZ['mean'],dnlZ['cov']).reshape(3,)

        obj_fcn = lambda x: objNLLGPC(x,self.inf_fcn)
        obj_fcn_grad = lambda x: objNLLGPCGrad(x,self.inf_fcn)
        hyp_x = np.concatenate((self.hyp['mean'],self.hyp['cov']),axis=0)
        # [fval,gradient] = object_fcn(hyp_x)
        # [fval_differ,gradient_differ] = differ(object_fcn,hyp_x)

        hyp_x = minimize(obj_fcn, hyp_x, method='SLSQP', jac=obj_fcn_grad,bounds=Bounds(-3,3),
                         options={'ftol': 1e-3, 'disp': False})
        hyp_x=hyp_x.x
        
        # hyp.mean = hyp_x[0];
        self.hyp['mean'] = [0]
        self.hyp['cov'] = hyp_x[1:]
        self.hyp['lik'] = []
        self.post,nlZ,dnlZ = self.infEP(self.hyp)

        return 0

    def predict(self,X_pred): 
        X_pred = np.array(X_pred).reshape((-1,self.vari_num))
        # predict function
        X_pred_nomlz = (X_pred - self.aver_X) / self.stdD_X
        pred_num = X_pred_nomlz.shape[0]
        ys = np.ones((pred_num,1))
        alpha = self.post['alpha']
        L = self.post['L']
        sW = self.post['sW']

        #verify whether L contains valid Cholesky decomposition or something different
        # Lchol = isnumeric(L) and np.all(np.logical_and(np.all(tril(L,- 1) == 0),np.transpose(np.diag(L))) > np.logical_and(0,np.transpose(True)))
        x_pred_num = X_pred_nomlz.shape[0]

        nperbatch = 1000
        nact = 0

        ymu = np.zeros((x_pred_num,1))
        ys2 = ymu
        miu_pred = ymu
        var_pred = ymu
        possibility = ymu

        while nact < x_pred_num-1:
            id = range(nact,np.minimum(nact + nperbatch,x_pred_num))
            kss,dkss = self.cov_fcn(self.hyp['cov'],X_pred_nomlz[id,:],'diag')
            Ks,dKs = self.cov_fcn(self.hyp['cov'],self.X_nomlz,X_pred_nomlz[id,:])
            ms = self.mean_fcn(self.hyp['mean'],X_pred_nomlz[id,:])
            Fmu = ms + np.matmul(np.transpose(Ks),alpha)
            miu_pred[id] = Fmu
            # if Lchol:
            V = np.linalg.solve(np.transpose(L),sW*Ks)
            var_pred[id] = kss - np.transpose(np.sum(V*V, 0)).reshape((V.shape[1],1))
            # else:
            #     if isnumeric(L):
            #         LKs = L * Ks
            #     else:
            #         LKs = L(Ks)
            #     var_pred[id] = kss + np.transpose(np.sum(np.multiply(Ks,LKs), 0))
            var_pred[id] = np.amax(var_pred[id],0)
            Fs2 = var_pred[id]
            Ys = ys
            Lp,Ymu,Ys2 = self.lik_fcn(self.hyp['lik'],Ys,Fmu,Fs2)
            
            for idx in range(len(id)):
                possibility[id[idx]] = Lp[idx]
                ymu[id[idx]] = Ymu[idx]
                ys2[id[idx]] = Ys2[idx]
                
            nact = id[-1]

        possibility = np.exp(possibility)
        Y_pred = np.ones((pred_num,1))
        Y_pred[possibility < 0.5] = - 1
        return Y_pred,possibility,miu_pred,var_pred


    def calCov(self,cov,X,Z =None): 
        # obtain covariance of x
        # cov: eta,len(equal to 1/len.^2)
        x_num=self.x_num
        vari_num=self.vari_num
        
        # k = eta*exp(-sum(x_dis*len)/vari_num);
        len = np.exp(cov[0])
        eta = np.exp(cov[1])
        # predict
        if Z is not None:
            if str(Z) == str('diag'):
                K = eta
            else:
                z_num,_ = Z.shape
                # initializate square of X inner distance/ vari_num
                K = np.zeros((x_num,z_num))
                for len_idx in range(vari_num):
                    K = K + (X[:,len_idx].reshape(x_num,1) - np.transpose(Z[:,len_idx]).reshape(1,z_num)) ** 2 * len / vari_num
                K = eta * np.exp(- K)
                
            dK_dcov=[]
        else:
            # initializate square of X inner distance sq
            sq_dis_v = np.zeros((x_num,x_num,vari_num))
            for len_idx in range(vari_num):
                matrix=(X[:,len_idx].reshape(self.x_num,1) - np.transpose(X[:,len_idx]).reshape(1,self.x_num)) ** 2 / vari_num
                sq_dis_v[:,:,len_idx] = matrix
            # exp of x__x with theta
            exp_dis = np.zeros((x_num,x_num))
            for len_idx in range(vari_num):
                exp_dis = exp_dis + sq_dis_v[:,:,len_idx] * len
            exp_dis = np.exp(- exp_dis)
            K = exp_dis * eta
            
            dK_dcov = list()
            dK_dlen = np.zeros((x_num,x_num))
            for len_idx in range(vari_num):
                dK_dlen = dK_dlen + sq_dis_v[:,:,len_idx]
            dK_dlen = np.multiply(- dK_dlen,K) * len
            dK_dcov.append(dK_dlen) 
            dK_dcov.append(K)

        return K,dK_dcov


    def infEP(self,hyp): 
        # Expectation Propagation approximation to the posterior Gaussian Process.
        # The fcn takes a specified covariance fcn (see covfcns.m) and
        # likelihood fcn (see likfcns.m),and is designed to be used with
        # gp.m. See also infMethods.m. In the EP algorithm,the sites are
        # updated in random order,for better performance when cases are ordered
        # according to the targets.

        # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2013-09-13.

        # See also INFMETHODS.M.
        tol = 0.0001
        max_sweep = 10
        min_sweep = 2

        inf = 'infEP'
        n = self.X_nomlz.shape[0]
        K,dK = self.cov_fcn(hyp['cov'],self.X_nomlz)
        m = self.mean_fcn(hyp['mean'],self.X_nomlz)

        # A note on naming: variables are given short but descriptive names in
        # accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
        # and s2 are mean and variance,nu and tau are natural parameters. A leading t
        # means tilde,a subscript _ni means "not i" (for cavity parameters),or _n
        # for a vector of cavity parameters. N(f|mu,Sigma) is the posterior.

        # marginal likelihood for ttau = tnu = zeros(n,1); equals n*log(2) for likCum*
        lZ,_,_=self.lik_fcn(hyp['lik'],self.Y,m,np.diag(K).reshape(self.x_num,1))
        nlZ0 = - sum(lZ)
        if self.last_ttau is None:
            ttau = np.zeros((n,1))
            tnu = np.zeros((n,1))
            Sigma = K
            mu = m
            nlZ = nlZ0
        else:
            ttau = self.last_ttau
            tnu = self.last_tnu
            Sigma,mu,L,alpha,nlZ = self.epComputeParams(K,self.Y,ttau,tnu,hyp,m)
            if nlZ > nlZ0:
                ttau = np.zeros((n,1))
                tnu = np.zeros((n,1))
                Sigma = K
                mu = m
                nlZ = nlZ0

        nlZ_old = float('inf')
        sweep = 0

        while (np.abs(nlZ - nlZ_old) > tol and sweep < max_sweep) or sweep < min_sweep:

            nlZ_old = nlZ
            sweep = sweep + 1
            # for i in random.sample(range(self.x_num),self.x_num):
            for i in range(self.x_num):
                tau_ni = (1 / Sigma[i,i] - ttau[i,0]).reshape(1,1)
                nu_ni = (mu[i,0] / Sigma[i,i] - tnu[i,0]).reshape(1,1)
                # compute the desired derivatives of the indivdual log partition fcn
                lZ,dlZ,d2lZ = self.lik_fcn(hyp['lik'],[[self.Y[i,0]]],nu_ni / tau_ni,1 / tau_ni)
                ttau_old = ttau[i,0]
                tnu_old = tnu[i,0]
                ttau[i,0] = - d2lZ / (1 + d2lZ / tau_ni)
                ttau[i,0] = np.amax(ttau[i,0],0)
                tnu[i,0] = (dlZ - nu_ni / tau_ni * d2lZ) / (1 + d2lZ / tau_ni)
                dtt = ttau[i,0] - ttau_old
                dtn = tnu[i,0] - tnu_old
                si = Sigma[:,i]
                si=si.reshape(len(si),1)
                ci = dtt / (1 + dtt * si[i,0])
                Sigma = Sigma - ci * si * np.transpose(si).reshape(1,len(si))
                mu = mu - (ci * (mu[i,0] + si[i,0] * dtn) - dtn) * si
            # recompute since repeated rank-one updates can destroy numerical precision
            Sigma,mu,L,alpha,nlZ = self.epComputeParams(K,self.Y,ttau,tnu,hyp,m)

        self.last_ttau = ttau
        self.last_tnu = tnu

        post={'alpha':alpha,'sW':np.sqrt(ttau),'L':L}

        dnlZ = hyp
        tau_n = 1.0 / np.diag(Sigma).reshape(ttau.shape) - ttau
        nu_n = mu / np.diag(Sigma).reshape(tnu.shape) - tnu
        sW = np.sqrt(ttau)
        sW_matrix=sW*np.eye(self.x_num)
        F = alpha * (np.transpose(alpha)) - sW*(np.linalg.solve(L,(np.linalg.solve(np.transpose(L),sW_matrix))))
        K,dK = self.cov_fcn(hyp['cov'],self.X_nomlz)
        for i in range(len(hyp['cov'])):
            dnlZ['cov'][i] = - sum(sum(F*dK[i])) / 2
        for i in range(len(hyp['lik'])):
            _,dlik,d2lZ = self.lik_fcn(hyp['lik'],self.Y,nu_n / tau_n,1.0 / tau_n,i)
            dnlZ['lik'][i] = - sum(dlik)
        lZ,dlZ,d2lZ = self.lik_fcn(hyp['lik'],self.Y,nu_n / tau_n,1.0 / tau_n)
        for i in range(len(hyp['mean'])):
            dm = self.mean_fcn(hyp['mean'],self.X_nomlz,i)
            dnlZ['mean'][i] = - np.matmul(np.transpose(dlZ),dm)[0,0]

        return post,nlZ,dnlZ


    def epComputeParams(self,K,y,ttau,tnu,hyp,m): 
        # fcn to compute the parameters of the Gaussian approximation,Sigma and
        # mu,and the negative log marginal likelihood,nlZ,from the current site
        # parameters,ttau and tnu. Also returns L (useful for predictions).
        n = len(y)
        sW = np.sqrt(ttau)
        matrix=(sW.reshape(self.x_num,1) * np.transpose(sW).reshape(1,self.x_num))
        matrix=np.eye(n) +matrix *K
        L = np.transpose(np.linalg.cholesky(matrix))
        V = np.linalg.solve(np.transpose(L),(sW.reshape(self.x_num,1)*K))
        Sigma = K - np.matmul(np.transpose(V) , V)
        alpha = tnu - sW*(np.linalg.solve(L,(np.linalg.solve(np.transpose(L),((sW*(np.matmul(K,tnu) + m)))))))
        mu = np.matmul(K,alpha) + m
        v = np.diag(Sigma).reshape(self.x_num,1)
        tau_n = 1.0 / v - ttau
        nu_n = mu / np.diag(Sigma).reshape(self.x_num,1) - tnu
        lZ,dlZ,d2lZ = self.lik_fcn(hyp['lik'],y,nu_n / tau_n,1.0 / tau_n)
        p = tnu - np.multiply(m,ttau)
        q = nu_n - np.multiply(m,tau_n)

        nlZ = sum(np.log(np.diag(L))) - sum(lZ) - np.matmul(np.transpose(p),np.matmul(Sigma , p ))/ 2 + np.matmul(np.transpose(v) , p ** 2) / 2 - np.matmul(np.transpose(q) , (np.multiply((np.multiply(ttau / tau_n,q) - 2 * p),v))) / 2 - sum(np.log(1 + ttau / tau_n)) / 2
        return Sigma,mu,L,alpha,nlZ


    def meanConst(self,hyp,x,i=None): 
        # Constant mean fcn. The mean fcn is parameterized as:
        # m(x) = c
        # The hyperparameter is:
        # hyp = [ c ]
        # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2010-08-04.
        # See also MEANfcnS.M.
        c = hyp
        if i is None:
            A = c * np.ones((len(x),1))
        else:
            if i == 1:
                A = np.ones((len(x),1))
            else:
                A = np.zeros((len(x),1))
                
        return A


    def likErf(self,hyp,y=None,mu=None,s2=None,i=None): 
        # likErf - Error fcn or cumulative Gaussian likelihood fcn for binary
        # classification or probit regression. The expression for the likelihood is
        #   likErf(t) = (1+erf(t/sqrt(2)))/2 = normcdf(t).

        # Several modes are provided,for computing likelihoods,derivatives and moments
        # respectively,see likfcns.m for the details. In general,care is taken
        # to avoid numerical issues when the arguments are extreme.

        # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2014-03-19.

        # See also LIKfcnS.M.

        # if mu is None:
        #     varargout = np.array(['0'])
        #     return varargout

        if y is not None:
            y = np.sign(y)
            y[y == 0] = 1
        else:
            y = 1

        # if np.asarray(y).size == 0:
        #     y = 1

        if i is None:
            z = mu / np.sqrt(1 + s2)
            dlZ = np.array([])
            d2lZ = np.array([])
            z = z*y
            lZ,n_p = self.logphi(z)
            dlZ = (y*n_p) / np.sqrt(1 + s2)
            d2lZ = (- n_p*(z + n_p)) / (1 + s2)
            return lZ,dlZ,d2lZ
        else:
            return []


    def logphi(self,z): 
        # Safe computation of logphi(z) = log(normcdf(z)) and its derivatives
        #                    dlogphi(z) = normpdf(x)/normcdf(x).
        # The fcn is based on idx 5725 in Hart et al. and gsl_sf_log_erfc_e.

        # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2013-11-13.
        lp = np.zeros(z.shape)
        B1 = (z*z < 0.0492)
        id1 = [i for i in range(len(z)) if B1[i]]
        lp0 = - z[id1,0] / (np.sqrt(2 * np.pi))
        c = np.array([[0.00048204],[- 0.00142906],[0.0013200243174],[0.0009461589032],[- 0.0045563339802],[0.00556964649138],[0.00125993961762116],[- 0.01621575378835404],[0.02629651521057465],[- 0.001829764677455021],[2 * (1 - np.pi / 3)],[(4 - np.pi) / 3],[1],[1]])
        f = 0
        for i in range(14):
            f = lp0*(c[i] + f)
        lp[id1,0] = - 2 * f - np.log(2)
        B2 = (z < - 11.3137)
        id2 = [i for i in range(len(z)) if B2[i]]
        uid2= [i for i in range(len(z)) if not B2[i]]
        r = np.array([[1.2753666447299659],[5.019049726784267],[6.160209853109631],[7.409740605964742],[2.978865626393993]])
        q = np.array([[2.260528520767327],[9.396034016235054],[12.048951927855128],[17.081440747466004],[9.608965327192788],[3.3690752069827528]])
        num = 0.5641895835477551
        for i in range(5):
            num = (- z[id2,0]*num) / np.sqrt(2) + r[i]
        den = 1.0
        for i in range(6):
            den = (- z[id2,0]*den) / np.sqrt(2) + q[i]

        e = num / den
        lp[id2,0] = np.log(e / 2) - z[id2,0] ** 2 / 2
        B3 = np.logical_and(~B2 ,~ B1 )
        id3 = [i for i in range(len(z)) if B3[i]]
        lp[id3,0] = np.log(special.erfc(- z[id3,0] / np.sqrt(2))/ 2) 

        dlp = np.zeros(z.shape)
        dlp[id2,0] = np.abs(den / num) * np.sqrt(2 / np.pi)
        
        dlp[uid2,0 ] = np.exp((- z[uid2,0]*z[uid2,0]) / 2 - lp[uid2,0]) / np.sqrt(2 * np.pi)

        return lp,dlp


class GPCMF():
    def __init__(self,X_HF:np.array,Y_HF,X_LF:np.array,Y_LF,hyp = {'mean':[0],'cov':np.zeros(5)}):
        # generate gaussian process classifier model
        # version 6,this version is assembly of gpml-3.6 EP method
        # X is x_num x vari_num matirx,Y is x_num x 1 matrix
        # low_bou,up_bou is 1 x vari_num matrix
        # only support binary classification,-1 and 1

        # input:
        # X,Y,hyp(mean,cov(len,eta))

        # abbreviation:
        # pred: predicted,nomlz: normalization,num: num
        # var: variance
        self.X_HF=np.array(X_HF)
        self.Y_HF=np.array(Y_HF).reshape(len(Y_HF),1)
        self.X_LF=np.array(X_LF)
        self.Y_LF=np.array(Y_LF).reshape(len(Y_LF),1)
        self.X=np.concatenate((X_HF,X_LF),axis=0)
        self.Y=np.concatenate((Y_HF,Y_LF),axis=0)
        self.hyp=hyp

        self.x_HF_num,self.vari_num = X_HF.shape
        self.x_LF_num,self.vari_num = X_LF.shape
        self.x_num=self.x_HF_num+self.x_LF_num
        
        # normalization data
        self.aver_X = np.mean(self.X,axis=0).reshape((1,self.vari_num))
        self.stdD_X = np.std(self.X,ddof =1,axis=0).reshape((1,self.vari_num))
        self.stdD_X[self.stdD_X==0]=1
        
        self.inf_fcn=self.infEP
        self.mean_fcn=self.meanConst
        self.cov_fcn=self.calCovMF
        self.lik_fcn=self.likErf
        
        self.last_ttau=None
        self.last_tnu=None
        
        self.post=None
        
        pass
    
    
    def train(self): 
        self.X_nomlz = (self.X - self.aver_X) / self.stdD_X
        self.X_HF_nomlz = (self.X_HF - self.aver_X) / self.stdD_X
        self.X_LF_nomlz = (self.X_LF - self.aver_X) / self.stdD_X
        
        def objNLLGPC(x,inf_fcn): 
            hyp_iter={'mean':[x[0]],'cov':x[1:],'lik':[]}
            _,nlZ,_ = inf_fcn(hyp_iter)
            return nlZ
    
        def objNLLGPCGrad(x,inf_fcn): 
            hyp_iter={'mean':[x[0]],'cov':x[1:],'lik':[]}
            __,_,dnlZ = inf_fcn(hyp_iter)
            return np.append(dnlZ['mean'],dnlZ['cov']).reshape(6,)

        obj_fcn = lambda x: objNLLGPC(x,self.inf_fcn)
        obj_fcn_grad = lambda x: objNLLGPCGrad(x,self.inf_fcn)
        hyp_x = np.concatenate((self.hyp['mean'],self.hyp['cov']),axis=0)
        # [fval,gradient] = object_fcn(hyp_x)
        # [fval_differ,gradient_differ] = differ(object_fcn,hyp_x)

        hyp_x = minimize(obj_fcn, hyp_x, method='SLSQP', jac=obj_fcn_grad,bounds=Bounds(-3,3),
                         options={'disp': False})
        hyp_x=hyp_x.x
        
        # hyp.mean = hyp_x[0];
        self.hyp['mean'] = [0]
        self.hyp['cov'] = hyp_x[1:]
        self.hyp['lik'] = []
        self.post,nlZ,dnlZ = self.infEP(self.hyp)

        return 0

    def predict(self,X_pred): 
        # predict function
        X_pred = np.array(X_pred).reshape((-1,self.vari_num))
        X_pred_nomlz = (X_pred - self.aver_X) / self.stdD_X
        pred_num = X_pred_nomlz.shape[0]
        ys = np.ones((pred_num,1))
        alpha = self.post['alpha']
        L = self.post['L']
        sW = self.post['sW']

        #verify whether L contains valid Cholesky decomposition or something different
        # Lchol = isnumeric(L) and np.all(np.logical_and(np.all(tril(L,- 1) == 0),np.transpose(np.diag(L))) > np.logical_and(0,np.transpose(True)))
        x_pred_num = X_pred_nomlz.shape[0]

        nperbatch = 1000
        nact = 0

        ymu = np.zeros((x_pred_num,1))
        ys2 = ymu
        miu_pred = ymu
        var_pred = ymu
        possibility = ymu

        while nact < x_pred_num-1:
            id = range(nact,np.minimum(nact + nperbatch,x_pred_num))
            kss,dkss = self.cov_fcn(self.hyp['cov'],X_pred_nomlz[id,:],'diag')
            Ks,dKs = self.cov_fcn(self.hyp['cov'],self.X_nomlz,X_pred_nomlz[id,:])
            ms = self.mean_fcn(self.hyp['mean'],X_pred_nomlz[id,:])
            Fmu = ms + np.matmul(np.transpose(Ks),alpha)
            miu_pred[id] = Fmu
            # if Lchol:
            V = np.linalg.solve(np.transpose(L),sW*Ks)
            var_pred[id] = kss - np.transpose(np.sum(V*V, 0)).reshape((V.shape[1],1))
            # else:
            #     if isnumeric(L):
            #         LKs = L * Ks
            #     else:
            #         LKs = L(Ks)
            #     var_pred[id] = kss + np.transpose(np.sum(np.multiply(Ks,LKs), 0))
            var_pred[id] = np.amax(var_pred[id],0)
            Fs2 = var_pred[id]
            Ys = ys
            Lp,Ymu,Ys2 = self.lik_fcn(self.hyp['lik'],Ys,Fmu,Fs2)
            
            for idx in range(len(id)):
                possibility[id[idx]] = Lp[idx]
                ymu[id[idx]] = Ymu[idx]
                ys2[id[idx]] = Ys2[idx]
                
            nact = id[-1]

        possibility = np.exp(possibility)
        Y_pred = np.ones((pred_num,1))
        Y_pred[possibility < 0.5] = - 1
        return Y_pred,possibility,miu_pred,var_pred


    def calCovMF(self,cov,X,Z = None): 
        # obtain covariance of x
        # cov: lenD,etaD,lenL,etaL,rho
        # len equal to 1/len_origin.^2

        # # k = eta*exp(-sum(x_dis*theta)/vari_num);
        x_num=self.x_num
        x_HF_num=self.x_HF_num
        x_LF_num=self.x_LF_num
        vari_num=self.vari_num
        
        # k = eta*exp(-sum(x_dis*len)/vari_num);
        lenD = np.exp(cov[0])
        etaD = np.exp(cov[1])
        lenL = np.exp(cov[2])
        etaL = np.exp(cov[3])
        rho = np.exp(cov[4])
        
        # predict
        if Z is not None:
            if str(Z) == str('diag'):
                K = rho * rho * etaL + etaD
            else:
                z_num,_ = Z.shape
                # initializate square of X inner distance sq/ vari_num
                sq_dis_v = np.zeros((x_num,z_num,vari_num))
                for len_idx in range(vari_num):
                    sq_dis_v[:,:,len_idx] = (X[:,len_idx].reshape(-1,1) - np.transpose(Z[:,len_idx]).reshape(1,-1)) ** 2 / vari_num
                # exp of x__x with D
                exp_disD = np.zeros((x_HF_num,z_num))
                for len_idx in range(vari_num):
                    exp_disD = exp_disD + sq_dis_v[:x_HF_num,:,len_idx] * lenD
                exp_disD = np.exp(- exp_disD)
                # exp of x__x with L
                exp_disL = np.zeros((x_num,z_num))
                for len_idx in range(vari_num):
                    exp_disL = exp_disL + sq_dis_v[:,:,len_idx] * lenL
                exp_disL = np.exp(- exp_disL)
                # covariance
                K = exp_disL
                K[:x_HF_num,:] = rho * rho * etaL * K[:x_HF_num,:] + etaD * exp_disD
                K[x_HF_num:x_num,:] = rho * etaL * K[x_HF_num:x_num,:]
                
            dK_dcov=[]
        else:
            # initializate square of X inner distance sq/ vari_num
            sq_dis_v = np.zeros((x_num,x_num,vari_num))
            for len_idx in range(vari_num):
                sq_dis_v[:,:,len_idx] = (X[:,len_idx].reshape((-1,1)) - np.transpose(X[:,len_idx]).reshape((1,-1))) ** 2 / vari_num
                
            # exp of x__x with H
            exp_disD = np.zeros((x_num,x_num))
            for len_idx in range(vari_num):
                exp_disD[:x_HF_num,:x_HF_num] = exp_disD[:x_HF_num,:x_HF_num] + sq_dis_v[:x_HF_num,:x_HF_num,len_idx] * lenD
            exp_disD[:x_HF_num,:x_HF_num] = np.exp(- exp_disD[:x_HF_num,:x_HF_num])
            KD = etaD * exp_disD
            
            # exp of x__x with L
            exp_disL = np.zeros((x_num,x_num))
            for len_idx in range(vari_num):
                exp_disL = exp_disL + sq_dis_v[:,:,len_idx] * lenL
            exp_disL = np.exp(- exp_disL)
            eta_exp_disL = etaL * exp_disL
            
            # times rho: HH to rho2,HL to rho,LL to 1
            KL = eta_exp_disL
            KL[:x_HF_num,:x_HF_num] = (rho * rho) * eta_exp_disL[:x_HF_num,:x_HF_num]
            KL[:x_HF_num,x_HF_num:x_num] = rho * eta_exp_disL[:x_HF_num,x_HF_num:x_num]
            KL[x_HF_num:x_num,:x_HF_num] = np.transpose(KL[:x_HF_num,x_HF_num:x_num]).reshape((-1,x_HF_num))
            K = KL + KD
            
            dK_dcov = list()
            # len D
            dK_dlenD = np.zeros((x_num,x_num))
            for len_idx in range(vari_num):
                dK_dlenD[:x_HF_num,:x_HF_num] = dK_dlenD[:x_HF_num,:x_HF_num] + sq_dis_v[:x_HF_num,:x_HF_num,len_idx]
            dK_dlenD[:x_HF_num,:x_HF_num] = np.multiply(- dK_dlenD[:x_HF_num,:x_HF_num],KD[:x_HF_num,:x_HF_num]) * lenD
            dK_dcov.append(dK_dlenD)
            # eta D
            dK_dcov.append(KD)
            
            # len L
            dK_dlenL = np.zeros((x_num,x_num))
            for len_idx in range(vari_num):
                dK_dlenL = dK_dlenL + sq_dis_v[:,:,len_idx]
            dK_dlenL = -dK_dlenL*KL*lenL
            dK_dcov.append(dK_dlenL)
            
            # eta L
            dK_dcov.append(KL)
            
            # rho
            dK_drho = np.zeros((x_num,x_num))
            dK_drho[:x_HF_num,:x_HF_num] = 2 * rho * rho * eta_exp_disL[:x_HF_num,:x_HF_num]
            dK_drho[:x_HF_num,x_HF_num:x_num] = rho * eta_exp_disL[:x_HF_num,x_HF_num:x_num]
            dK_drho[x_HF_num:x_num,:x_HF_num] = np.transpose(dK_drho[:x_HF_num,x_HF_num:x_num])
            dK_dcov.append(dK_drho)

        return K,dK_dcov


    def infEP(self,hyp): 
        # Expectation Propagation approximation to the posterior Gaussian Process.
        # The fcn takes a specified covariance fcn (see covfcns.m) and
        # likelihood fcn (see likfcns.m),and is designed to be used with
        # gp.m. See also infMethods.m. In the EP algorithm,the sites are
        # updated in random order,for better performance when cases are ordered
        # according to the targets.

        # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2013-09-13.

        # See also INFMETHODS.M.
        tol = 0.0001
        max_sweep = 10
        min_sweep = 2

        inf = 'infEP'
        n = self.X_nomlz.shape[0]
        K,dK = self.cov_fcn(hyp['cov'],self.X_nomlz)
        m = self.mean_fcn(hyp['mean'],self.X_nomlz)

        # A note on naming: variables are given short but descriptive names in
        # accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
        # and s2 are mean and variance,nu and tau are natural parameters. A leading t
        # means tilde,a subscript _ni means "not i" (for cavity parameters),or _n
        # for a vector of cavity parameters. N(f|mu,Sigma) is the posterior.

        # marginal likelihood for ttau = tnu = zeros(n,1); equals n*log(2) for likCum*
        lZ,_,_=self.lik_fcn(hyp['lik'],self.Y,m,np.diag(K).reshape(self.x_num,1))
        nlZ0 = - sum(lZ)
        if self.last_ttau is None:
            ttau = np.zeros((n,1))
            tnu = np.zeros((n,1))
            Sigma = K
            mu = m
            nlZ = nlZ0
        else:
            ttau = self.last_ttau
            tnu = self.last_tnu
            Sigma,mu,L,alpha,nlZ = self.epComputeParams(K,self.Y,ttau,tnu,hyp,m)
            if nlZ > nlZ0:
                ttau = np.zeros((n,1))
                tnu = np.zeros((n,1))
                Sigma = K
                mu = m
                nlZ = nlZ0

        nlZ_old = float('inf')
        sweep = 0

        while (np.abs(nlZ - nlZ_old) > tol and sweep < max_sweep) or sweep < min_sweep:

            nlZ_old = nlZ
            sweep = sweep + 1
            # for i in random.sample(range(self.x_num),self.x_num):
            for i in range(self.x_num):
                tau_ni = (1 / Sigma[i,i] - ttau[i,0]).reshape(1,1)
                nu_ni = (mu[i,0] / Sigma[i,i] - tnu[i,0]).reshape(1,1)
                # compute the desired derivatives of the indivdual log partition fcn
                lZ,dlZ,d2lZ = self.lik_fcn(hyp['lik'],[[self.Y[i,0]]],nu_ni / tau_ni,1 / tau_ni)
                ttau_old = ttau[i,0]
                tnu_old = tnu[i,0]
                ttau[i,0] = - d2lZ / (1 + d2lZ / tau_ni)
                ttau[i,0] = np.amax(ttau[i,0],0)
                tnu[i,0] = (dlZ - nu_ni / tau_ni * d2lZ) / (1 + d2lZ / tau_ni)
                dtt = ttau[i,0] - ttau_old
                dtn = tnu[i,0] - tnu_old
                si = Sigma[:,i]
                si=si.reshape(len(si),1)
                ci = dtt / (1 + dtt * si[i,0])
                Sigma = Sigma - ci * si * np.transpose(si).reshape(1,len(si))
                mu = mu - (ci * (mu[i,0] + si[i,0] * dtn) - dtn) * si
            # recompute since repeated rank-one updates can destroy numerical precision
            Sigma,mu,L,alpha,nlZ = self.epComputeParams(K,self.Y,ttau,tnu,hyp,m)

        self.last_ttau = ttau
        self.last_tnu = tnu

        post={'alpha':alpha,'sW':np.sqrt(ttau),'L':L}

        dnlZ = hyp
        tau_n = 1.0 / np.diag(Sigma).reshape(ttau.shape) - ttau
        nu_n = mu / np.diag(Sigma).reshape(tnu.shape) - tnu
        sW = np.sqrt(ttau)
        sW_matrix=sW*np.eye(self.x_num)
        F = alpha * (np.transpose(alpha)) - sW*(np.linalg.solve(L,(np.linalg.solve(np.transpose(L),sW_matrix))))
        K,dK = self.cov_fcn(hyp['cov'],self.X_nomlz)
        for i in range(len(hyp['cov'])):
            dnlZ['cov'][i] = - sum(sum(F*dK[i])) / 2
        for i in range(len(hyp['lik'])):
            _,dlik,d2lZ = self.lik_fcn(hyp['lik'],self.Y,nu_n / tau_n,1.0 / tau_n,i)
            dnlZ['lik'][i] = - sum(dlik)
        lZ,dlZ,d2lZ = self.lik_fcn(hyp['lik'],self.Y,nu_n / tau_n,1.0 / tau_n)
        for i in range(len(hyp['mean'])):
            dm = self.mean_fcn(hyp['mean'],self.X_nomlz,i)
            dnlZ['mean'][i] = - np.matmul(np.transpose(dlZ),dm)[0,0]

        return post,nlZ,dnlZ


    def epComputeParams(self,K,y,ttau,tnu,hyp,m): 
        # fcn to compute the parameters of the Gaussian approximation,Sigma and
        # mu,and the negative log marginal likelihood,nlZ,from the current site
        # parameters,ttau and tnu. Also returns L (useful for predictions).
        n = len(y)
        sW = np.sqrt(ttau)
        matrix=(sW.reshape(self.x_num,1) * np.transpose(sW).reshape(1,self.x_num))
        matrix=np.eye(n) +matrix *K
        L = np.transpose(np.linalg.cholesky(matrix))
        V = np.linalg.solve(np.transpose(L),(sW.reshape(self.x_num,1)*K))
        Sigma = K - np.matmul(np.transpose(V) , V)
        alpha = tnu - sW*(np.linalg.solve(L,(np.linalg.solve(np.transpose(L),((sW*(np.matmul(K,tnu) + m)))))))
        mu = np.matmul(K,alpha) + m
        v = np.diag(Sigma).reshape(self.x_num,1)
        tau_n = 1.0 / v - ttau
        nu_n = mu / np.diag(Sigma).reshape(self.x_num,1) - tnu
        lZ,dlZ,d2lZ = self.lik_fcn(hyp['lik'],y,nu_n / tau_n,1.0 / tau_n)
        p = tnu - np.multiply(m,ttau)
        q = nu_n - np.multiply(m,tau_n)

        nlZ = sum(np.log(np.diag(L))) - sum(lZ) - np.matmul(np.transpose(p),np.matmul(Sigma , p ))/ 2 + np.matmul(np.transpose(v) , p ** 2) / 2 - np.matmul(np.transpose(q) , (np.multiply((np.multiply(ttau / tau_n,q) - 2 * p),v))) / 2 - sum(np.log(1 + ttau / tau_n)) / 2
        return Sigma,mu,L,alpha,nlZ


    def meanConst(self,hyp,x,i=None): 
        # Constant mean fcn. The mean fcn is parameterized as:
        # m(x) = c
        # The hyperparameter is:
        # hyp = [ c ]
        # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2010-08-04.
        # See also MEANfcnS.M.
        c = hyp
        if i is None:
            A = c * np.ones((len(x),1))
        else:
            if i == 1:
                A = np.ones((len(x),1))
            else:
                A = np.zeros((len(x),1))
                
        return A


    def likErf(self,hyp,y=None,mu=None,s2=None,i=None): 
        # likErf - Error fcn or cumulative Gaussian likelihood fcn for binary
        # classification or probit regression. The expression for the likelihood is
        #   likErf(t) = (1+erf(t/sqrt(2)))/2 = normcdf(t).

        # Several modes are provided,for computing likelihoods,derivatives and moments
        # respectively,see likfcns.m for the details. In general,care is taken
        # to avoid numerical issues when the arguments are extreme.

        # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2014-03-19.

        # See also LIKfcnS.M.

        # if mu is None:
        #     varargout = np.array(['0'])
        #     return varargout

        if y is not None:
            y = np.sign(y)
            y[y == 0] = 1
        else:
            y = 1

        # if np.asarray(y).size == 0:
        #     y = 1

        if i is None:
            z = mu / np.sqrt(1 + s2)
            dlZ = np.array([])
            d2lZ = np.array([])
            z = z*y
            lZ,n_p = self.logphi(z)
            dlZ = (y*n_p) / np.sqrt(1 + s2)
            d2lZ = (- n_p*(z + n_p)) / (1 + s2)
            return lZ,dlZ,d2lZ
        else:
            return []


    def logphi(self,z): 
        # Safe computation of logphi(z) = log(normcdf(z)) and its derivatives
        #                    dlogphi(z) = normpdf(x)/normcdf(x).
        # The fcn is based on idx 5725 in Hart et al. and gsl_sf_log_erfc_e.

        # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2013-11-13.
        lp = np.zeros(z.shape)
        B1 = (z*z < 0.0492)
        id1 = [i for i in range(len(z)) if B1[i]]
        lp0 = - z[id1,0] / (np.sqrt(2 * np.pi))
        c = np.array([[0.00048204],[- 0.00142906],[0.0013200243174],[0.0009461589032],[- 0.0045563339802],[0.00556964649138],[0.00125993961762116],[- 0.01621575378835404],[0.02629651521057465],[- 0.001829764677455021],[2 * (1 - np.pi / 3)],[(4 - np.pi) / 3],[1],[1]])
        f = 0
        for i in range(14):
            f = lp0*(c[i] + f)
        lp[id1,0] = - 2 * f - np.log(2)
        B2 = (z < - 11.3137)
        id2 = [i for i in range(len(z)) if B2[i]]
        uid2= [i for i in range(len(z)) if not B2[i]]
        r = np.array([[1.2753666447299659],[5.019049726784267],[6.160209853109631],[7.409740605964742],[2.978865626393993]])
        q = np.array([[2.260528520767327],[9.396034016235054],[12.048951927855128],[17.081440747466004],[9.608965327192788],[3.3690752069827528]])
        num = 0.5641895835477551
        for i in range(5):
            num = (- z[id2,0]*num) / np.sqrt(2) + r[i]
        den = 1.0
        for i in range(6):
            den = (- z[id2,0]*den) / np.sqrt(2) + q[i]

        e = num / den
        lp[id2,0] = np.log(e / 2) - z[id2,0] ** 2 / 2
        B3 = np.logical_and(~B2 ,~ B1 )
        id3 = [i for i in range(len(z)) if B3[i]]
        lp[id3,0] = np.log(special.erfc(- z[id3,0] / np.sqrt(2))/ 2) 

        dlp = np.zeros(z.shape)
        dlp[id2,0] = np.abs(den / num) * np.sqrt(2 / np.pi)
        
        dlp[uid2,0 ] = np.exp((- z[uid2,0]*z[uid2,0]) / 2 - lp[uid2,0]) / np.sqrt(2 * np.pi)

        return lp,dlp


if __name__ == '__main__':
    data=io.loadmat('C_30.mat')
    up_bou = data['up_bou'].reshape(2,)
    low_bou = data['low_bou'].reshape(2,)
    
    model_GPC=GPC(data['X'],data['Y'])
    model_GPC.train()
    
    draw_X, draw_Y = np.meshgrid(np.linspace(
        low_bou[0], up_bou[0], 21), np.linspace(low_bou[1], up_bou[1], 21))
    draw_Point = np.concatenate(
        (draw_X.reshape((441, 1)), draw_Y.reshape((441, 1))), axis=1)
    draw_Z,possibility,miu_pred,var_pred = model_GPC.predict(draw_Point)
    possibility=possibility.reshape((21, 21))
    
    plt.contour(draw_X, draw_Y, possibility)
    plt.show()
    
    # data=io.loadmat('CMF_34.mat')
    # low_bou = data['low_bou'].reshape(2,)
    # up_bou = data['up_bou'].reshape(2,)
    
    # X_HF=data['XHF']
    # Y_HF=data['YHF']
    # X_LF=data['XLF']
    # Y_LF=data['YLF']
    
    # model_GPCMF=GPCMF(data['XHF'],data['YHF'],data['XLF'],data['YLF'])
    # model_GPCMF.train()
    
    # draw_X, draw_Y = np.meshgrid(np.linspace(
    #     low_bou[0], up_bou[0], 21), np.linspace(low_bou[1], up_bou[1], 21))
    # draw_Point = np.concatenate(
    #     (draw_X.reshape((441, 1)), draw_Y.reshape((441, 1))), axis=1)
    # draw_Z,possibility,miu_pred,var_pred = model_GPCMF.predict(draw_Point)
    # possibility=possibility.reshape((21, 21))
    
    # plt.contour(draw_X, draw_Y, possibility)
    # plt.show()
    
    