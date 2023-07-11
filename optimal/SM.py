import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class RBF():
    '''
    radial basis fcn interp pre model fcn
    input initial data X,Y,which are real data
    X,Y are x_num x vari_num matrix
    aver_X,stdD_X is 1 x x_num matrix
    output is a radial basis model,include X,Y,base_fcn
    and predict_fcn

    Copyright 2023 Adel
    '''

    def __init__(self, X, Y, basis_fcn=lambda r: r**3):
        self.X = np.array(X)
        self.basis_fcn = basis_fcn

        self.x_num, self.vari_num = self.X.shape
        self.Y = np.array(Y).reshape(self.x_num, 1)

        self.aver_X = np.mean(X, axis=0).reshape((1, self.vari_num))
        self.stdD_X = np.std(X, ddof=1, axis=0).reshape((1, self.vari_num))
        self.stdD_X[self.stdD_X == 0] = 1

        self.aver_Y = np.mean(Y, axis=0).reshape((1, 1))
        self.stdD_Y = np.std(Y, ddof=1, axis=0).reshape((1, 1))
        self.stdD_Y[self.stdD_Y == 0] = 1
        
        pass

    def train(self):
        x_num = self.x_num
        vari_num = self.vari_num

        # normalize data
        self.X_nomlz = (self.X - self.aver_X) / self.stdD_X
        self.Y_nomlz = (self.Y - self.aver_Y) / self.stdD_Y

        # initialization distance of all X
        X_dis = np.zeros((x_num, x_num))
        for vari_idx in range(vari_num):
            X_dis = X_dis + (self.X_nomlz[:, vari_idx].reshape(x_num,1) -
                             np.transpose(self.X_nomlz[:, vari_idx]).reshape(1,x_num)) ** 2

        X_dis = np.sqrt(X_dis)

        self.rdibas_matrix = self.basis_fcn(X_dis)
        # stabilize matrix
        self.rdibas_matrix = self.rdibas_matrix + np.eye(self.x_num) * 1e-09
        # get inverse matrix
        self.inv_rdibas_matrix = np.linalg.solve(
            self.rdibas_matrix, np.eye(self.x_num))
        # solve beta
        self.beta = np.dot(self.inv_rdibas_matrix, self.Y_nomlz)
        
        return 0

    def predict(self, X_pred):
        X_pred = np.array(X_pred).reshape((-1,self.vari_num))

        x_pred_num, __ = X_pred.shape
        # normalize data
        X_pred_nomlz = (X_pred - self.aver_X) / self.stdD_X
        # calculate distance
        X_dis_pred = np.zeros((x_pred_num, self.x_num))
        for vari_idx in range(self.vari_num):
            X_dis_pred = X_dis_pred + \
                (X_pred_nomlz[:, vari_idx].reshape(x_pred_num, 1) -
                 np.transpose(self.X_nomlz[:, vari_idx]).reshape(1, self.x_num)) ** 2

        X_dis_pred = np.sqrt(X_dis_pred)
        # predict variance
        Y_pred_nomlz = np.dot(self.basis_fcn(X_dis_pred), self.beta)
        # normalize data
        Y_pred = Y_pred_nomlz * self.stdD_Y + self.aver_Y
        
        return Y_pred


class RBFMF():
    '''
    radial basis fcn interp pre model fcn
    input initial data X,Y,which are real data
    X,Y are x_num x vari_num matrix
    aver_X,stdD_X is 1 x x_num matrix
    output is a radial basis model,include X,Y,base_fcn
    and predict_fcn

    Copyright 2023 Adel
    '''

    def __init__(self, X_HF, Y_HF, basis_fcn_HF=lambda r: r**3, X_LF=None, Y_LF=None, basis_fcn_LF=lambda r: r**3):
        self.X_HF = np.array(X_HF)
        self.basis_fcn_HF = basis_fcn_HF
        self.x_HF_num, self.vari_num = self.X_HF.shape
        self.Y_HF = np.array(Y_HF).reshape(self.x_HF_num, 1)
        
        self.X_LF = np.array(X_LF)
        self.basis_fcn_LF = basis_fcn_LF
        self.Y_LF = np.array(Y_LF).reshape(self.X_LF.shape[0], 1)

        self.aver_X = np.mean(X_HF, axis=0).reshape((1, self.vari_num))
        self.stdD_X = np.std(X_HF, ddof=1, axis=0).reshape((1, self.vari_num))
        self.stdD_X[self.stdD_X == 0] = 1

        self.aver_Y = np.mean(Y_HF, axis=0).reshape((1, 1))
        self.stdD_Y = np.std(Y_HF, ddof=1, axis=0).reshape((1, 1))
        self.stdD_Y[self.stdD_Y == 0] = 1
        
        pass

    def train(self):
        # train LF
        LF_model=RBF(self.X_LF,self.Y_LF,self.basis_fcn_LF)
        LF_model.train()
        self.LF_model=LF_model
        
        x_HF_num = self.x_HF_num
        vari_num = self.vari_num

        # normalize data
        self.X_nomlz = (self.X_HF - self.aver_X) / self.stdD_X
        self.Y_nomlz = (self.Y_HF - self.aver_Y) / self.stdD_Y
        YHF_pred=LF_model.predict(self.X_HF)
        YHF_pred_nomlz = (YHF_pred - self.aver_Y) / self.stdD_Y
        
        # initialization distance of all X
        X_dis = np.zeros((x_HF_num, x_HF_num))
        for vari_idx in range(vari_num):
            X_dis = X_dis + (self.X_nomlz[:, vari_idx].reshape(x_HF_num,1) -
                             np.transpose(self.X_nomlz[:, vari_idx]).reshape(1,x_HF_num)) ** 2

        X_dis = np.sqrt(X_dis)

        self.rdibas_matrix = self.basis_fcn_HF(X_dis)
        # add low fildelity value
        self.H = np.concatenate((self.rdibas_matrix*YHF_pred_nomlz,self.rdibas_matrix),axis=1)
        self.H_hessian = np.dot(self.H , np.transpose(self.H))
        # stabilize matrix
        self.H_hessian = self.H_hessian + np.eye(x_HF_num) * 1e-09
        # get inv matrix
        self.inv_H_hessian = np.linalg.solve(self.H_hessian,np.eye(x_HF_num))
        # solve omega
        self.omega = np.dot(np.transpose(self.H) , np.dot(self.inv_H_hessian , self.Y_nomlz))
        self.alpha = self.omega[:x_HF_num]
        self.beta = self.omega[x_HF_num:]
        
        return 0

    def predict(self, X_pred):
        X_pred = np.array(X_pred).reshape((-1,self.vari_num))

        x_pred_num, __ = X_pred.shape
        # normalize data
        X_pred_nomlz = (X_pred - self.aver_X) / self.stdD_X
        # calculate distance
        X_dis_pred = np.zeros((x_pred_num, self.x_HF_num))
        for vari_idx in range(self.vari_num):
            X_dis_pred = X_dis_pred + \
                (X_pred_nomlz[:, vari_idx].reshape(x_pred_num, 1) -
                 np.transpose(self.X_nomlz[:, vari_idx]).reshape(1, self.x_HF_num)) ** 2

        X_dis_pred = np.sqrt(X_dis_pred)
        
        # predict low fildelity value
        Y_pred_LF = self.LF_model.predict(X_pred)
        # nomalizae
        Y_pred_LF_nomlz = (Y_pred_LF - self.aver_Y) / self.stdD_Y
        
        # combine two matrix
        rdibas_matrix_pred = self.basis_fcn_HF(X_dis_pred)
        H_pred = np.concatenate((rdibas_matrix_pred*Y_pred_LF_nomlz,rdibas_matrix_pred),axis=1)
        # predict data
        Y_pred_nomlz = np.dot(H_pred , self.omega)
        # normalize data
        Y_pred = Y_pred_nomlz * self.stdD_Y + self.aver_Y
        
        return Y_pred


def interpVisualize(model,low_bou,up_bou):
    draw_X, draw_Y = np.meshgrid(np.linspace(low_bou[0], up_bou[0], 21), np.linspace(low_bou[1], up_bou[1], 21))
    draw_Point = np.concatenate(
        (draw_X.reshape((441, 1)), draw_Y.reshape((441, 1))), axis=1)
    draw_Z = model.predict(draw_Point)
    draw_Z=draw_Z.reshape((21, 21))
    
    plt.contour(draw_X, draw_Y, draw_Z)
    plt.show()


if __name__ == '__main__':
    from scipy import io
    
    # data = io.loadmat('PK.mat')
    # low_bou = data['low_bou'].reshape(2,)
    # up_bou = data['up_bou'].reshape(2,)

    # surrogate = RBF(data['X'], data['Y'])
    # surrogate.train()

    # draw_X, draw_Y = np.meshgrid(np.linspace(
    #     low_bou[0], up_bou[0], 21), np.linspace(low_bou[1], up_bou[1], 21))
    # draw_Point = np.concatenate(
    #     (draw_X.reshape((441, 1)), draw_Y.reshape((441, 1))), axis=1)
    # draw_Z = surrogate.predict(draw_Point).reshape((21, 21))

    # fig = plt.figure()
    # axe = Axes3D(fig, auto_add_to_figure=False)
    # fig.add_axes(axe)
    # axe.plot_surface(draw_X, draw_Y, draw_Z)

    # axe.set_xlabel(r'$x$', fontsize=18)
    # axe.set_ylabel(r'$y$', fontsize=18)
    # axe.set_zlabel(r'$z$', fontsize=18)
    # axe.view_init(30, -120)
    # plt.axis('auto')
    # plt.show()
    
    
    data = io.loadmat('Forrester.mat')

    surrogate = RBFMF(data['XHF'], data['YHF'],X_LF=data['XLF'],Y_LF=data['YLF'])
    surrogate.train()

    draw_X=np.linspace(0,1,20).reshape(20,1)
    draw_Y = surrogate.predict(draw_X)

    plt.plot(draw_X,draw_Y)
    plt.plot(data['x'],data['y_real'],'-')
    plt.plot(data['x'],data['y_real_low'],'--')
    plt.axis('auto')
    plt.show()

