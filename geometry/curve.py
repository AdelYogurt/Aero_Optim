import numpy as np
from math import comb

class CurveBezier():

    def __init__(self, point_list, max_order=None, u_list=[]):
        # get Bezier

        # input:
        # node_list(Bezier will across), order, u_list(default is linspace(0,1,node_number))
        # control_list

        point_list=np.array(point_list)
        
        if max_order:
            # fitting Bezier by node point
            # input node_list, max_order
            # fitting method is least square method
            node_list = point_list
            node_number, dimension = node_list.shape
            
            if max_order > node_number - 1:
                raise Exception(
                    'CurveBezier: curve order large than node number-1')
                
            self.node_number = node_number
            self.node_list = node_list
            self.dimension = dimension
            self.order = max_order
            
            if len(u_list) > 0:
                u_list = np.array(u_list)
            else:
                u_list = np.linspace(0, 1, node_number)
                
            # generate fitting matrix
            FM = np.zeros((node_number, max_order + 1))
            for node_index in np.arange(1, node_number+1).reshape(-1):
                FM[node_index, :] = self.baseFunction(u_list(node_index))
            control_list = np.linalg.solve(FM, node_list)
            
            control_number = max_order + 1
            self.control_list = control_list
            self.control_number = control_number
        else:
            # generate Bezier by control point
            # input control_list
            control_list = point_list
            control_number, dimension = control_list.shape
            if control_number < 2:
                raise Exception('CurveBezier: control_number less than 2')
            order = control_number - 1
            self.control_list = control_list
            self.control_number = control_number
            self.order = order
            self.dimension = dimension    

        return

    def interpPoint(self, u_x_list):
        # according u_x calculate point
        # u_x_list is u_x_number x 1 matrix
        # point_list is point_number x dimension matrix

        u_x_list=np.array(u_x_list).reshape((-1,))
        u_x_number = len(u_x_list)
        point_list = np.zeros((u_x_number, self.dimension))
        for u_x_index in range(0, u_x_number):
            u_x = u_x_list[u_x_index]
            P_list = self.baseFunction(u_x)
            point_list[u_x_index, :] = np.dot(P_list,self.control_list)

        return point_list

    def baseFunction(self=None, u_x=None):
        max_order = self.order
        b_list = np.zeros((1, max_order + 1))
        x_list = np.zeros((1, max_order + 1))
        middle_index = int(np.ceil((max_order + 1) / 2))
        # calculate binom
        for b_index in range(0, middle_index):
            b_list[0,b_index] = comb(max_order, b_index)

        # symmetry
        b_list[0,middle_index:] = np.flip(b_list[0,:int(np.floor((max_order + 1) / 2))].reshape((1,-1)),axis=1)
        # calculate x
        for x_index in range(0, max_order+1):
            x_list[0,x_index] = u_x ** (x_index) * \
                (1 - u_x) ** (max_order - x_index)

        P_list = np.multiply(b_list, x_list)
        return P_list


def importdata(str_file):
    with open(str_file) as file_data:
        str_data = file_data.readlines()
    data = [data_unit.split()[0].split(',') for data_unit in str_data]
    return data


if __name__ == "__main__":

    curve_data = importdata('coord_data/airfoil_shape_up.txt')
    curve_data=[[float(data_unit) for data_unit in data] for data in curve_data]
    curve_data=np.array(curve_data)
    curve=CurveBezier(curve_data)
    print(curve.interpPoint(0.5))
    