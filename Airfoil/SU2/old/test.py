import sys
import csv
import math

'''
add air center to head distance to all vsc
'''
def addAirCenter(flow_velocity_list=[],AOA_list=[]):
    Cl_index=9
    Mz_index=13
    end_index=17
    
    for flow_velocity in flow_velocity_list:
        for AOA in AOA_list:
            filename_csv = 'tur_inc_history_'+str(AOA)+'_'+str(flow_velocity)+'.csv'
            
            # load data
            csv_data = open(filename_csv, mode='r', newline='', errors='ignore')
            data_list = list(csv.reader(csv_data))
            csv_data.close()
            
            # write data
            csv_data = open(filename_csv, mode='w', newline='', errors='ignore')
            csv_writer = csv.writer(csv_data)
    
            data_row = data_list[0]
            csv_writer.writerow(data_row[0:end_index+1]+['air_center'])

            for data_row in data_list[1:]:
                csv_writer.writerow(data_row[0:end_index+1] + [str(0-float(data_row[Mz_index])/float(data_row[Cl_index]))])
                
            csv_data.close()

'''
calculate all CL CD CM to centroid
average last N time
'''
def calAerodynamic(filename_output,centroid_distance,average_time=10,flow_velocity_list=[],AOA_list=[]):
    Cd_index=8
    Cl_index=9
    Mz_index=13
    air_center_index=18
    
    Cd_matrix=[]
    Cl_matrix=[]
    My_matrix=[]
    air_center_matrix=[]
    
    for flow_velocity in flow_velocity_list:
        Cd_list=[]
        Cl_list=[]
        My_list=[]
        air_center_list=[]
        
        for AOA in AOA_list:
            # load data
            filename_csv = 'tur_inc_history_'+str(AOA)+'_'+str(flow_velocity)+'.csv'
            csv_data = open(filename_csv, mode='r', newline='', errors='ignore')
            data_list = list(csv.reader(csv_data))
            csv_data.close()
            
            # average data
            Cd_average=0
            Cl_average=0
            My_average=0
            air_center_average=0
            for average_index in range(0,average_time):
                Cd_average=Cd_average+(float(data_list[-1-average_index][Cd_index]))
                Cl_average=Cl_average+(float(data_list[-1-average_time][Cl_index]))
                My_average=My_average+(float(data_list[-1-average_time][Mz_index]))
                air_center_average=air_center_average+(float(data_list[-1-average_time][air_center_index]))
                
            Cd_average=Cd_average/average_time
            Cl_average=Cl_average/average_time
            My_average=My_average/average_time
            air_center_average=air_center_average/average_time
            
            # ouput data
            Cd_list.append(Cd_average)
            Cl_list.append(Cl_average)
            My_list.append(My_average)
            air_center_list.append(air_center_average)
            
        Cd_matrix.append(Cd_list)
        Cl_matrix.append(Cl_list)
        My_matrix.append(My_list)
        air_center_matrix.append(air_center_list)
    
    file_output=open(filename_output,mode='w',newline='', errors='ignore')
    
    # output Cd
    file_output.write('Cd,0,1,2,3,5,10\n')
    for flow_velocity_index in range(0,len(flow_velocity_list)):
        file_output.write(str(flow_velocity_list[flow_velocity_index]))
        for AOA_index in range(0,len(AOA_list)):
            file_output.write(','+str(Cd_matrix[flow_velocity_index][AOA_index]))
        file_output.write('\n')
    file_output.write('\n')
        
    # output Cl
    file_output.write('Cl,0,1,2,3,5,10\n')
    for flow_velocity_index in range(0,len(flow_velocity_list)):
        file_output.write(str(flow_velocity_list[flow_velocity_index]))
        for AOA_index in range(0,len(AOA_list)):
            file_output.write(','+str(Cl_matrix[flow_velocity_index][AOA_index]))
        file_output.write('\n')
    file_output.write('\n')
        
    # output My
    file_output.write('My,0,1,2,3,5,10\n')
    for flow_velocity_index in range(0,len(flow_velocity_list)):
        file_output.write(str(flow_velocity_list[flow_velocity_index]))
        for AOA_index in range(0,len(AOA_list)):
            file_output.write(','+str(My_matrix[flow_velocity_index][AOA_index]))
        file_output.write('\n')
    file_output.write('\n')
        
    # output air_center
    file_output.write('air_center,0,1,2,3,5,10\n')
    for flow_velocity_index in range(0,len(flow_velocity_list)):
        file_output.write(str(flow_velocity_list[flow_velocity_index]))
        for AOA_index in range(0,len(AOA_list)):
            file_output.write(','+str(air_center_matrix[flow_velocity_index][AOA_index]))
        file_output.write('\n')
    file_output.write('\n')
    
    file_output.close()


if __name__ == '__main__':
    
    flow_velocity_list = [50]
    AOA_list = [0,1,2,3,5,10]
    
    # flow_velocity_list=[10]
    # AOA_list=[0]
    
    filename_output='aerodynamic_data.csv'
    centroid_distance=0.0
    
    addAirCenter(flow_velocity_list=flow_velocity_list,AOA_list=AOA_list)
    calAerodynamic(filename_output,centroid_distance,flow_velocity_list=flow_velocity_list,AOA_list=AOA_list)
                
