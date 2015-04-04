"""
post process wind turbine Cp data
"""


import csv
import math
import matplotlib.pyplot as plt

def normalize_x(x):
    
    
    minimum = min(x)
    
    x = [i+abs(minimum) for i in x]
    
    maximum = max(x)
   
    for i in range(0,len(x),1):
        x[i] = (x[i]) / maximum
        
    return x


def main():
  
    with open('pressure_data.csv','rb') as csvfile:
        file_data = csv.reader(csvfile)
        data = list(file_data)
        data = data[1:]
    
    x  = []
    y  = []
    z  = []
    cp = []
    
    for i in range(0,len(data),1):
        x.append(data[i][0])
        y.append(data[i][1])
        z.append(data[i][2])
        cp.append(data[i][3])
    
    x_30 = []
    cp_30 = []
    
    for i in range(0,len(z),1):
        
        if float(z[i]) > 1.45 and float(z[i]) < 1.6 :
            x_30.append(float(x[i])/0.711)
            cp_30.append(float(cp[i]))
    
    x_30 = normalize_x(x_30)
    
    # experimental data
    xu_30  = [ 1, 0.8, 0.68, 0.56, 0.36, 0.2, 0.1, 0.06, 0.04, 0.02, 0.01, 0.005, 0, 0, 0.005, 0.01, 0.02, 0.06, 0.14, 0.28, 0.44, 0.68, 0.92, 1]
    cpu_30 = [ 0.167175, -0.006089, -0.16294, -0.466726, -0.969941, -1.134659, -1.220121, -1.430707, -1.626495, -1.815081, -1.780514, -2.114277, -2.752586, -2.752586, 0.576019, 0.918158, 1, 0.855112, 0.506722, -0.008466, -0.223934, 0.172738, 0.292909, 0.167175]
    
    
    plt.figure()
    plt.plot(x_30,cp_30,'r-',marker='v',label='panel method')
    plt.plot(xu_30,cpu_30,'bo',label='NREL Expt')
    plt.ylim([-3,2])
    plt.legend(loc=1)
    plt.xlabel('x/c')
    plt.ylabel('Cp')
    plt.title('30% span, U = 7 m/s')
    plt.gca().invert_yaxis()
    #plt.show()
    plt.savefig('cp_7ms_30.eps')
    plt.close()
    
    
    x_47 = []
    cp_47 = []
    
    for i in range(0,len(z),1):
        
        if float(z[i]) > 2.2 and float(z[i]) < 2.4 :
            x_47.append(float(x[i])/0.627)
            cp_47.append(float(cp[i]))
    
    x_47 = normalize_x(x_47)
    
    # experimental data
    xu_47 = [0, 0.005, 0.01, 0.02, 0.06, 0.14, 0.28, 0.44, 0.68, 0.92, 1,
             1, 0.8, 0.68, 0.56, 0.36, 0.2, 0.1, 0.06, 0.04, 0.02, 0.01, 0.005, 0]
    cpu_47 = [-2.575825, 0.332103, 0.829754, 1, 0.852922, 0.496046, -0.019364, -0.195252, 0.179667, 0.307379, 0.102799,
              0.102799, -0.047678, -0.204188, -0.466682, -1.053567, -1.288154, -1.469653, -1.655028, -1.843953, -2.274761, -2.467014, -2.321913, -2.575825]
    
    plt.figure()
    plt.plot(x_47,cp_47,'r-',marker='v',label='panel method')
    plt.plot(xu_47,cpu_47,'bo',label='NREL Expt')
    plt.ylim([-3,2])
    plt.gca().invert_yaxis()
    plt.legend(loc=1)
    plt.xlabel('x/c')
    plt.ylabel('Cp')
    plt.title('47% span, U = 7 m/s')
    #plt.show()
    plt.savefig('cp_7ms_47.eps')
    plt.close()
    
    
    x_63 = []
    cp_63 = []
    
    for i in range(0,len(z),1):
        
        if float(z[i]) > 3.2 and float(z[i]) < 3.25 :
            x_63.append(float(x[i])/0.543)
            cp_63.append(float(cp[i]))
    
    x_63 = normalize_x(x_63)
    
    # experimental data
    xu_63 = [0, 0.005, 0.01, 0.02, 0.06, 0.14, 0.28, 0.44, 0.68, 0.92, 1,
             1, 0.8, 0.68, 0.56, 0.36, 0.2, 0.1, 0.06, 0.04, 0.02, 0.01, 0.005, 0]
    cpu_63 = [-1.919118, 0.663699, 0.953923, 1, 0.792005, 0.402968, -0.091503, -0.278885, 0.182597, 0.352175, 0.174389,
              0.174389, -0.058211, -0.259285, -0.605071, -1.068455, -1.229872, -1.365889, -1.512394, -1.724331, -1.915798, -1.995034, -1.840502, -1.919118]
    
    
    plt.figure()
    plt.plot(x_63,cp_63,'r-',marker='v',label='panel method')
    plt.plot(xu_63,cpu_63,'bo',label='NREL Expt')
    plt.ylim([-3,2])
    plt.gca().invert_yaxis()
    plt.legend(loc=1)
    plt.xlabel('x/c')
    plt.ylabel('Cp')
    plt.title('63% span, U = 7 m/s')
    #plt.show()
    plt.savefig('cp_7ms_63.eps')
    plt.close()
    
    
    x_80 = []
    cp_80 = []
    
    for i in range(0,len(z),1):
        
        if float(z[i]) > 3.8 and float(z[i]) < 4.1 :
            x_80.append(float(x[i])/0.457)
            cp_80.append(float(cp[i]))
    
    x_80 = normalize_x(x_80)
    
    # experimental data
    xu_80 = [ 0, 0.005, 0.01, 0.02, 0.06, 0.14, 0.28, 0.44, 0.68, 0.92, 1,
             1, 0.8, 0.68, 0.56, 0.36, 0.2, 0.1, 0.06, 0.04, 0.02, 0.01, 0.005, 0]
    cpu_80 = [-1.460808, 0.814826, 1, 0.981287, 0.724081, 0.3118, -0.160666, -0.348032, 0.126752, 0.346738, 0.211074,
              0.211074, -0.068508, -0.255429, -0.596876, -1.079419, -1.13962, -1.241882, -1.286436, -1.385176, -1.591168, -1.71867, -1.55582, -1.460808]
    
    plt.figure()
    plt.plot(x_80,cp_80,'r-',marker='v',label='panel method')
    plt.plot(xu_80,cpu_80,'bo',label='NREL Expt')
    plt.ylim([-3,2])
    plt.gca().invert_yaxis()
    plt.legend(loc=1)
    plt.xlabel('x/c')
    plt.ylabel('Cp')
    plt.title('80% span, U = 7 m/s')
    #plt.show()
    plt.savefig('cp_7ms_80.eps')
    plt.close()
    
    
    x_95 = []
    cp_95 = []
    
    for i in range(0,len(z),1):
        
        if float(z[i]) > 4.9 :
            x_95.append(float(x[i])/0.358)
            cp_95.append(float(cp[i]))
    
    x_95 = normalize_x(x_95)
    
    # experimental data
    xu_95 = [0, 0.005, 0.01, 0.02, 0.06, 0.14, 0.28, 0.44, 0.68, 0.92, 1,
             1, 0.8, 0.68, 0.56, 0.36, 0.2, 0.1, 0.06, 0.04, 0.02, 0.01, 0.005, 0]
    
    cpu_95 = [-0.469789, 0.916363, 1, 0.906636, 0.591274, 0.197514, -0.342239, -0.479616, 0.107823, 0.312772, 0.22502,
              0.22502, -0.031503, -0.184658, -0.390161, -0.838815, -0.874251, -0.942373, -0.988781, -1.066789, -1.132228, -0.998347, -0.722063, -0.469789]
    
    plt.figure()
    plt.plot(x_95,cp_95,'r-',marker='v',label='panel method')
    plt.plot(xu_95,cpu_95,'bo',label='NREL Expt')
    plt.ylim([-3,2])
    plt.gca().invert_yaxis()
    plt.legend(loc=1)
    plt.xlabel('x/c')
    plt.ylabel('Cp')
    plt.title('95% span, U = 7 m/s')
    #plt.show()
    plt.savefig('cp_7ms_95.eps')
    plt.close()
    
    

if __name__ == "__main__":
    main()
