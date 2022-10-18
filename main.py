import numpy as np
import matplotlib.pyplot as plt 
import Mesh
import math
import MES 

N = 10
L = 100 #[nm]
a = L/(2*N)
m = 0.067
a_b = 0.05292
dxi = 0.1
omega = 10. / 27211.6   #Hartree units
global_numbers = np.zeros((2*N+1, 2*N+1))
overlap_matrix = np.zeros((4,4))
kinetic_matrix = np.zeros((4,4))
global_numbers = Mesh.create_global_numbers(N,global_numbers)
for i in range(1,4*N**2 + 1):
    for j in range(1,5):
        #print((i,j,Mesh.find_global_number(i,j,N)))
        number = Mesh.find_global_number(i,j,N)
        #print(find_real_space_coordinate(number,N,a))
        

PSI = MES.MES(dxi,N,a, m, omega, a_b)
MES.local_overlap_matrix(overlap_matrix,a/a_b)  #length scaled to Hartree units, so the overlap matrix is in Hartee units
MES.local_kinetic_energy_matrix(kinetic_matrix, m)  #in Hartree units
print(overlap_matrix*9*4/(a/a_b)**2)
print(kinetic_matrix*2*m*6)

text_file = open("data.txt", "w")
for i in range((2*N)*int(2/dxi)) :
    for j in range((2*N)*int(2/dxi)) :
        print(PSI[i,j,0], PSI[i,j,1], PSI[i,j,2], file=text_file)
text_file.close()