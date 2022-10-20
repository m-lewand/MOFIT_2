from lib2to3.pgen2 import driver
import numpy as np
import matplotlib.pyplot as plt 
import Mesh
import math
import MES 
import Potentials
from scipy.linalg import eigh 



N = 2
L = 100 #[nm]
a = L/(2*N)
m = 0.067
a_b = 0.05292
dxi = 0.1
omega = 10. / 27211.6   #Hartree units
global_numbers = np.zeros((2*N+1, 2*N+1))
overlap_matrix = np.zeros((4,4))
kinetic_matrix = np.zeros((4,4))
S_matrix = np.zeros(((2*N+1)**2, (2*N+1)**2))
H_matrix = np.zeros(((2*N+1)**2, (2*N+1)**2))
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

#Global matrices merging 
for element in range(1, (2*N)**2 + 1):
    for i in range(1,5):
        for j in range(1,5):
            S_matrix[Mesh.find_global_number(element,i,N) - 1, Mesh.find_global_number(element,j,N) - 1] +=  overlap_matrix[i-1,j-1]
            H_matrix[Mesh.find_global_number(element,i,N) - 1, Mesh.find_global_number(element,j,N) - 1] += kinetic_matrix[i-1,j-1] + Potentials.Parabolic(element,i,j,m,omega,a/a_b,N)


#Boundary conditions
for i in range(1, (2*N+1)**2 + 1):
    #for j in range((2*N+1)**2):
    if i <= 2*N + 1:    #bottom row 
        S_matrix[i - 1,:] = 0 
        S_matrix[:,i - 1] = 0
        S_matrix[i - 1,i - 1] = 1

        H_matrix[i - 1,:] = 0 
        H_matrix[:,i - 1] = 0
        H_matrix[i - 1,i - 1] = -1410
    elif i <= (2*N+1)**2  and i >= (2*N+1)**2 - 2*N: #top row 
        S_matrix[i - 1,:] = 0 
        S_matrix[:,i - 1] = 0
        S_matrix[i - 1,i - 1] = 1

        H_matrix[i - 1,:] = 0 
        H_matrix[:,i - 1] = 0
        H_matrix[i - 1,i - 1] = -1410
    elif i % (2*N + 1) == 0:    #right column
        S_matrix[i - 1,:] = 0 
        S_matrix[:,i - 1] = 0
        S_matrix[i - 1,i - 1] = 1

        H_matrix[i - 1,:] = 0 
        H_matrix[:,i - 1] = 0
        H_matrix[i - 1,i - 1] = -1410
    elif i % (2*N + 1) == 1:    #left column
        S_matrix[i - 1,:] = 0 
        S_matrix[:,i - 1] = 0
        S_matrix[i - 1,i - 1] = 1

        H_matrix[i - 1,:] = 0 
        H_matrix[:,i - 1] = 0
        H_matrix[i - 1,i - 1] = -1410


#Solving for eigenvalues
Eigenvalues, Eigenvectors = eigh(H_matrix, S_matrix, type = 1, overwrite_a = True, overwrite_b = True)
print(Eigenvalues * 27211.6)


print(Eigenvectors[:,-1])
