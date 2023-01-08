from lib2to3.pgen2 import driver
import numpy as np
import matplotlib.pyplot as plt 
import Mesh
import math
import MES 
import Potentials
from scipy.linalg import eigh 



N = 10
L = 80 #[nm]
a = L/(2*N)
m = 0.067
a_b = 0.05292
dxi = 0.1
omega = 10. / 27211.6   #Hartree units
#for gaussian confinement potential
sigma = 5 # [nm]
gaussian_deg = 1
gaussian_depth = -50 / 27211.6    #Hartree units 

time_evolution = False

global_numbers = np.zeros((2*N+1, 2*N+1))
overlap_matrix = np.zeros((4,4))
kinetic_matrix = np.zeros((4,4))
S_matrix = np.zeros(((2*N+1)**2, (2*N+1)**2))
H_matrix = np.zeros(((2*N+1)**2, (2*N+1)**2))
global_numbers = Mesh.create_global_numbers(N,global_numbers)

text_file = open("Z1.txt", "w")
print('Numer elementu / Lokalny numer wezla / Globalny numer wezla / (x,y)', file = text_file)
for i in range(1,4*N**2 + 1):
    for j in range(1,5):
        number = Mesh.find_global_number(i,j,N)
        print((i,j,Mesh.find_global_number(i,j,N), Mesh.find_real_space_coordinate(number,N,a)), file = text_file)
        #number = Mesh.find_global_number(i,j,N)
        #print(find_real_space_coordinate(number,N,a))
text_file.close()  

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
            #H_matrix[Mesh.find_global_number(element,i,N) - 1, Mesh.find_global_number(element,j,N) - 1] += kinetic_matrix[i-1,j-1] + Potentials.Parabolic(element,i,j,m,omega,a/a_b,N)
            H_matrix[Mesh.find_global_number(element,i,N) - 1, Mesh.find_global_number(element,j,N) - 1] += kinetic_matrix[i-1,j-1] + Potentials.generalized_gaussian(element,i,j,m,gaussian_depth,a/a_b,N, sigma/a_b, gaussian_deg)

            if element == 11 and i == j:
                #print(Potentials.Parabolic(element,i,j,m,omega,a/a_b,N)/overlap_matrix[i-1,j-1] * 27211.6, end = ' ')
                print('Diagonal element',i, i,' = ', Potentials.Parabolic(element,i,j,m,omega,a/a_b,N)/overlap_matrix[i-1,j-1] * 27211.6)
                print('Diagonal gaussian element',i, i,' = ', Potentials.generalized_gaussian(element,i,j,m,gaussian_depth,a/a_b,N, sigma/a_b, gaussian_deg)/overlap_matrix[i-1,j-1] * 27211.6)

                if j ==4:
                    print('\n')



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

printed = 1 
indexes = []
text_file = open("Eigenvalues.txt", "w")
for i in range(len(Eigenvalues)):
    if Eigenvalues[i] > 0 and printed <= 15:
        print('Eigenvalue no.', printed, ' = ', Eigenvalues[i] * 27211.6, file = text_file)
        indexes.append(i)
        printed += 1
text_file.close()
print(indexes)

#Normalisation <c_n|S|c_n> = 1
# for i in range(15):
#     normalisation_constant = np.matmul(Eigenvectors[indexes[i]], np.matmul(S_matrix, Eigenvectors[indexes[i]]))
#     Eigenvectors[indexes[i]] = Eigenvectors[indexes[i]] / normalisation_constant


for eigen_number in range(6):
    #print(indexes[eigen_number])
    PSI_from_MES = MES.MES_from_vector(dxi,N,a, m, omega, a_b, Eigenvectors[:,int(indexes[eigen_number])])
    #text_file = open("Eigenstates" + str(eigen_number) + ".txt", "w")
    with open("Eigenstates" + str(eigen_number) + ".txt", "w") as text_file:
        for i in range((2*N)*int(2/dxi)) :
            for j in range((2*N)*int(2/dxi)) :
                print(PSI_from_MES[i,j,0], PSI_from_MES[i,j,1], PSI_from_MES[i,j,2], file=text_file)
    #text_file.close()


########################## TIME EVOLUTION ###########################
if time_evolution:
    theoretical_period = 2*np.pi/(Eigenvalues[int(indexes[1])] - Eigenvalues[int(indexes[0])])
    print('Theoretical period is:',theoretical_period, '[atomic units]')
    dt = 100
    total_time = 5e4
    t_steps = int(total_time/dt)
    X_matrix = np.zeros(((2*N+1)**2, (2*N+1)**2))
    x_t = np.zeros(t_steps, dtype=complex)
    D = np.zeros((t_steps, (2*N+1)**2), dtype=complex)
    D[0,:] =  Eigenvectors[:,indexes[0]] + Eigenvectors[:,indexes[5]]
    time_factor = complex(0, -dt/2)

    R_matrix = S_matrix + time_factor*H_matrix
    L_matrix = S_matrix - time_factor*H_matrix
    L_inverse_matrix = np.linalg.inv(L_matrix)

    LR_matrix = np.matmul(L_inverse_matrix, R_matrix)
    with open("Eigenstates_time.dat", "w") as text_file:
        for t in range(t_steps-1):
            D[t+1,:] = np.matmul(LR_matrix, D[t,:])

            if not(t%10) and t*dt <  theoretical_period:
                PSI_from_MES = MES.MES_from_vector(dxi,N,a, m, omega, a_b, np.abs(D[t,:])**2)
                #text_file = open("Eigenstates" + str(eigen_number) + ".txt", "w")
                for i in range((2*N)*int(2/dxi)) :
                    for j in range((2*N)*int(2/dxi)) :
                        print(PSI_from_MES[i,j,0], PSI_from_MES[i,j,1], PSI_from_MES[i,j,2], file=text_file)
                print("\n\n", file=text_file)
                
            


    for element in range(1, (2*N)**2 + 1):
        for i in range(1,5):
            for j in range(1,5):
                X_matrix[Mesh.find_global_number(element,i,N) - 1, Mesh.find_global_number(element,j,N) - 1] += MES.x_operator(element, i, j,a/a_b,N)


    number_of_zeros = 0
    print_period_flag = True 
    with open("X.dat", "w") as text_file:     
        for i in range(t_steps):
            x_t[i] = np.matmul(np.conj(D[i,:]), np.matmul(X_matrix, D[i,:]))
            if x_t[i]*x_t[i-1] < 0 and number_of_zeros == 0:
                time_first = i*dt 
                number_of_zeros += 1
            elif x_t[i]*x_t[i-1] < 0:
                number_of_zeros += 1 

            if number_of_zeros == 3 and print_period_flag:
                print('Period from simulation: ', i*dt - time_first, '[atomic units]')
                print_period_flag = False 
            print( i *dt, x_t[i].real*a_b, file=text_file)        
