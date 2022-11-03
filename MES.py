import numpy as np 
import Mesh


def MES(dxi,N,a,m,omega, a_b):
    n_steps = int(2/dxi) #+ 1
    psi_realspace = np.zeros((n_steps*(2*N) , n_steps*(2*N),3))
    for element_row in range(1,2*N+1):  #iterates over element rows 
        for j in range(n_steps):    #iterates over dxi_2
            for element_in_row in range(1,2*N+1):  #iterates over element in current element row
                element_number = element_in_row + (element_row-1)*2*N   
                for i in range(n_steps):    #iterates over dxi_1
                    xi_1 = -1 + i*dxi #local coordinates
                    xi_2 = -1 + j*dxi
                    x = Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,1,N),N,a)[0]*(1-xi_1)/2 + \
                        Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,2,N),N,a)[0]*(1+xi_1)/2
                    y =  Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,1,N),N,a)[1]*(1-xi_2)/2 + \
                        Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,3,N),N,a)[1]*(1+xi_2)/2
                    #print(x,y)
                    psi = 0
                    for local in range(1,5):
                        psi += np.exp(-m*omega/2*((Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,local,N),N,a)[0]/a_b)**2 + \
                            (Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,local,N),N,a)[1]/a_b)**2))*shape_function(xi_1,xi_2,local)
                    psi_realspace[i+(element_in_row-1)*n_steps, j + (element_row-1)*n_steps,0] = x
                    psi_realspace[i+(element_in_row-1)*n_steps, j + (element_row-1)*n_steps,1] = y
                    psi_realspace[i+(element_in_row-1)*n_steps, j + (element_row-1)*n_steps,2] = psi
    return psi_realspace

def MES_from_vector(dxi,N,a,m,omega, a_b, Eigenvector):
    n_steps = int(2/dxi) #+ 1
    psi_realspace = np.zeros((n_steps*(2*N) , n_steps*(2*N),3))
    for element_row in range(1,2*N+1):  #iterates over element rows 
        for j in range(n_steps):    #iterates over dxi_2
            for element_in_row in range(1,2*N+1):  #iterates over element in current element row
                element_number = element_in_row + (element_row-1)*2*N   
                for i in range(n_steps):    #iterates over dxi_1
                    xi_1 = -1 + i*dxi #local coordinates
                    xi_2 = -1 + j*dxi
                    x = Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,1,N),N,a)[0]*(1-xi_1)/2 + \
                        Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,2,N),N,a)[0]*(1+xi_1)/2
                    y =  Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,1,N),N,a)[1]*(1-xi_2)/2 + \
                        Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,3,N),N,a)[1]*(1+xi_2)/2
                    #print(x,y)
                    psi = 0
                    for local in range(1,5):                    
                        psi += Eigenvector[Mesh.find_global_number(element_number,local,N) - 1]*shape_function(xi_1,xi_2,local)                       
                            
                    psi_realspace[i+(element_in_row-1)*n_steps, j + (element_row-1)*n_steps,0] = x
                    psi_realspace[i+(element_in_row-1)*n_steps, j + (element_row-1)*n_steps,1] = y
                    psi_realspace[i+(element_in_row-1)*n_steps, j + (element_row-1)*n_steps,2] = psi
    return psi_realspace


def f(xi,function_number):
    if function_number == 1:
        return (1-xi)/2
    elif function_number == 2:
        return (1+xi)/2
    else:
        print('Error in shape function f_1')


def shape_function(xi_1, xi_2,local):
    if local == 1:
        return f(xi_1,1)*f(xi_2,1)
    elif local == 2:
        return f(xi_1,2)*f(xi_2,1)
    elif local == 3:
        return f(xi_1,1)*f(xi_2,2)
    elif local == 4:
        return f(xi_1,2)*f(xi_2,2)
    else:
        
        print('\n !!!!!!!!!!!!!!!!!! ERROR IN SHAPE FUNCTIONS !!!!!!!!!!!!!!!!!! \n')

def local_overlap_matrix(overlap_matrix,a):
    w_gauss = [5/9 , 8/9 , 5/9]
    p_gauss = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
    for i in range(0,4):    #iterates over 4x4 matrix elements
        for j in range(0,4):    #iterates over 4x4 matrix elements
            for l in range(0,3):    #iterates over points in Gauss quadrature
                for n in range(0,3):    #iterates over points in Gauss quadrature
                    overlap_matrix[i,j] +=  a**2/4 * w_gauss[l] * w_gauss[n] * shape_function(p_gauss[l],p_gauss[n],j+1)* \
                                            shape_function(p_gauss[l],p_gauss[n],i+1)

def local_kinetic_energy_matrix(kinetic_matrix,m):
    dxi = (1. - np.sqrt(3./5))/5
    w_gauss = [5./9. , 8./9. , 5./9.]
    p_gauss = [-np.sqrt(3./5), 0, np.sqrt(3./5)]
    for i in range(0,4):    #iterates over 4x4 matrix elements
        for j in range(0,4):    #iterates over 4x4 matrix elements
            for l in range(0,3):    #iterates over points in Gauss quadrature
                for n in range(0,3):    #iterates over points in Gauss quadrature
                    kinetic_matrix[i,j] += w_gauss[l]*w_gauss[n]/(2*m)* \
                                    ( (shape_function(p_gauss[l]+dxi,p_gauss[n],j+1) - shape_function(p_gauss[l]-dxi,p_gauss[n],j+1))/(2*dxi) * \
                                    (shape_function(p_gauss[l]+dxi,p_gauss[n],i+1) - shape_function(p_gauss[l]-dxi,p_gauss[n],i+1))/(2*dxi) + \
                                    (shape_function(p_gauss[l],p_gauss[n]+dxi,j+1) - shape_function(p_gauss[l],p_gauss[n]-dxi,j+1))/(2*dxi) * \
                                    (shape_function(p_gauss[l],p_gauss[n]+dxi,i+1) - shape_function(p_gauss[l],p_gauss[n]-dxi,i+1))/(2*dxi) )
