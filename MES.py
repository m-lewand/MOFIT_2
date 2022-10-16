import numpy as np 
import Mesh

def MES(dxi,N,a,m,omega):
    n_steps = int(2/dxi) + 1
    psi_realspace = np.zeros((n_steps*(2*N + 1),n_steps(2*N + 1)))
    for element_row in range(1,2*N+1):  #iterates over element rows 
        for j in range(n_steps):    #iterates over dxi_2
            for element_in_row in range(1,2*N+1):  #iterates over element in current element row
                element_number = element_in_row + (element_row-1)*2*N   
                for i in range(n_steps):    #iterates over dxi_1
                    xi_1 = -1 + i*dxi
                    xi_2 = -1 + j*dxi
                    x = Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,1,N),N,a)[0]*(1-xi_1)/2 + \
                        Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,2,N),N,a)[0]*(1+xi_1)/2
                    y =  Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,1,N),N,a)[1]*(1-xi_2)/2 + \
                        Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,3,N),N,a)[1]*(1+xi_2)/2
                    #print(x,y)
                    psi = 0
                    for local in range(4):
                        psi += np.exp(m*omega/2*(Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,local,N),N,a)[0]**2 + \
                            Mesh.find_real_space_coordinate(Mesh.find_global_number(element_number,local,N),N,a)[1]**2))*shape_function(xi_1,xi_2,local)
                    psi_realspace[i+(element_in_row-1)*n_steps, j + (element_row-1)*n_steps] = (x,y,psi)
                    print('chuj')
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
        return f(xi_1,1)+f(xi_2,1)
    elif local == 2:
        return f(xi_1,2)*f(xi_2,1)
    elif local == 3:
        return f(xi_1,1)*f(xi_2,2)
    elif local == 4:
        return f(xi_1,2)*f(xi_2,2)
    else:
        print('Error in shape functions')


