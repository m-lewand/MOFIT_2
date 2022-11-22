import Mesh 
import numpy as np 
import MES 

def Parabolic(element, i, j, m, omega, a, N):
    v = 0
    w_gauss = [5./9. , 8./9. , 5./9.]
    p_gauss = [-np.sqrt(3./5), 0, np.sqrt(3./5)]
    for l in range(0,3):
        for n in range(0,3):
            v += a**2/4.*m*omega**2/2. * w_gauss[l]*w_gauss[n]* MES.shape_function(p_gauss[l],p_gauss[n],j)* \
            MES.shape_function(p_gauss[l],p_gauss[n],i)*((Mesh.find_real_space_coordinate(Mesh.find_global_number(element, 1,N), N,a)[0]/2 * (1-p_gauss[l]) + Mesh.find_real_space_coordinate(Mesh.find_global_number(element, 2,N), N,a)[0]/2 * (1+p_gauss[l]) )**2 \
                + (Mesh.find_real_space_coordinate(Mesh.find_global_number(element, 1,N), N,a)[1]/2 * (1-p_gauss[n]) + Mesh.find_real_space_coordinate(Mesh.find_global_number(element, 2,N), N,a)[1]/2 * (1+p_gauss[n]) )**2)
    return v 