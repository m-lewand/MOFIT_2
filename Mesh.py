import math

def create_global_numbers(N,global_numbers):
    #global element number on real mesh 
    for i in range(2*N+1):
        for j in range(2*N+1):
            global_numbers[i,j] = i+j
    return global_numbers

def find_global_number(element_number, local_number, N):
    #element_number and local_number numerated from 1 
    row = math.floor((element_number -1)/ (2*N))
    if local_number == 1:
        global_number = element_number + row 
    elif local_number == 2:
        global_number = element_number + row + 1
    elif local_number == 3:
        global_number = element_number + row + (2*N + 1) + 1
    elif local_number == 4:
        global_number = element_number + row + (2*N + 1)
    else:
        print('Invalid local number')
    return global_number

def find_real_space_coordinate(global_number,N, a):
    y = math.floor((global_number-1)/ (2*N+1))
    x = (global_number + y) % (2*N+1 + 1) - 1
    return ((x-N)*a, (y-N)*a)