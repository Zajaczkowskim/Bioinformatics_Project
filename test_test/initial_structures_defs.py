import numpy as np
from utils import *
import random
from hilbert import decode, encode

def line(n):
    points = []
    for i in range(n):
        points.append([i, 0, 0])
    return np.array(points)

def random_walk(n):
    points=[[0,0,0] for _ in range(n) ]
    for i in range(1,n):
        val_x = random.randint(-1,1)
        val_y = random.randint(-1,1)
        val_z = random.randint(-1,1)
        base = points[i-1]
        points[i] = [base[0]+val_x,base[1]+val_y,base[2]+val_z]
    return points

def random_walk2(n):
    steps_x = np.random.normal(loc=0, scale=1, size=n)
    steps_y = np.random.normal(loc=0, scale=1, size=n)
    steps_z = np.random.normal(loc=0, scale=1, size=n)
    positions_x = np.cumsum(steps_x)
    positions_y = np.cumsum(steps_y)
    positions_z = np.cumsum(steps_z)
    return np.column_stack((positions_x,positions_y,positions_z))


def hilbert_curve2d(n):
    points=[[0,0,0]] # starting point

    instructions_dict = {"A": ["D","A","A","B"],
                         "B": ["C","B","B","A"],
                         "C": ["B","C","C","D"],
                         "D": ["A","D","D","C"]} # production rules from Wikipedia
    
    n_iter = np.ceil(n/4) #number of blocks needed to connect given number of points (n)
    stack = ["A"]  # list of blocks
    
    # defining blocks
    while len(stack) < n_iter:
        instructions_list=[]
        while len(stack) > 0:
            instructions_list = np.concatenate((instructions_list,np.array(instructions_dict[stack.pop()])))
        stack = instructions_list[::-1].tolist()

    i = 1
    while i <= n_iter:
        current_stage = stack.pop() #get current block (stage)
        base = points[-1] #get the base point, we will start drawing from start's coordinates

        if current_stage == 'A': #draw A
            points.append([base[0], base[1]+1,base[2]])
            points.append([base[0]+1, base[1]+1 ,base[2]])
            points.append([base[0]+1, base[1],base[2]])
            base = points[-1]
            next_stage = stack[-1] if len(stack) >0 else None
            if next_stage  == "A" or next_stage == "D": # prepare base for next stage
                points.append([base[0]+1, base[1],base[2]])
            else:
                points.append([base[0], base[1]-1,base[2]])

        if current_stage == 'B':
            points.append([base[0]-1, base[1],base[2]])
            points.append([base[0]-1, base[1]-1,base[2]])
            points.append([base[0], base[1]-1,base[2]])
            base = points[-1]
            next_stage = stack[-1] if len(stack) >0 else None
            if next_stage  == "B" or next_stage == "C":
                points.append([base[0], base[1]-1,base[2]])
            else:
                points.append([base[0]+1, base[1],base[2]])

        if current_stage == 'C':
            points.append([base[0], base[1]-1,base[2]])
            points.append([base[0]-1, base[1]-1,base[2]])
            points.append([base[0]-1, base[1],base[2]])
            base = points[-1]
            next_stage = stack[-1] if len(stack) >0 else None
            if next_stage  == "C" or next_stage == "B":
                points.append([base[0]-1, base[1],base[2]])
            else:
                points.append([base[0], base[1]+1,base[2]])

        if current_stage == 'D':
            points.append([base[0]+1, base[1],base[2]])
            points.append([base[0]+1, base[1]+1,base[2]])
            points.append([base[0], base[1]+1,base[2]])
            base = points[-1]
            next_stage = stack[-1] if len(stack) >0 else None
            if next_stage  == "D" or next_stage == "A":
                points.append([base[0], base[1]+1,base[2]])
            else:
                points.append([base[0]-1, base[1],base[2]])
        
        i +=1
        
    while len(points) > n:
        points.pop() # delete redundant points
    return points

def hilbert_curve3d(n):
    num_dims = 3

    num_bits = int(np.ceil(np.log2(n)/num_dims)) # number of bits per dimention
    
    max_h = 2**(num_bits*num_dims) # this is the maximum number of nodes in our cube. It wont allways be full 
    
    hilberts = np.arange(n)
    locs = decode(hilberts, num_dims, num_bits)

    return locs    


def helisa(n, radius = 1, step_z = 0.1, n_rotations = 5):
    steps_z = np.repeat(step_z,n)
    angles = np.linspace(0,n_rotations* 2*np.pi, n, endpoint= False)
    x = radius * np.cos(angles)
    y = radius * np.sin(angles)
    z = np.cumsum(steps_z)

    return np.column_stack((x,y,z))

def sphere(n, radius = 1):
    angles_xz = np.random.uniform(0,np.pi, n)
    angles_xy = np.random.uniform(0,2*np.pi, n)
    x = radius * np.cos(angles_xy) * np.cos(angles_xz)
    y = radius * np.cos(angles_xy) * np.sin(angles_xz)
    z = radius * np.sin(angles_xy)

    return np.column_stack((x,y,z))

def cube(n, edge_len = 1):
    x = np.random.uniform(0,edge_len,n)
    y = np.random.uniform(0,edge_len,n)
    z = np.random.uniform(0,edge_len,n)

    return np.column_stack((x,y,z))

def cube_outer(n, edge_len = 1):
    walls = np.random.randint(1,7,n)
    points = []
    for wall in walls:
        x = 0 if wall == 1 else edge_len if wall == 3 else np.random.uniform(0,edge_len)
        y = 0 if wall == 2 else edge_len if wall == 4 else np.random.uniform(0,edge_len)
        z = 0 if wall == 5 else edge_len if wall == 6 else np.random.uniform(0,edge_len)
        points.append([x,y,z])
    return points