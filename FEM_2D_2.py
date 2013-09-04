# Unsteady analysis 2D FEM program.
# In this version, add a function of unsteady analysis.
# SI Units

from numpy import array, zeros
from numpy.linalg import solve, inv
from openacoustics.gmsh import *
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
from matplotlib import animation
import fem
import FEM_2D_1


def right_member(x, y):
    return 1

def triangularElement(node0, node1, node2):
    b01 = node0[1]-node1[1]
    b12 = node1[1]-node2[1]
    b20 = node2[1]-node0[1]
    a10 = node1[0]-node0[0]
    a21 = node2[0]-node1[0]
    a02 = node0[0]-node2[0]

    mean_coord = (node0 + node1 + node2)/3

    delta = (a10*b20-b01*a02)/2

    Kn = array([[b12*b12+a21*a21, b12*b20+a21*a02, b12*b01+a21*a10],
                [b20*b12+a02*a21, b20*b20+a02*a02, b20*b01+a02*a10],
                [b01*b12+a10*a21, b01*b20+a10*a02, b01*b01+a10*a10]]) \
                    /(4*delta)

    Mn = array([[1/6., 1/12., 1/12.],
                [1/12., 1/6., 1/12.],
                [1/12., 1/12., 1/6.]])

    bn = array([1, 1, 1])*delta/3*right_member(mean_coord[0], mean_coord[1])

    return Kn, Mn, bn

def totalmatrix(nodes, elements, heat):
    # Create total matrix.
    numberOfNodes = nodes.shape[0]
    K = zeros((numberOfNodes, numberOfNodes))
    M = zeros((numberOfNodes, numberOfNodes))
    b = zeros(numberOfNodes)

    for element in elements:
        node0 = nodes[element[0], :]
        node1 = nodes[element[1], :]
        node2 = nodes[element[2], :]

        # Calculate elementary matrices.
        Kn, Mn, bn = triangularElement(node0, node1, node2)
        Kn = heat*Kn

        # Assemble global system matrices.
        for rowIndex in range(Kn.shape[0]):
            for columnIndex in range(Kn.shape[1]):
                K[element[rowIndex], element[columnIndex]]\
                    += Kn[rowIndex, columnIndex]

        for rowIndex in range(Mn.shape[0]):
            for columnIndex in range(Mn.shape[1]):
                M[element[rowIndex], element[columnIndex]]\
                    += Mn[rowIndex, columnIndex]

        for Index in range(bn.shape[0]):
            b[element[Index]] += bn[Index]

    return K, M, b

def neumann():
    pass

def dirichlet():
    # Take account of Dirichlet boundary condition.
    # Return a list : diris = [[node number,value], [...]]
    # example : diri = [[0,0]]
    i = 0
    diris = []
    for node in nodes:
        #if node[1]==0:
        if (node[1] == 0 or node[0] == 0 or node[1] == 0.5 or node[0] == 0.5):
        #if (node[0]==0 or node[1]==0):
            diris.append([i, 0])
        i += 1

    return diris

def unsteady(u):
    pass

def color_value(slv):
    # Create color value.
    val = []
    for element in elements:
        p0 = slv[element[0]]
        p1 = slv[element[1]]
        p2 = slv[element[2]]
        val.append((p0+p1+p2)/3.)
    val = array(val)

    return val

def postprocessor(nodes, val):
    x = nodes[:, 0]
    y = nodes[:, 1]

    fig = plt.figure()
    plt.gca().set_aspect('equal')
    plt.tripcolor(x,y,elements,facecolors=val,edgecolors='k')
    plt.colorbar()
    plt.title("Poisson's equation")
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.savefig("FEM.png")
    plt.show()


if __name__ == '__main__':
    # Preprocessor
    nodes, elements = FEM_2D_1.preprocessor(True)


    # Solver (steady analysis)
    heat = 1 # Heat conductivity
    numberOfNodes = nodes.shape[0]
    K, M, b = totalmatrix(nodes, elements, heat)
    diris = dirichlet() # Take account of Dirichlet boundary condition.
    #K, b = fem.dis_to_force(diris, numberOfNodes, K, b)
    K, b = fem.dis_to_force2(diris, K, b)
    u = solve(K, b)

    """
    # Solver (unsteady analysis)
    dt = 0.1
    Md = M*(1/dt)
    Kd = K*(1/2)
    for i in range(20):
        u = dot(inv(Md+Kd), b+dot(Md-Kd, u))
        #postprocessor(nodes, u)
    """

    # Postprocessor
    val = color_value(u)
    postprocessor(nodes, val) # Visualize solution.


    # Print number of elements.
    print "number of elements", elements.shape[0]
