# Unsteady analysis 2D FEM program.
# In this version, add a function of animation of unsteady analysis.
# Can't do animation, do viualize last slv
# SI Units

from numpy import array, zeros
from numpy.linalg import solve, inv
from openacoustics.gmsh import *
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
from matplotlib import animation
import fem
import FEM_2D_1
from FEM_2D_2 import right_member, triangularElement, totalmatrix
from FEM_2D_2 import neumann, unsteady, postprocessor


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



    # Solver (unsteady analysis)
    dt = 0.1
    Md = M*(1/dt)
    Kd = K*(1/2)
    """
    for i in range(20):
        u = dot(inv(Md+Kd), b+dot(Md-Kd, u))
        #postprocessor(nodes, u)
    """


    # Postprocessor
    val = color_value(u)
    fig = plt.figure()
    ims = []
    x = nodes[:, 0]
    y = nodes[:, 1]
    for i in range(10):
        u = dot(inv(Md+Kd), b+dot(Md-Kd, u))
        temp = plt.tripcolor(x, y, elements, facecolors=val, adgecolors='k')
        ims.append(temp)
    im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000,
                                       blit=True)
    plt.show()
    #postprocessor(nodes, val) # Visualize solution.


    # Print number of elements.
    print "number of elements", elements.shape[0]

