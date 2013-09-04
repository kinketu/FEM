# Steady analysis 2D FEM program.
# In this version, add a function of right_member.
# SI Units

from numpy import array, zeros
from numpy.linalg import solve
from openacoustics.gmsh import *
import matplotlib.pyplot as plt
import fem


def preprocessor(debug=False):
    if debug == False:
        # Create nodes and elements : Method 1
        # Create an instance of class gmsh.
        gmsh = GMsh2D()

        # Load the .geo file.
        geo = 'steady2dfem.geo'
        gmsh.loadGeo(geo, 0.1, order=1)

        # Load nodes indices and coordinates.
        nodes = gmsh.getNodes()

        # Get the sorce node number.
        source = gmsh.getPoints("S")

        # Load node numbers of triangular elements.
        elements = gmsh.getTriangles("", order=1)

    elif debug == True:
        # Create nodes and element : Method 2
        # Coodinates of all nodes
        # example: nodes = array([[0., 0.], [0.5, 0.], [0.5, 0.5], [0., 0.5]])
        nodes = array([[0., 0.], [0.5, 0.], [0., 0.5], [0.5, 0.5]])

        # Elements
        # example: elements = array([[0, 1, 3], [1, 2, 3]])
        elements = array([[0, 1, 3], [0, 3, 2]])

    return nodes, elements

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
                [b01*b12+a10*a21, b01*b20+a10*a02, b01*b01+a10*a10]])/(4*delta)

    bn = array([1, 1, 1])*delta/3*right_member(mean_coord[0], mean_coord[1])

    return Kn, bn

def totalmatrices(nodes, elements):
    # Create total matrix.
    numberOfNodes = nodes.shape[0]
    K = zeros((numberOfNodes, numberOfNodes))
    b = zeros(numberOfNodes)

    for element in elements:
        node0 = nodes[element[0], :]
        node1 = nodes[element[1], :]
        node2 = nodes[element[2], :]

        # Calculate elementary matrices.
        Kn, bn = triangularElement(node0, node1, node2)

        # Assemble global system matrices.
        for rowIndex in range(Kn.shape[0]):
            for columnIndex in range(Kn.shape[1]):
                K[element[rowIndex], element[columnIndex]]\
                    += Kn[rowIndex, columnIndex]

        for Index in range(bn.shape[0]):
            b[element[Index]] += bn[Index]

    return K, b

def neumann():
    pass

def dirichlet():
    # Take account of Dirichlet boundary condition.
    # Return a list : diris = [[node number, value], [...]]
    # example : diri = [[0, 0]]
    i = 0
    diris = []
    for node in nodes:
        #if node[1]==0:
        #if (node[1] == 0 or node[0] == 0 or node[1] == 0.5 or node[0] == 0.5):
        if (node[0]==0 or node[1]==0):
            diris.append([i, 0])
        i += 1

    return diris

def postprocessor(nodes, slv):
    x = nodes[:, 0]
    y = nodes[:, 1]

    # Color value
    val=[]
    for element in elements:
        p0 = slv[element[0]]
        p1 = slv[element[1]]
        p2 = slv[element[2]]
        val.append((p0+p1+p2)/3.)
    val = array(val)

    plt.figure()
    plt.gca().set_aspect('equal')
    plt.tripcolor(x, y, elements, facecolors=val, edgecolors='k')
    plt.colorbar()
    plt.title("Poisson's equation")
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.savefig("FEM.png")
    plt.show()
    print "val\n", val


if __name__ == '__main__':
    # Preprocessor
    nodes, elements = preprocessor(True)


    # Solver
    numberOfNodes = nodes.shape[0]
    K, b = totalmatrices(nodes, elements)
    diris = dirichlet() # Take account of Dirichlet boundary condition.
    #K, b = fem.dis_to_force(diris, numberOfNodes, K, b)
    K, b = fem.dis_to_force2(diris, K, b)
    u = solve(K, b)


    # Postprocessor
    postprocessor(nodes, u) # Visualize solution.


    # Print number of elements.
    print "number of elements", elements.shape[0]

