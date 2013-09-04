# This program is to solve Linear elastic problem using FEM solver.
# Use Triangular Element.
# Asumption plane strain analysis.
# Units of pressure to use MPa
# Units of force to use MN

from numpy import array, zeros, dot
from numpy.linalg import solve
import matplotlib.pyplot as plt
from openacoustics.gmsh import *
import fem

def preprocessor(debug=True):
    if debug == True:
        # Create nodes and element.
        nodes = array([[0., 0.], [0., 1.], [1., 0.], [1., 1.]])
        elements = array([[0, 2, 1], [3, 1, 2]])

    elif debug == False:
        # Create nodes and element.
        # Create an instance of class gmsh.
        gmsh = GMsh2D()

        # Load the .geo file.
        geo = 'ElasticProblem.geo'
        gmsh.loadGeo(geo, 0.2, order=1)

        # Load nodes indices and coordinates.
        nodes = gmsh.getNodes()

        # Get the sorce node number.
        source = gmsh.getPoints("S")

        # Load node numbers of triangular elements.
        elements = gmsh.getTriangles("", order=1)

    return nodes, elements

def element_matrix(node0, node1, node2):
    t = 1 # t is thickness.
    b01 = node0[1]-node1[1]
    b12 = node1[1]-node2[1]
    b20 = node2[1]-node0[1]
    a10 = node1[0]-node0[0]
    a21 = node2[0]-node1[0]
    a02 = node0[0]-node2[0]

    detA =(a10*b20-b01*a02)

    Be = array([[b12, 0., b20, 0., b01, 0.],
                [0., a21, 0., a02, 0., a10],
                [a21, b12, a02, b20, a10, b01]])/detA

    Ke = t*dot(dot(Be.T, D), Be)*detA/2

    return Ke, Be

def node_matrix(d, element):
    dn = []
    for node in element:
        de.append([node*2, node*2+1])
    dn = array(dn).flatten()

    return dn

def createDmatrix(PlaneType="Strain"):
    if PlaneType == "Strain":
        # Plane Strain
        D = array([[1-Poisson, Poisson, 0.],
                   [Poisson, 1-Poisson, 0.],
                   [0., 0., (1-2*Poisson)/2.]])\
                       *Young/((1+Poisson)*(1-2*Poisson))

    elif PlaneType == "Stress":
        # Plane Stress
        D = array([[1., Poisson, 0.],
                   [Poisson, 1., 0.],
                   [0., 0., (1-Poisson)/2]])*Young/(1-Poisson**2)
    

    return D

def totalmatrices(nodes, elements):
    # Create total matrix.
    numberOfNodes = nodes.shape[0]
    numberOfElements = elements.shape[0]
    K = zeros((numberOfNodes*2, numberOfNodes*2))
    B = []

    for element in elements:
        node0 = nodes[element[0], :]
        node1 = nodes[element[1], :]
        node2 = nodes[element[2], :]

        # Calculate elementary matrices.
        Ke, Be = element_matrix(node0, node1, node2)
        B.append(Be)

        # Assemble global system matrices.
        # Convert local number into global number.
        list1 = []
        list2 = []
        for i in range(numberOfNodes):
            list1.append([i*2, i*2+1])
        for i in element:
            list2.append(list1[i])
        list2 = array(list2).flatten()
        #print list2
        
        for rowIndex in range(Ke.shape[0]):
            for columnIndex in range(Ke.shape[1]):
                K[list2[rowIndex], list2[columnIndex]]\
                    += Ke[rowIndex, columnIndex]

    B = array(B)

    return K, B


def total_strain():
    eps = []
    for element in elements:
        element_number = list(elements).index(list(element))
        dn = node_matrix(d, element)
        epse = dot(B[element_number], dn)
        eps.append(epse)

    eps = array(eps)

    return eps


def mizes_stress(stress):
    sigx = stress[0]
    sigy = stress[1]
    sigz = stress[2]
    tauxy = stress[3]
    tauyz = stress[4]
    tauzx = stress[5]

    mizes = sqrt(((sigx-sigy)**2+(sigy-sigz)**2+(sigz-sigx)**2\
                      +6*(tauxy**2+tauyz**2+tauzx**2))/2)

    return mizes

def postprocessor(nodes, saveimage=False):
    pass

def plot_displacement(old_nodes, new_nodes, saveimage=False):
    x1 = old_nodes[:, 0]
    y1 = old_nodes[:, 1]
    x2 = new_nodes[:, 0]
    y2 = new_nodes[:, 1]

    plt.figure()
    plt.gca().set_aspect('equal')
    plt.triplot(x1, y1, elements, 'ko-')
    plt.triplot(x2, y2, elements, 'bo--')
    plt.title("Elastic Problem")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.1, 1.1)

    if saveimage == True:
        plt.savefig("Elastic_Problem.png")
        
    plt.show()

def plot_mizes(nodes, elements, mizes, saveimage=False):
    x = nodes[:, 0]
    y = nodes[:, 1]

    plt.figure()
    plt.gca().set_aspect('equal')
    plt.tripcolor(x, y, elements, facecolors=mizes, edgecolors='k')
    plt.colorbar()
    plt.title("von-Mizes stress")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('tight')

    if saveimage == True:
        plt.savefig("mizes.png")
        
    plt.show()


if __name__ == '__main__':
    # Preprocessor
    nodes, elements = preprocessor()


    # Solver
    # Young's modulus and Poisson's ration of soft iron
    Young = 200*10**3 # MPa
    Poisson = 0.3
    D = createDmatrix()
    K, B = totalmatrices(nodes, elements)

    unval = 100 # unkown value
    f = array([unval, unval, unval, 0, 10, unval, 10, 0])
    # "discon" is a list of displcement constraint.
    # ex) discon = [[index number in the K matrix, value], ...]
    discon = [[0, 0], [1, 0], [2, 0], [5, 0]]

    K, f = fem.dis_to_force2(discon, K, f)
    d = solve(K, f)
    d = d.reshape((d.shape[0]/2, 2))
    dis_scale = 500 # displacement scale
    new_nodes = nodes + d * dis_scale

    #eps = total_strain()


    # Print Debug
    #print "matrixK\n", K
    #print "displacement\n", d


    # Postprocessor
    #plot_displacement(nodes, new_nodes)
    print d
