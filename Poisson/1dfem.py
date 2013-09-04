# 1D FEM program.

from numpy import array, zeros, mgrid
from numpy.linalg import solve
import numpy as np
import matplotlib.pyplot as plt
import fem


def preprocessor():
    # Prepare nodes and elements that are used in analysis.

    # Coodinates of all nodes
    # example: nodes = array([[0], [0.5], [1]])
    nodes = mgrid[0:1.1:0.1]

    # Elements
    # example: elements = array([[0, 1], [1, 2]])
    elements = []
    for i in range(len(nodes) - 1):
        elements.append([i, i + 1])
    array(elements)

    return nodes, elements

def barElement(node0, node1):
    # Calculate element matrix.
    L = node1 - node0
    Kn = array([[-L / 3 + 1 / L, -L / 6 - 1 / L],
                [-L / 6 - 1 / L, -L / 3 + 1 / L]])
    bn = array([2 * node0 + node1, node0 + 2 * node1]) * L / 6

    return Kn, bn

def totalmatrices(nodes, elements):
    # Create total matrix.
    numberOfNodes = nodes.shape[0]
    K = zeros([numberOfNodes, numberOfNodes])
    b = zeros(numberOfNodes)

    for element in elements:
        node0 = nodes[element[0]]
        node1 = nodes[element[1]]

        # Calculate elementary matrices.
        Kn, bn = barElement(node0, node1)

        # Assemble global system matrices.
        for rowIndex in range(Kn.shape[0]):
            for columnIndex in range(Kn.shape[1]):
                K[element[rowIndex], element[columnIndex]]\
                    += Kn[rowIndex, columnIndex]

        for Index in range(bn.shape[0]):
            b[element[Index]] += bn[Index]

    return K, b

def postprocessor(nodes):
    # Function of visualizing solution.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    x1 = nodes
    y1 = zeros(len(nodes))

    plt.plot(x1, u, 'k')

    x2 = mgrid[0:1:0.01]
    y2 = np.sin(x2) / np.cos(1) - x2

    plt.plot(x2, y2, 'k--')
    plt.title('1D FEM')
    plt.xlabel('x')
    plt.ylabel('u')
    #plt.savefig("1D_Fig3.png")
    plt.show()


if __name__ == '__main__':
    # Preprocessor
    nodes, elements = preprocessor()


    # Solver
    numberOfNodes = nodes.shape[0]
    K, b = totalmatrices(nodes, elements)

    # Take account of Dirichlet boundary condition.
    # Assumption Dirichlet boundary condition u0 = 0
    # Assumption diri = [[node number, value]] example: diri = [[0, 0]]
    diris = [[0, 0]]

    #K, b = fem.dis_to_force(diris, numberOfNodes, K, b)
    K, b = fem.dis_to_force2(diris, K, b)

    u = solve(K, b)


    # Postprocessor
    postprocessor(nodes) # Visualize solution.
