# Linear elastic problem with FEM program.
# Experimental code

from numpy import array, zeros, dot
from numpy.linalg import solve
import matplotlib.pyplot as plt
import fem
from ElasticProblem import mizes_stress, plot_mizes

# Preprocessor
nodes = array([[0., 0.], [0., 1.], [1., 0.], [1., 1.]])
elements = array([[0, 2, 1], [3, 1, 2]])

B0 = array([[-1, 0, 1, 0, 0, 0],
            [0, -1, 0, 0, 0, 1],
            [-1, -1, 0, 1, 1, 0]])

B1 = array([[1, 0, -1, 0, 0, 0],
            [0, 1, 0, 0, 0, -1],
            [1, 1, 0, -1, -1, 0]])

D = array([[7, 3, 0],
           [3, 7, 0],
           [0, 0, 2]])*50000/1.3

"""
K0 = array([[9, 5, -7, -2, -2, -3],
            [5, 9, -3, -2, -2, -7],
            [-7, -3, 7, 0, 0, 3],
            [-2, -2, 0, 2, 2, 0],
            [-2, -2, 0, 2, 2, 0],
            [-3, -7, 3, 0, 0, 7]])*50000/2.6

K1 = array([[9, 5, -7, -2, -2, -3],
            [5, 9, -3, -2, -2, -7],
            [-7, -3, 7, 0, 0, 3],
            [-2, -2, 0, 2, 2, 0],
            [-2, -2, 0, 2, 2, 0],
            [-3, -7, 3, 0, 0, 7]])*50000/2.6
"""

# Solver: Method 1
K = array([[9, 5, -2, -3, -7, -2, 0, 0],
           [5, 9, -2, -7, -3, -2, 0, 0],
           [-2, -2, 9, 0, 0, 5, -7, -3],
           [-3, -7, 0, 9, 5, 0, -2, -2],
           [-7, -3, 0, 5, 9, 0, -2, -2],
           [-2, -2, 5, 0, 0, 9, -3, -7],
           [0, 0, -7, -2, -2, -3, 9, 5],
           [0, 0, -3, -2, -2, -7, 5, 9]])*50000/2.6

unval = 100 # unknown value
f = array([unval, unval, unval, 0, 10, unval, 10, 0])

# "discon" is a list of displcement constraint.
# discon = [[index number in the K matrix, value], ...]
discon = [[0, 0], [1, 0], [2, 0], [5, 0]]
K, f = fem.dis_to_force2(discon, K, f)
d1 = solve(K, f)
d1 = d1.reshape((4, 2))
new_nodes1 = nodes + d1 * 500

"""
# Solver: Method 2
# Convert displacement constraint into force constraint.
# NewK is a converted matrix.
newK = array([[1, 0, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0, 0, 0],
              [0, 0, 0, 9, 5, 0, -2, -2],
              [0, 0, 0, 5, 9, 0, -2, -2],
              [0, 0, 0, 0, 0, 1, 0, 0],
              [0, 0, 0, -2, -2, 0, 9, 5],
              [0, 0, 0, -2, -2, 0, 5, 9]])*50000/2.6

newf = array([0, 0, 0, 0, 10, 0, 10, 0])

d2 = solve(newK, newf)
# Must change following. Reshape automatically.
d2 = d2.reshape((4, 2))
new_nodes2 = nodes + d2 * 500
"""

d = d1.flatten()
de0 = zeros(6)
de1 = zeros(6)
de0[0:2] = d[0:2]
de0[2:4] = d[4:6]
de0[4:6] = d[2:4]
de1[0:2] = d[6:8]
de1[2:4] = d[2:4]
de1[4:6] = d[4:6]

eps0 = dot(B0, de0)
eps1 = dot(B1, de1)
print "epsilon0\n", eps0
print "epsilon1\n", eps1

sgm0 = dot(D, eps0)
sgm1 = dot(D, eps1)
print "sigma0\n", sgm0
print "sigma1\n", sgm1

stress0 = zeros(6)
stress1 = zeros(6)
stress0[0:2] = sgm0[0:2]
stress0[3] = sgm0[2]
stress1[0:2] = sgm1[0:2]
stress1[3] = sgm1[2]

mizes0 = mizes_stress(stress0)
mizes1 = mizes_stress(stress1)
print "mizes stress0\n", mizes0
print "mizes stress1\n", mizes1
mizes = array([mizes0, mizes1])
print mizes, mizes.shape

visualize_mizes(new_nodes1, elements, mizes)

# Postprocessor
def postprocessor(nodes):
    x = nodes[:, 0]
    y = nodes[:, 1]

    plt.figure()
    plt.gca().set_aspect('equal')
    plt.triplot(x, y, elements, 'ko-')
    plt.title("Elastic Problem")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.1, 1.1)
    #plt.savefig("FEM.png")
    plt.show()

if __name__ == '__main__':
    print "-"*10+"Method 2" + "-" * 10
    print "displacement\n", d1
    print "coordinates\n", new_nodes1
    """
    print "-"*10+"Method 2" + "-" * 10
    print "displacement\n", d2
    print "coordinates\n", new_nodes2
    """
    #postprocessor(new_nodes)
