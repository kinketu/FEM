# FEM module


def dis_to_force(diris, numberOfNodes, K, b):
    """Convert displacement constraint into force constraint.
    Assumption diris = [[node number, value], ...]
    example: diris = [[0, 0]]
    Assumption matrix K is array type.
    """

    for diri in diris:
        i = diri[0]
        for j in range(numberOfNodes):
            if j != i:
                b[j] = b[j] - K[j, i] * diri[1]
            else:
                b[j] = diri[1]
            K[i, j] = 0
            K[j, i] = 0
        K[i, i] = 1
    return K, b

def dis_to_force2(diris, K, b):
    """Convert displacement constraint into force constraint.
    Assumption diris = [[node number, value], ...]
    example: diris = [[0, 0]]
    Assumption matrix K is array type.
    """

    numberOfNodes = K.shape[0]

    for diri in diris:
        i = diri[0]
        for j in range(numberOfNodes):
            if j != i:
                b[j] = b[j] - K[j, i] * diri[1]
            else:
                b[j] = diri[1]
            K[i, j] = 0
            K[j, i] = 0
        K[i, i] = 1
    return K, b
