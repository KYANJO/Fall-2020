def weights(z, x, m):
    """
    weights(z, x, m)

    Calculates finite difference weights of up to order m.

    Implements Fornberg's algorithm.

    ARGS:
        z : Location where approximations are to be accurate
        x : Vector with x-coordinates for the grid points
        m : Highest derivative that we want to find weights for

    RETURNS:
        c : Array of size [m+1, len(x)] containing (as output) in successive rows the weights for derivatives 0, 1, ..., m.

    EXAMPLE:
        To generate the 2nd order centered FD formula for the zeroth, first and second derivative, we make the following call to weights:
        c = weights(0, [-1, 0, 1], 2)

    (c) Translated by Andrew Jones from the original source by Fornberg
    """
    import numpy as np

    n = len(x)
    c = np.zeros((m+1, n))
    c1, c4 = 1, x[0] - z
    c[0,0] = 1
    for i in range(1,n):
        mn = min(i+1, m+1)
        c2, c4, c5 = 1, x[i]-z, c4
        for j in range(0, i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j==i-1:
                c[1:mn,i] = c1/c2 *(np.arange(1,mn)*c[0:mn-1,i-1] - c5*c[1:mn, i-1])
                c[0,i] = -c1*c5/c2 * c[0,i-1]
            c[1:mn,j] = (c4*c[1:mn,j] - np.arange(1,mn)*c[0:mn-1,j])/c3
            c[0,j] *= c4/c3
        c1 = c2
    return c
