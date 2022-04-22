
import numpy as np

def initialConditions(dr, Ntot, Fr):
    '''Returns the position of the membrne which is the solution
    of the static problem for the membrane. 
    Since A_prime determines the position of the membrane when 
    there are no contact points, this position vector eta should satisfy
    A*eta = eta + R
    
    Parameters:
    ----------
    dr: float
        Spatial grid spacing of the discretized version of the problem.
    Ntot: int
        Number of non-trivial spatial points in the system
    Fr: float
        Frobenius Constant of the system, equal to the density of the membrane, 
        times the radius, times gravity, divided by tension (mu * rS * g) / Tm
    
    Outputs:
    ----------
    eta: np.array, dtype = 'float64'
        Vector representing the position of the membrane and which is 
        in numerical equilibrium in this discretized version (i.e. when passed
        as argument for euler or trapecio, it will return the same vector as output)
    '''

    # Building A_prime

    # First entries in the tridiagonal matrix
    A_prime =  np.eye(Ntot) / (dr**2)
    A_prime[0, 0] = 2/(dr**2)
    A_prime[0, 1] = -2/(dr**2)

    # All the other entries
    for i in range(1, Ntot-1):
        A_prime[i, i-1] = -(1/(2*dr**2)) * (2*i-1)/(2*i)
        A_prime[i, i+1] = -(1/(2*dr**2)) * (2*i+1)/(2*i)
    
    # Last entry
    A_prime[-1, -2] = -1/(2*dr**2) * (2*Ntot - 3)/(2*Ntot-2)
    # Building RHS of the system
    R = -Fr * np.ones(Ntot, 1)
    # Solving the system
    eta = np.linalg.solve(A_prime, R)/2

    return eta