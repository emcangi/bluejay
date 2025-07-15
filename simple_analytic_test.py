import numpy as np

# parameters
n_species = 2
n_alt = 2
n_horiz = 2

# transport coefficients
# tforwards from column i to i+1; same for all altitudes and species for simplicity
Tf = 10.0
Tb = 5.0

# initial densities (species-first, altitude-second, horizontal-third)
# order: sp1-alt1-col1, sp2-alt1-col1, sp1-alt2-col1, sp2-alt2-col1,
#        sp1-alt1-col2, sp2-alt1-col2, sp1-alt2-col2, sp2-alt2-col2
n = np.arange(1, n_species*n_alt*n_horiz + 1, dtype=float)

# function to compute derivative and Jacobian

def derivative_and_jacobian(n_vec):
    d = np.zeros_like(n_vec)
    J = np.zeros((len(n_vec), len(n_vec)))
    idx = 0
    for sp in range(n_species):
        for alt in range(n_alt):
            # column 1 index
            i1 = sp + n_species*alt + 0*n_species*n_alt
            # column 2 index
            i2 = sp + n_species*alt + 1*n_species*n_alt
            n1 = n_vec[i1]
            n2 = n_vec[i2]
            d[i1] = Tb*n2 - Tf*n1
            d[i2] = Tf*n1 - Tb*n2
            # jacobian entries
            J[i1, i1] = -Tf
            J[i1, i2] = Tb
            J[i2, i1] = Tf
            J[i2, i2] = -Tb
    return d, J

if __name__ == "__main__":
    d, J = derivative_and_jacobian(n)
    np.set_printoptions(suppress=True)
    print("n =", n)
    print("d/dt =", d)
    print("Jacobian:\n", J)
