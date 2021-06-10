class RecursiveLS:

    '''
    Compute the planner's allocation by solving Bellman
    equation.
    '''

    def __init__(self,
                 pref,
                 x_grid,
                 π=0.5*np.ones((2, 2)),
                 g=np.array([0.1, 0.2])):

        self.π, self.g, self.S = π, g, len(π)
        self.pref, self.x_grid = pref, x_grid

        bounds = np.empty((self.S, 2))

        # bound for n
        bounds[0] = 0, 1

        # bound for xprime
        for s in range(self.S-1):
            bounds[s+1] = x_grid.min(), x_grid.max()

        self.bounds = bounds

        # initialization of time 1 value function
        self.V = None

    def time1_allocation(self, V=None, tol=1e-7):
        '''
        Solve the optimal time 1 allocation problem
        by iterating Bellman value function.
        '''

        π, g, S = self.π, self.g, self.S
        pref, x_grid, bounds = self.pref, self.x_grid, self.bounds

        # initial guess of value function
        if V is None:
            V = np.zeros((len(x_grid), S))

        # initial guess of policy
        z = np.empty((len(x_grid), S, S+2))

        # guess of n
        z[:, :, 1] = 0.5

        # guess of xprime
        for s in range(S):
            for i in range(S-1):
                z[:, s, i+2] = x_grid

        while True:
            # value function iteration
            V_new, z_new = T(V, z, pref, π, g, x_grid, bounds)

            if np.max(np.abs(V - V_new)) < tol:
                break

            V = V_new
            z = z_new

        self.V = V_new
        self.z1 = z_new
        self.c1 = z_new[:, :, 0]
        self.n1 = z_new[:, :, 1]
        self.xprime1 = z_new[:, :, 2:]

        return V_new, z_new

    def time0_allocation(self, b0, s0):
        '''
        Find the optimal time 0 allocation by maximization.
        '''

        if self.V is None:
            self.time1_allocation()

        π, g, S = self.π, self.g, self.S
        pref, x_grid, bounds = self.pref, self.x_grid, self.bounds
        V, z1 = self.V, self.z1

        x = 1. # x is arbitrary
        res = nelder_mead(obj_V,
                          z1[0, s0, 1:-1],
                          args=(x, s0, V, pref, π, g, x_grid, b0),
                          bounds=bounds,
                          tol_f=1e-10)

        n0, xprime0 = IC(res.x, x, s0, b0, pref, π, g)
        c0 = n0 - g[s0]
        z0 = np.array([c0, n0, *xprime0])

        self.z0 = z0
        self.n0 = n0
        self.c0 = n0 - g[s0]
        self.xprime0 = xprime0

        return z0

    def τ(self, c, n):
        '''
        Computes τ given c, n
        '''
        pref = self.pref
        uc, ul = pref.Uc(c, 1-n), pref.Ul(c, 1-n)

        return 1 - ul / uc

    def simulate(self, b0, s0, T, sHist=None):
        '''
        Simulates Ramsey plan for T periods
        '''
        pref, π = self.pref, self.π
        Uc = pref.Uc

        if sHist is None:
            sHist = self.mc.simulate(T, s0)

        cHist, nHist, Bhist, τHist, xHist = np.empty((5, T))
        RHist = np.zeros(T-1)

        # Time 0
        self.time0_allocation(b0, s0)
        cHist[0], nHist[0], xHist[0] = self.c0, self.n0, self.xprime0[s0]
        τHist[0] = self.τ(cHist[0], nHist[0])
        Bhist[0] = b0

        # Time 1 onward
        for t in range(1, T):
            s, x = sHist[t], xHist[t-1]
            cHist[t] = interp(self.x_grid, self.c1[:, s], x)
            nHist[t] = interp(self.x_grid, self.n1[:, s], x)

            τHist[t] = self.τ(cHist[t], nHist[t])

            Bhist[t] = x / Uc(cHist[t], 1-nHist[t])

            c, n = np.empty((2, self.S))
            for sprime in range(self.S):
                c[sprime] = interp(x_grid, self.c1[:, sprime], x)
                n[sprime] = interp(x_grid, self.n1[:, sprime], x)
            Euc = π[sHist[t-1]] @ Uc(c, 1-n)
            RHist[t-1] = Uc(cHist[t-1], 1-nHist[t-1]) / (self.pref.β * Euc)

            gHist = self.g[sHist]
            yHist = nHist
            
            if t < T-1:
                sprime = sHist[t+1]
                xHist[t] = interp(self.x_grid, self.xprime1[:, s, sprime], x)

        return [cHist, nHist, Bhist, τHist, gHist, yHist, xHist, RHist]

# Helper functions

@njit(parallel=True)
def T(V, z, pref, π, g, x_grid, bounds):
    '''
    One step iteration of Bellman value function.
    '''

    S = len(π)

    V_new = np.empty_like(V)
    z_new = np.empty_like(z)

    for i in prange(len(x_grid)):
        x = x_grid[i]
        for s in prange(S):
            res = nelder_mead(obj_V,
                              z[i, s, 1:-1],
                              args=(x, s, V, pref, π, g, x_grid),
                              bounds=bounds,
                              tol_f=1e-10)

            # optimal policy
            n, xprime = IC(res.x, x, s, None, pref, π, g)
            z_new[i, s, 0] = n - g[s]        # c
            z_new[i, s, 1] = n               # n
            z_new[i, s, 2:] = xprime         # xprime
            
            V_new[i, s] = res.fun

    return V_new, z_new

@njit
def obj_V(z_sub, x, s, V, pref, π, g, x_grid, b0=None):
    '''
    The objective on the right hand side of the Bellman equation.
    z_sub contains guesses of n and xprime[:-1].
    '''

    S = len(π)
    β, U = pref.β, pref.U

    # find (n, xprime) that satisfies implementability constraint
    n, xprime = IC(z_sub, x, s, b0, pref, π, g)
    c, l = n-g[s], 1-n

    # if xprime[-1] violates bound, return large penalty
    if (xprime[-1] < x_grid.min()):
        return -1e9 * (1 + np.abs(xprime[-1] - x_grid.min()))
    elif (xprime[-1] > x_grid.max()):
        return -1e9 * (1 + np.abs(xprime[-1] - x_grid.max()))

    # prepare Vprime vector
    Vprime = np.empty(S)
    for sprime in range(S):
        Vprime[sprime] = interp(x_grid, V[:, sprime], xprime[sprime])

    # compute the objective value
    obj = U(c, l) + β * π[s] @ Vprime

    return obj

@njit
def IC(z_sub, x, s, b0, pref, π, g):
    '''
    Find xprime[-1] that satisfies the implementability condition
    given the guesses of n and xprime[:-1].
    '''

    β, Uc, Ul = pref.β, pref.Uc, pref.Ul

    n = z_sub[0]
    xprime = np.empty(len(π))
    xprime[:-1] = z_sub[1:]

    c, l = n-g[s], 1-n
    uc = Uc(c, l)
    ul = Ul(c, l)

    if b0 is None:
        diff = x
    else:
        diff = uc * b0

    diff -= uc * (n - g[s]) - ul * n + β * π[s][:-1] @ xprime[:-1]
    xprime[-1] = diff / (β * π[s][-1])

    return n, xprime
