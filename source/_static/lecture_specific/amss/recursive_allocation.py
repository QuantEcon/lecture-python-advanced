class AMSS:
    # WARNING: THE CODE IS EXTREMELY SENSITIVE TO CHOCIES OF PARAMETERS.
    # DO NOT CHANGE THE PARAMETERS AND EXPECT IT TO WORK

    def __init__(self, pref, β, Π, g, x_grid, bounds_v):
        self.β, self.Π, self.g = β, Π, g
        self.x_grid = x_grid
        self.n = x_grid[0][2]
        self.S = len(Π)
        self.bounds = bounds_v
        self.pref = pref

        self.T_v, self.T_w = bellman_operator_factory(Π, β, x_grid, g,
                                                      bounds_v)

        self.V_solved = False
        self.W_solved = False

    def compute_V(self, V, σ_v_star, tol_vfi, maxitr, print_itr):

        T_v = self.T_v

        self.success = False

        V_new = np.zeros_like(V)

        Δ = 1.0
        for itr in range(maxitr):
            T_v(V, V_new, σ_v_star, self.pref)

            Δ = np.max(np.abs(V_new - V))

            if Δ < tol_vfi:
                self.V_solved = True
                print('Successfully completed VFI after %i iterations'
                      % (itr+1))
                break

            if (itr + 1) % print_itr == 0:
                print('Error at iteration %i : ' % (itr + 1), Δ)

            V[:] = V_new[:]

        self.V = V
        self.σ_v_star = σ_v_star

        return V, σ_v_star

    def compute_W(self, b_0, W, σ_w_star):
        T_w = self.T_w
        V = self.V

        T_w(W, σ_w_star, V, b_0, self.pref)

        self.W = W
        self.σ_w_star = σ_w_star
        self.W_solved = True
        print('Succesfully solved the time 0 problem.')

        return W, σ_w_star

    def solve(self, V, σ_v_star,  b_0, W, σ_w_star, tol_vfi=1e-7,
              maxitr=1000, print_itr=10):
        print("===============")
        print("Solve time 1 problem")
        print("===============")
        self.compute_V(V, σ_v_star, tol_vfi, maxitr, print_itr)
        print("===============")
        print("Solve time 0 problem")
        print("===============")
        self.compute_W(b_0, W, σ_w_star)

    def simulate(self, s_hist, b_0):
        if not (self.V_solved and self.W_solved):
            msg = "V and W need to be successfully computed before simulation."
            raise ValueError(msg)

        pref = self.pref
        x_grid, g, β, S = self.x_grid, self.g, self.β, self.S
        σ_v_star, σ_w_star = self.σ_v_star, self.σ_w_star

        T = len(s_hist)
        s_0 = s_hist[0]

        # Pre-allocate
        n_hist = np.zeros(T)
        x_hist = np.zeros(T)
        c_hist = np.zeros(T)
        τ_hist = np.zeros(T)
        b_hist = np.zeros(T)
        g_hist = np.zeros(T)

        # Compute t = 0
        l_0, T_0 = σ_w_star[s_0]
        c_0 = (1 - l_0) - g[s_0]
        x_0 = (-pref.Uc(c_0, l_0) * (c_0 - T_0 - b_0) +
               pref.Ul(c_0, l_0) * (1 - l_0))

        n_hist[0] = (1 - l_0)
        x_hist[0] = x_0
        c_hist[0] = c_0
        τ_hist[0] = 1 - pref.Ul(c_0, l_0) / pref.Uc(c_0, l_0)
        b_hist[0] = b_0
        g_hist[0] = g[s_0]

        # Compute t > 0
        for t in range(T - 1):
            x_ = x_hist[t]
            s_ = s_hist[t]
            l = np.zeros(S)
            T = np.zeros(S)
            for s in range(S):
                x_arr = np.array([x_])
                l[s] = eval_linear(x_grid, σ_v_star[s_, :, s], x_arr)
                T[s] = eval_linear(x_grid, σ_v_star[s_, :, S+s], x_arr)

            c = (1 - l) - g
            u_c = pref.Uc(c, l)
            Eu_c = Π[s_] @ u_c

            x = u_c * x_ / (β * Eu_c) - u_c * (c - T) + pref.Ul(c, l) * (1 - l)

            c_next = c[s_hist[t+1]]
            l_next = l[s_hist[t+1]]

            x_hist[t+1] = x[s_hist[t+1]]
            n_hist[t+1] = 1 - l_next
            c_hist[t+1] = c_next
            τ_hist[t+1] = 1 - pref.Ul(c_next, l_next) / pref.Uc(c_next, l_next)
            b_hist[t+1] = x_ / (β * Eu_c)
            g_hist[t+1] = g[s_hist[t+1]]

        return c_hist, n_hist, b_hist, τ_hist, g_hist, n_hist


def obj_factory(Π, β, x_grid, g):
    S = len(Π)

    @njit
    def obj_V(σ, state, V, pref):
        # Unpack state
        s_, x_ = state

        l = σ[:S]
        T = σ[S:]

        c = (1 - l) - g
        u_c = pref.Uc(c, l)
        Eu_c = Π[s_] @ u_c
        x = u_c * x_ / (β * Eu_c) - u_c * (c - T) + pref.Ul(c, l) * (1 - l)

        V_next = np.zeros(S)

        for s in range(S):
            V_next[s] = eval_linear(x_grid, V[s], np.array([x[s]]))

        out = Π[s_] @ (pref.U(c, l) + β * V_next)

        return out

    @njit
    def obj_W(σ, state, V, pref):
        # Unpack state
        s_, b_0 = state
        l, T = σ

        c = (1 - l) - g[s_]
        x = -pref.Uc(c, l) * (c - T - b_0) + pref.Ul(c, l) * (1 - l)

        V_next = eval_linear(x_grid, V[s_], np.array([x]))

        out = pref.U(c, l) + β * V_next

        return out

    return obj_V, obj_W


def bellman_operator_factory(Π, β, x_grid, g, bounds_v):
    obj_V, obj_W = obj_factory(Π, β, x_grid, g)
    n = x_grid[0][2]
    S = len(Π)
    x_nodes = nodes(x_grid)

    @njit(parallel=True)
    def T_v(V, V_new, σ_star, pref):
        for s_ in prange(S):
            for x_i in prange(n):
                state = (s_, x_nodes[x_i])
                x0 = σ_star[s_, x_i]
                res = optimize.nelder_mead(obj_V, x0, bounds=bounds_v,
                                           args=(state, V, pref))

                if res.success:
                    V_new[s_, x_i] = res.fun
                    σ_star[s_, x_i] = res.x
                else:
                    print("Optimization routine failed.")

    bounds_w = np.array([[-9.0, 1.0], [0., 10.]])

    def T_w(W, σ_star, V, b_0, pref):
        for s_ in prange(S):
            state = (s_, b_0)
            x0 = σ_star[s_]
            res = optimize.nelder_mead(obj_W, x0, bounds=bounds_w,
                                       args=(state, V, pref))

            W[s_] = res.fun
            σ_star[s_] = res.x

    return T_v, T_w
