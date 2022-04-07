import sympy as sym
from sympy.abc import F, z, s

# Throughout this mini-project we follow the notations of ["An O(s^r)-resolution ODE framework for understanding discrete-time algorithms...", Haihao Lu, 2021]
# F     -   vector field
# s     -   timestep
# z     -   variable

# F = sym.Function('F')(z)

class DTA:
    def __init__(self, upd):
        super().__init__()
        if isinstance(upd, str):
            self.upd = sym.sympify(upd)
        # assert isinstance(upd, sym.Function)
        self._ODE_coeff_funs = None
        self._upd_derivs = None
        self._aux_coeff_funs = None

    def __str__(self):
        return self.upd.__str__()

    def step(self, z0, field, steps=1):
        raise NotImplementedError

    def O1_ODE(self):
        raise NotImplementedError



    def Os_ODE(self):
        raise NotImplementedError

    def Osr_ODE(self, r):
        ODE_coeff_funs = [None] * (r+1) # f_i(Z) in the paper Eq.(13)
        upd_derivs = [None] * (r+2) # g_i(z, s) in the paper Eq.(10)
        aux_coeff_funs = [ [None] * (r+2) for _ in range(r+2) ] # h_{ji}(Z) in the paper Eq.(12)
        # TODO: adjust the sizes so the code looks more serious (I only need r+1 in one of the dimensions I think)
        # TODO: check that the equations in the paper are correct (I blindly copied them into code)

        for j in range(r+2):
            upd_derivs[j] = sym.diff(self.upd, s, j)
            upd_derivs[j] = upd_derivs[j].subs(s, 0)
            upd_derivs[j] = upd_derivs[j] / sym.factorial(j)
            upd_derivs[j] = sym.simplify(upd_derivs[j])

        # j = 0
        aux_coeff_funs[0][0] = z
        for i in range(1, r+1):
            aux_coeff_funs[0][i] = 0

        # # j = 1
        # for i in range(1, r+1):
        #     aux_coeff_funs[1][i] = ODE_coeff_funs[i] # impossible
        # # j >= 2
        # for j in range(2, r+2):
        #     for i in range(r+1):
        #         aux_coeff_funs[j][i] = 0
        #         for l in range(i+1):
        #             aux_coeff_funs[j][i] += sym.diff(aux_coeff_funs[j-1][l], z) * aux_coeff_funs[1][i-l]
        #         aux_coeff_funs[j][i] = sym.simplify(aux_coeff_funs[j][i])
        
        # for i in range(r+1):
        #     ODE_coeff_funs[i] = upd_derivs[i+1]
        #     for l in range(2, i+2):
        #         ODE_coeff_funs[i] -= aux_coeff_funs[l][i+1-l] / sym.factorial(l)
        #     ODE_coeff_funs[i] = sym.simplify(ODE_coeff_funs[i])

        ## next we compute f_i and h_{l, i+1-l} (l=2:(i+1)) jointly, recursively for all i, by Eqs.(13) and (12)
        for i in range(0, r+1):
            for l in range(2, i+2):
                aux_coeff_funs[l][i+1-l] = 0
                for k in range(i+1-l + 1):
                    aux_coeff_funs[l][i+1-l] += sym.diff(aux_coeff_funs[l-1][k], z) * aux_coeff_funs[1][i+1-l - k]
                    # print(i, l, k)
                    # print(aux_coeff_funs[l-1][k])
                aux_coeff_funs[l][i+1-l] = sym.simplify(aux_coeff_funs[l][i+1-l])
            ODE_coeff_funs[i] = upd_derivs[i+1]
            for l in range(2, i+2):
                ODE_coeff_funs[i] -= aux_coeff_funs[l][i+1-l] / sym.factorial(l)
            ODE_coeff_funs[i] = sym.simplify(ODE_coeff_funs[i])
            aux_coeff_funs[1][i] = ODE_coeff_funs[i]
            
        if self._ODE_coeff_funs is None or len(self._ODE_coeff_funs) < r+1:
            self._ODE_coeff_funs = ODE_coeff_funs
        if self._upd_derivs is None or len(self._upd_derivs) < r+1:
            self._upd_derivs = upd_derivs
        if self._aux_coeff_funs is None or len(self._aux_coeff_funs) < r+2:
            self._aux_coeff_funs = aux_coeff_funs

        return ODE_coeff_funs

# Note: in the paper they use -F instead of F
GDA = DTA("z + s * F(z)")
EGM = DTA("z + s * F( z + s * F(z) )")


