import math
import itertools
import numpy as np
import sympy


def all_partitions(n):
    """Returns all partitions of a positive integer n"""
    answer = set()
    answer.add((n, ))
    for i in range(1, n):
        for j in all_partitions(n - i):
            answer.add(tuple(sorted((i, ) + j)))
    return answer


def num_inverse(n, coeff_mod):
    if coeff_mod == 0:
        return sympy.Rational(n)**(-1)
    return pow(n, -1, coeff_mod)


def prod(array):
    """Take a product of a list of algebraic objects.
    Annoyingly, math.prod missing from some versions of python...

    :param array: list
    :return: obj
    """
    if len(array) == 0:
        return 1
    output = array.pop(0)
    while len(array) > 0:
        output *= array.pop(0)
    return output


class WPS(object):
    """Weighted projective space"""

    def __init__(self, weights):
        self.weights = weights
        self.dim = len(weights) - 1
        self.lbk = {k: self._get_lbk(k) for k in range(1, len(weights))}
        self.coefs = [self.lbk[1] * self.lbk[k] // self.lbk[k+1] for k in range(1,self.dim)]

    def _get_lbk(self, k):
        subweights = list(itertools.combinations(self.weights, k+1))
        numbers = [math.prod(sw)//math.gcd(*sw) for sw in subweights]
        return math.lcm(*numbers)

    def get_lbk(self, k):
        """Returns l^{b}_{k} in the notation of Kawasaki"""
        return self.lbk[k]

    def get_coeffs(self):
        """Returns coefficients c_{k} satisfying gamma_1 * gamma_k = c_k * gamma_{k+1}
        in the notation of Kawasaki"""
        return self.coefs


class MultiProj(object):
    """Product of projective spaces"""

    def __init__(self, dims):
        self.dims = dims
        self.total_dim = sum(self.dims)
        self.vars = [sympy.Symbol(f"x_{i}") for i in range(len(dims))]
        # total chern class and volume form
        self.chern = 1
        self.volume = 1
        for i, d in enumerate(self.dims):
            self.chern *= ((1 + self.vars[i])**(d + 1)).as_poly()
            self.volume *= self.vars[i] ** d
        self.chern = self.truncate(self.chern)

    def c_to_dict(self, c):
        """Turns cohomology class into dict"""
        return sympy.Poly(c, *self.vars).as_dict()

    def dict_to_c(self, d):
        """Turns dict into cohomology class"""
        return sympy.Poly.from_dict(d, *self.vars)

    def get_mth(self, c, m):
        """Returns grading m summand of an inhomogeneous cohomology class"""
        d = self.c_to_dict(c)
        output_d = {k: v for k, v in d.items() if sum(k)==m}
        return self.dict_to_c(output_d)

    def truncate(self, c):
        """Truncates a polynomial in self.vars to eliminate high order terms which vanish in cohomology"""
        c_dict = self.c_to_dict(c)
        c_out = dict()
        for k, v in c_dict.items():
            if all([k[i] <= self.dims[i] for i in range(len(self.dims))]):
                c_out[k] = v
        return self.dict_to_c(c_out)

    def get_line(self, deg):
        """Returns Chern class of a line bundle specified by its (multi)degree"""
        output = 1
        for i, d in enumerate(deg):
            output += d * self.vars[i]
        return sympy.Poly(output)

    def get_line_inv(self, deg):
        """Returns inverse of total Chern class of a line bundle specified by its (multi)degree"""
        output = 1
        c_1_line = sum([d*self.vars[i] for i, d in enumerate(deg)])
        for i in range(1, self.total_dim + 1):
            output += (-c_1_line) ** i
        return self.truncate(output.as_poly())

    def integrate(self, c):
        """From a polynomial c, get the coeff of the volume form.
        Interpret as integration of c over the space"""
        output = sympy.Poly(c).coeff_monomial(self.volume)
        if not output:
            output = 0
        return output


class CompIntersection(object):
    """Complete intersection in a product of projective spaces"""

    def __init__(self, mp, degs):
        """To view a multiprojective space as a complete intersection, set degs = []
        In general, degs will be list of list, each having len=(# factors of the multiprojective space)
        """
        self.mp = mp
        self.degs = degs
        self.total_dim = mp.total_dim - len(degs)
        self._set_chern()

    def __str__(self):
        return f"CompIntersection in P^({self.mp.dims}) of deg={self.degs}"

    def _set_chern(self):
        # initially set chern class as that of ambient projective space
        self.chern = self.mp.chern
        self.chern_normal = 1
        # update by diving self.chern by chern classes of lines cutting out the variety if interest
        for deg in self.degs:
            self.chern *= self.mp.get_line_inv(deg)
            self.chern_normal *= self.mp.get_line(deg)
        self.chern = self.mp.truncate(self.chern)
        self.chern_normal = self.mp.truncate(self.chern_normal)
        self.chern_numbers = dict()

    def get_sub_intersection(self, deg):
        """Get another complete intersection determined by intersection with a hypersurface"""
        return CompIntersection(mp=self.mp, degs=self.degs + [deg])

    def get_c1_intersection(self):
        """Get the intersection giving by cutting by a line whose c_1 is the first chern class of self
        The result will necessarily have vanishing c1
        """
        deg = self.get_mth_chern(1).coeffs()
        return self.get_sub_intersection(deg)

    def get_product(self, other_ci):
        """Realize a product of complete intersection as a complete intersection"""
        out_dims = list(self.mp.dims) + list(other_ci.mp.dims)
        out_degs = [d + [0]*len(other_ci.mp.dims) for d in self.degs] + \
                   [[0]*len(self.mp.dims) + d for d in other_ci.degs]
        return CompIntersection(MultiProj(out_dims), out_degs)

    def get_mth_chern(self, m):
        """Get the mth Chern class, expressed as a class coming from the ambient projective space"""
        return self.mp.get_mth(self.chern, m)

    def get_chern_number(self, partition):
        """Get Chern number associated to a partition of self.total_dim"""
        total_degree = sum(partition)
        if total_degree != self.total_dim:
            raise ValueError(f"Variety of total_dim={self.total_dim} does not have {partition}th Chern number")
        if partition in self.chern_numbers:
            return self.chern_numbers[partition]
        chern_part = 1
        for p in partition:
            chern_part *= self.get_mth_chern(p)
        output = self.integrate(chern_part)
        self.chern_numbers[partition] = output
        return output

    def get_all_chern_numbers(self):
        """Get all Chern numbers"""
        all_parts = all_partitions(self.total_dim)
        for part in all_parts:
            self.get_chern_number(part)
        return self.chern_numbers

    def integrate(self, c):
        """Integrate a polynomial over the complete intersection"""
        big_c = c * self.chern_normal
        return self.mp.integrate(big_c)

    def get_chern_numbers_branched(self, branch_deg, branch_order, verbose=False):
        """Use Izawas formula to compute Chern numbers of the variety obtained by branching the complete intersection
        along a subvariety specified by a branch degreee, with specified branch order. Follows the article...

        `Chern number formula for ramified coverings` by Takeshi Izawa

        For the case when the branch locus is a smooth complete sub-intersection.

        When verbose=True, print the polynomial from the article"""
        if any([bd % branch_order != 0 for bd in branch_deg]):
            raise ValueError(f"Branch order {branch_order} does not divide the branch degree {branch_deg}")
        # sets the chern numbers if they're not computed
        self.get_all_chern_numbers()
        # first chern class of the line associated to branch locus
        branch_line_c1 = self.mp.get_mth(c=self.mp.get_line(branch_deg), m=1)
        branch_locus = self.get_sub_intersection(branch_deg)
        l = sympy.Symbol("l", commutative=True)
        mu = sympy.Symbol("m", commutative=True)
        # chern_symbols[i] = "c_i"
        chern_symbols = [sympy.Symbol(f"c{i}") for i in range(0, self.total_dim + 1)]
        output = dict()
        for part in self.chern_numbers.keys():
            # compute N_{i} indexed by c_{i}
            n = [0 for _ in range(self.total_dim+1)]
            for p in part:
                n[p] += 1
            # now n[i] = N_{i} in Izawa's notation
            # compute the product piece in H^(N)_{xi}
            big_prod_left = 1
            for i in range(1, self.total_dim+1):
                n_i = n[i]
                big_prod_left *= (chern_symbols[i] + l*chern_symbols[i-1])**n_i
            big_prod_right = big_prod_left.subs({l: 0})
            # H^{N1...Nn}_{xi} from Izawa's paper
            h = ((big_prod_left - big_prod_right)*(l**(-1))).expand()
            # get the coefficients of the polynomial P_{\alpha} determined by H as a list
            p_alpha = [h.coeff(l, k) for k in range(self.total_dim)]
            delta_poly = sum([(1 - mu**(i+1)) * mu**(-i) *p_alpha[i] * (l**i) for i in range(len(p_alpha))])
            if verbose:
                print("Izawa delta polynomial for partition", part, "\n", delta_poly)
            # Substitutions to get a cohomology class in the ambient multiprojective space
            delta = delta_poly.subs({mu: branch_order})
            delta = delta.subs({l: branch_line_c1.as_expr()})
            chern_subs = {chern_symbols[i]: branch_locus.get_mth_chern(i).as_expr() for i in range(self.total_dim)}
            delta = delta.subs(chern_subs)
            # Integrate Chern numbers along subvariety
            delta_number = branch_locus.integrate(delta)
            # Now add the contributions to Izawa's difference function
            c_part = branch_order * self.get_chern_number(part) + delta_number
            output[part] = c_part
        return output

def get_z_kernel(m, transpose = False):
    """For an integer matrix m provided as a list of lists, compute the kernel, returned as a list of lists
    We need this extra function because sympy computes kernels over the rationals

    The linear algebra here is borrowed from...
    https://stackoverflow.com/questions/53921654/fast-computation-of-integer-basis-for-kernel-of-a-matrix-using-gpu
    """
    mat = sympy.Matrix(m)
    if transpose:
        mat = mat.T
    rational_nullspace = mat.nullspace()
    output = []
    for v in rational_nullspace:
        lcm_denom = sympy.lcm([x.denominator for x in v])
        v_int = [int(lcm_denom * x) for x in v]
        assert sympy.gcd([x for x in v_int]) == 1
        output.append(v_int)
    return output


def get_additive_cob_gens(n):
    """Get a set of complete intersections which additively spans Omega^{U}_{2n}.
    This will contain lots of redundancies!
    We take all products of projective spaces and milnor hypersurfaces of dim > 2
    Recall the milnor hypersurface of dim = i + j - 1 is the deg=[1,1] hypersurface in (P^i)x(P^j)
    """
    output = [CompIntersection(MultiProj(p), []) for p in all_partitions(n)]
    if n <= 2:
        return output
    # get all of the milnor hypersurfaces of dim<=n
    # indexed by (i,j) with i<= j so that dim=i+j-1
    milnor_indices = dict()
    for milnor_dim in range(3,n+1):
        milnor_indices[milnor_dim] = set()
        i = 1
        while i <= milnor_dim:
            j = milnor_dim-i+1
            if i <= j:
                milnor_indices[milnor_dim].add((i,j))
            i += 1
    for total_milnor_dim in range(3,n+1):
        milnor_dim_parts = [p for p in all_partitions(total_milnor_dim) if all([x > 2 for x in p])]
        milnor_multi_dims = list()
        for p in milnor_dim_parts:
             milnor_multi_dims += itertools.product(*[milnor_indices[x] for x in p])
        milnor_multi_dims = list(set(milnor_multi_dims))
        proj_dims = all_partitions(n - total_milnor_dim)
        for md in milnor_multi_dims:
            for pd in proj_dims:
                mfld = CompIntersection(MultiProj([md[0][0],md[0][1]]), [[1,1]])
                for i in range(1, len(md)):
                    mfld = mfld.get_product(
                        CompIntersection(MultiProj([md[i][0],md[i][1]]), [[1,1]])
                    )
                if not total_milnor_dim == n:
                    projective_part = CompIntersection(MultiProj(pd),[])
                    mfld = mfld.get_product(projective_part)
                output.append(mfld)
    return output


def get_euler_only(n, su=False):
    """Find the GCD of c_n(X) over all X in Omega^{U}_{2n} such that all other Chern numbers are zero.
    If su==True, restrict to image of Omega^{SU}_{2n} in Omega^{U}_{2n} instead.
    In the latter case, we apply Theorem 5.11 of https://arxiv.org/pdf/1903.07178
    """
    assert n > 1
    parts = all_partitions(n)
    if not su:
        _EULER_INDEX = (n,)
        gens = get_additive_cob_gens(n)
        chern_numbers = {str(g): g.get_all_chern_numbers() for g in gens}
        # put these into a matrix with each row corresponding to a chern number other than c_n
        # each column corresponds to a product of projective spaces
        m = [[chern_numbers[str(g)][ind] for ind in parts if not ind==_EULER_INDEX] for g in gens]
        nullspace = get_z_kernel(m, transpose=True)
        eulers = [sum(n[i] * chern_numbers[str(gens[i])][_EULER_INDEX] for i in range(len(gens))) for n in nullspace]
    else:
        if 2*n % 8 != 4:
            # In this case, we look at W_n, the set of all cobordism classes in Omega^{U}_{2n}
            # all of whose chern numbers of the form c1*c1*... equal to zero
            # then Im(Omega^{SU}_{2n} in Omega^{U}_{2n}) is the subspace of W_n such that all classes of the form
            # c1*... are zero. So in this dim, Im(Omega^{SU}_{2n} in Omega^{U}_{2n}) is already contained in
            # the subspace of Omega^{U}_{2n} such that every chern number of the form c1*... is zero.
            return get_euler_only(n, su=False)
        if 2*n % 8 == 4:
            # In this case Im(Omega^{SU}_{2n} in Omega^{U}_{2n}) is the image of W_{n+1} under the image of the map
            # which sends a manifold X to Y the Poincare dual of its c1.
            # Since cI(Y) = (c1*cI)(X) we want to find all of the Y whose chern classes of the form
            # c1*cI = 0, except possibly c1*cn. Then the c1*cn(Y) will give us the desired Euler number
            # find the generators of Omega^{U}_{2n+2}
            _EULER_INDEX = (1,n)
            upper_gens = get_additive_cob_gens(n+1)
            chern_numbers = {str(g): g.get_all_chern_numbers() for g in upper_gens}
            upper_parts = chern_numbers[str(upper_gens[0])].keys()
            vanishing_parts = [c for c in upper_parts if c[0]==1 and not c == _EULER_INDEX]
            # find which linear combos of them have all chern numbers of the form c1*c1*... equal to zero
            m = [
                [chern_numbers[str(g)][ind] for ind in vanishing_parts]
                for g in upper_gens
            ]
            # the kernel will consist of those Y we are interested in
            nullspace = get_z_kernel(m, transpose=True)
            eulers = [
                sum(v[i] * chern_numbers[str(upper_gens[i])][_EULER_INDEX] for i in range(len(upper_gens)))
                for v in nullspace
            ]
    return sympy.gcd(eulers)
