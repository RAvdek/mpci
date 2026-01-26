import math
from itertools import combinations
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
        subweights = list(combinations(self.weights, k+1))
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
        In general, degs will be list of list, each having len=(# factors of the multiprojective space)"""
        self.mp = mp
        self.degs = degs
        self.total_dim = mp.total_dim - len(degs)
        self.chern = mp.chern
        self.chern_normal = 1
        for deg in degs:
            self.chern *= self.mp.get_line_inv(deg)
            self.chern_normal *= self.mp.get_line(deg)
        self.chern = self.mp.truncate(self.chern)
        self.chern_numbers = dict()

    def get_sub_intersection(self, deg):
        """Get another complete intersection determined by intersection with a hypersurface"""
        return CompIntersection(mp=self.mp, degs=self.degs + [deg])

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


def get_euler_only(n):
    """Find the generator of Omega^{U}_{2n} such that all Chern numbers are zero except for c_{n}

    The output is a pair (int_generator, euler) where...
    - int_generator a dict {partition: coeff} describing a linear combination of products of projective space
    - euler is the euler characteristic of the generator

    The linear algebra here is borrowed from...
    https://stackoverflow.com/questions/53921654/fast-computation-of-integer-basis-for-kernel-of-a-matrix-using-gpu
    """
    _EULER_INDEX = (n,)
    assert n > 1
    parts = [p for p in all_partitions(n)]
    # compute all of the chern numbers of all products of projective spaces having total
    # dimension n. These give the generators of the complex bordism group Omega^{U}_{n}.
    chern_numbers = {p: CompIntersection(MultiProj(list(p)), []).get_all_chern_numbers() for p in parts}
    # put these into a matrix with each row corresponding to a chern number other than c_n
    # each column corresponds to a product of projective spaces
    m = [[chern_numbers[var][ind] for ind in parts if not ind==_EULER_INDEX] for var in parts]
    # the kernel will have rank 1
    nullspace = sympy.Matrix(m).T.nullspace()
    assert len(nullspace) == 1
    # get the generator as a list of sympy fractions
    rational_generator = list(nullspace[0].T)
    # get the lcm of the denomenators
    lcm_denom = sympy.lcm([x.denominator for x in rational_generator])
    int_generator = [int(lcm_denom * x) for x in rational_generator]
    # ensure that the gcd is 1
    assert sympy.gcd(int_generator) == 1
    # convert to a dictionary mapping partitions (indexing products of projective spaces) to coefficients
    int_generator = {parts[i]: int_generator[i] for i in range(len(parts))}
    # compute chern number of the generator
    special_cherns = {ind: sum([int_generator[p] * chern_numbers[p][ind] for p in parts]) for ind in parts}
    for p in parts:
        if not p == _EULER_INDEX:
            assert special_cherns[p] == 0
    euler = special_cherns[_EULER_INDEX]
    return int_generator, euler
