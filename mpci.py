import itertools
import logging
import math
import pickle
import sympy

LOG_FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
LOGGER = logging.getLogger(__name__)
_SESSION_CHERN_MEMORY = dict()
_SESSION_PARTITIONS = dict()
_DB_FILENAME = "./_MPCI_MEMORY.pk"

# These internal functions store computations of chern numbers

def _get_chern_from_db(dims, degs):
    return _SESSION_CHERN_MEMORY.get((tuple(dims), tuple([tuple(d) for d in degs])))


def _set_entry_in_db(dims, degs, chern_numbers=None):
    assert isinstance(chern_numbers, dict)
    if chern_numbers is not None:
        _SESSION_CHERN_MEMORY[(tuple(dims), tuple([tuple(d) for d in degs]))] = chern_numbers
    else:
        # create an empty dictionary otherwise
        _SESSION_CHERN_MEMORY[(tuple(dims), tuple([tuple(d) for d in degs]))] = dict()


def load_memory():
    try:
        with open(_DB_FILENAME, "rb") as f:
            memory = pickle.load(f)
        global _SESSION_CHERN_MEMORY
        global _SESSION_PARTITIONS
        _SESSION_CHERN_MEMORY = memory["chern"]
        _SESSION_PARTITIONS = memory["parts"]
    except OSError:
        LOGGER.warning("No database file found. Use save_memory() at the end of your session to create one.")


def save_memory():
    with open(_DB_FILENAME, "wb") as f:
        pickle.dump({"chern": _SESSION_CHERN_MEMORY, "parts": _SESSION_PARTITIONS}, f)


def e_2pii(num, denom):
    frac = sympy.Rational(num, denom) * 2 * sympy.pi
    return sympy.cos(frac) + sympy.I * sympy.sin(frac)


def get_bernoulli(n):
    """Return the nth Bernoulli number according to `topologists' conventions
    as described in Milnor-Stasheff Characteristic classes"""
    if n == 0:
        return 1
    return int(-((-1)**n)) * sympy.bernoulli(2 * n)


def l_series_func(n):
    if n == 0:
        return 1
    coeff = sympy.bernoulli(2 * n)
    coeff *= 2 ** (2 * n)
    coeff /= sympy.factorial(2 * n)
    return coeff


def todd_series_func(n):
    coeff = sympy.bernoulli(n)
    coeff /= sympy.factorial(n)
    return coeff


def l_genus_top_coeff(n):
    """Return the coefficient of p_{n} in L_{n}"""
    output = get_bernoulli(n)
    output *= 2**(2 * n)
    output *= 2**((2 * n) - 1) - 1
    output /= sympy.factorial(2 * n)
    return output

def all_partitions(n):
    """Returns all partitions of a positive integer n"""
    if n in _SESSION_PARTITIONS.keys():
        return _SESSION_PARTITIONS[n]
    answer = set()
    answer.add((n, ))
    for i in range(1, n):
        for j in all_partitions(n - i):
            answer.add(tuple(sorted((i, ) + j)))
    _SESSION_PARTITIONS[n] = answer
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
    # make a copy so we don't corrupt a mutable object
    array = list(array)[:]
    if len(array) == 0:
        return 1
    output = array.pop(0)
    while len(array) > 0:
        output *= array.pop(0)
    return output


class Genus(object):

    def __init__(self, series, n_vars):
        self.series = series
        self.n_vars = n_vars
        self._n_vars_names = [sympy.Symbol(f"x{i}") for i in range(1, self.n_vars + 1)]
        self._symmetric_poly = None
        self._symmetric_degs = None
        self._set_evaluation()
        self._set_symmetric_degrees()

    def _set_evaluation(self):
        # build polynomial f up to n variables
        poly = 1
        for v in self._n_vars_names:
            poly *= self.series.subs({'t': v})
        symmetrized = sympy.polys.polyfuncs.symmetrize(poly, formal=True)
        self._symmetric_poly = symmetrized[0].as_poly()

    def _set_symmetric_degrees(self):
        # this is very stupid and assumes that symmetric functions are called sn
        degs = dict()
        for g in self._symmetric_poly.gens:
            degs[g] = int(g.name.split('s')[1])
        self._symmetric_degs = degs

    def get_symmetric_poly(self, deg=None):
        poly = self._symmetric_poly
        if not deg:
            return poly
        if not deg > self.n_vars:
            # We get errors which I currently don't understand when deg=n_vars
            raise ValueError("Degree of symmetric poly must be less than degree")
        gens = poly.gens
        poly_dict = poly.as_dict()
        output = 0
        for k in poly_dict.keys():
            n_keys = len(k)
            d = sum([(i + 1)* k[i] for i in range(n_keys)])
            if d == deg:
                output += prod([gens[i] ** k[i] for i in range(n_keys)]) * poly_dict[k]
        return output

    def get_symmetric_degrees(self):
        return self._symmetric_degs

    @classmethod
    def from_series_function(cls, series_func, n_vars):
        t = sympy.Symbol("t")
        sf_zero = series_func(0)
        if sf_zero != 1:
            raise ValueError(f"series_func(0) = f{sf_zero} != 1")
        series = sum([series_func(k)*(t**k) for k in range(n_vars)])
        return cls(series, n_vars)

    @classmethod
    def L(cls, n_vars):
        """L genus from the Hirzebruch signature theorem"""
        return cls.from_series_function(l_series_func, n_vars)

    @classmethod
    def Todd(cls, n_vars):
        """The Todd genus"""
        return cls.from_series_function(todd_series_func, n_vars)


class BrieskornFiber(object):

    def __init__(self, coeffs):
        assert isinstance(coeffs, list)
        assert len(coeffs) > 1
        assert all([isinstance(x, int) for x in coeffs])
        assert all([x > 1 for x in coeffs])
        # store coeffs as tuple so its immutable
        self.coeffs = tuple(coeffs)
        self.dim = len(coeffs) - 1
        self._sigma = None
        self._monodromy_eigenvalues = None
        self._monodromy_polynomial = None
        self._monodromy_polynomial_var = sympy.Symbol("z")

    def get_rank_middle_hom(self):
        """Compute the rank of the middle dimensional homology"""
        return prod([x - 1 for x in self.coeffs])

    def get_chi(self):
        """Return Euler characteristic"""
        return 1 + (self.get_rank_middle_hom() * (-1)** self.dim)

    def get_monodromy_eigenvalues(self):
        """Return a list of monodromy eigenvalues (with repetition)"""
        # Return if already computed
        if self._monodromy_eigenvalues is not None:
            return self._monodromy_eigenvalues
        # Compute from https://webhomes.maths.ed.ac.uk/~v1ranick/slides/singexot.pdf
        # which yields the same determinant as https://www.numdam.org/item/SB_1966-1968__10__13_0.pdf
        # each eigenvalue is of the form
        # e^{2pi * i * k_{0} / c_{0}} * ... * e^{2pi * i * k_{n} / c_{n}} = e^{2pi * i * X}
        # x = (sum_{i=0}^{n} (k_{i} prod_{j \neq i}c_{j}))/ (prod_{0}^{n} c_{j})
        # where c_{j} = self.coeffs[j] and 0 < k_{j} < c_{j}
        # so we need to compute all possible X  = x * (prod c_{j})
        # start by enumerating the summands
        x_summands = {i: list() for i in range(len(self.coeffs))}
        for i in range(len(self.coeffs)):
            coeff = self.coeffs[i]
            prod_all_others = prod([self.coeffs[k] for k in range(len(self.coeffs)) if k != i])
            for j in range(1, coeff):
                x_summands[i].append( j * prod_all_others)
        big_xs = [0]
        indices = list(range(len(self.coeffs)))
        while len(indices) > 0:
            i = indices.pop()
            new_bigxs = [x + y for x in big_xs for y in x_summands[i]]
            big_xs = new_bigxs
        prod_all_coeffs = prod(self.coeffs)
        eigenvalues = [e_2pii(big_x, prod_all_coeffs) for big_x in big_xs]
        # a quick sanity check
        assert len(eigenvalues) == self.get_rank_middle_hom()
        self._monodromy_eigenvalues = eigenvalues
        return eigenvalues

    def get_monodromy_polynomial(self, val=None):
        """Compute determinant of the homological milnor monodromy"""
        if self._monodromy_polynomial is None:
            # compute if not already available
            poly = 1
            for v in self.get_monodromy_eigenvalues():
                poly *= (self._monodromy_polynomial_var - v)
            self._monodromy_polynomial = sympy.simplify(poly)
        if val is None:
            return self._monodromy_polynomial
        output = self._monodromy_polynomial.subs({self._monodromy_polynomial_var: val})
        return sympy.simplify(output)

    def boundary_is_homotopy_sphere(self):
        """Returns True or False using the fact that the boundary is a homotopy sphere
        iff det(1-monodromy) = \pm 1"""
        return self.get_monodromy_polynomial(1) in [-1,1]

    def get_kervaire(self):
        """Get the kervaire invariant if the complex dimension is odd"""
        if self.dim % 2 == 0:
            raise ValueError("Kervaire invariant only applicable in odd complex dimension")
        v = self.get_monodromy_polynomial(val=-1) % 8
        if v in [3,5]:
            return 1
        return 0

    def get_bp4m_boundary(self):
        """Returns x, m where
        - the boundary is x times the generator of bP_{4m},
        - m is the order of bP_{4m}
        Here 4m is the real dim of the Milnor fiber and bP_{4m} is the group of homotopy
        spheres bounding parallelizable manifolds"""
        if self.dim % 2 != 0:
            raise ValueError("Complex dimension not divisible by 2")
        if not self.boundary_is_homotopy_sphere():
            raise ValueError("Boundary is not a homotopy sphere")
        # compute order of bP_{4m}, using real dim / 4 = complex dim / 2
        v = get_bernoulli(self.dim / 2)
        v *= sympy.Rational(4,self.dim)
        v = v.numerator
        order_bp4m = v * (2**(self.dim - 2)) * (2**(self.dim - 1) - 1)
        assert self._sigma % 8 == 0
        return (self.get_sigma() // 8) % order_bp4m, order_bp4m

    def get_sigma(self):
        """Compute the signature as described in Hirzebruch's Bourbaki lecture"""
        if self._sigma is not None:
            return self._sigma
        if self.dim % 2 != 0:
            return 0
        # generate lists of rationals r of the form
        # sum_{0}^{n} a_{k}/coeff_{k}
        # with a_{k} = 1,...,coeff_{k}-1
        # Inspection on some examples indicates that the code looks good...
        # what should we get for [5,3,2,2,2]
        #
        # rationals come from sequences
        # (1,1,1,1,1), (1,2,1,1,1) -> 1/5 + 1/3 + 3/2, 1/5 + 2/3 + 3/2 -> (6 + 10 + 45)/30, (6 + 20 + 45)/30
        # (2,1,1,1,1), (2,2,1,1,1) -> 2/5 + 1/3 + 3/2, 2/5 + 2/3 + 3/2 -> (12 + 10 + 45)/30, (12 + 20 + 45)/30
        # (3,1,1,1,1), (3,2,1,1,1) -> 3/5 + 1/3 + 3/2, 3/5 + 2/3 + 3/2 -> (18 + 10 + 45)/30, (18 + 20 + 45)/30
        # (4,1,1,1,1), (4,2,1,1,1) -> 4/5 + 1/3 + 3/2, 4/5 + 2/3 + 3/2 -> (24 + 10 + 45)/30, (24 + 20 + 45)/30
        #
        # 61/30, 71/30 -> 2.03, 2.36
        # 67/30, 77/30 -> 2.23, 2.56
        # 73/30, 83/30 -> 2.43, 2.76
        # 79/30, 89/30 -> 2.63, 2.96
        coeffs = list(self.coeffs)[:]
        rationals = [0]
        while len(coeffs) > 0:
            coeff = coeffs.pop()
            new_rationals = list()
            for c in range(1,coeff):
                summand = sympy.Rational(c,coeff)
                new_rationals += [r + summand for r in rationals]
            rationals = new_rationals
        # now we need to understand how many of the rationals r have either...
        # (0 < r < 1 mod 2) or (1 < r < 2 mod 2)
        # 0 < a/b < 1 + 2k iff 0 < a < b + 2bk iff 0 < a < b mod 2b
        # 1 < a/b < 2k iff b < a < 2bk
        sigma_plus = 0
        sigma_minus = 0
        for r in rationals:
            top = r.numerator
            bottom = r.denominator
            if top % 2*bottom < bottom:
                sigma_plus += 1
            elif bottom < top:
                sigma_minus += 1
        self._sigma = sigma_plus - sigma_minus
        return self._sigma

    def get_gompf_boundary(self):
        """Compute the Generalized Gompf invariant of the boundary... To appear in some article
        This is tricky because I'm not sure that the sign of
        signature will match the orientation correctly!"""
        return 2*l_genus_top_coeff(self.dim)*self.get_chi() - self.get_sigma()


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
        """set total chern class and associated vars as that of ambient projective space"""
        self.chern = self.mp.chern
        self.chern_normal = 1
        # update by diving self.chern by chern classes of lines cutting out the variety if interest
        for deg in self.degs:
            self.chern *= self.mp.get_line_inv(deg)
            self.chern_normal *= self.mp.get_line(deg)
        self.chern = self.mp.truncate(self.chern)
        self.chern_normal = self.mp.truncate(self.chern_normal)

    def _has_key_in_memory(self, partition = None):
        """check the database _SESSION_CHERN_MEMORY to see if data is available"""
        db_data = _get_chern_from_db(self.mp.dims, self.degs)
        if db_data is None:
            # If there is no database entry, create it
            _set_entry_in_db(self.mp.dims, self.degs, dict())
            return False
        else:
            if partition is None:
                # we have a database entry and are not caring about finding a partition
                return True
        # now we're assuming we've found a DB entry and are looking for a particular partition
        if db_data.get(partition) is None:
            return False
        return True

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
        # do not compute if we've already done so!
        if self._has_key_in_memory(partition):
            return _get_chern_from_db(self.mp.dims, self.degs)[partition]
        chern_part = 1
        for p in partition:
            chern_part *= self.get_mth_chern(p)
        output = self.integrate(chern_part)
        return output

    def get_all_chern_numbers(self, log=False):
        """Get all Chern numbers"""
        if log:
            LOGGER.info(f"Get all Chern numbers for {self}")
        if self._has_key_in_memory():
            data = _get_chern_from_db(self.mp.dims, self.degs)
            if data is not None:
                # this is a bad pattern, but I need a quick fix
                if len(data.keys()) > 0:
                    return data
        all_parts = all_partitions(self.total_dim)
        all_chern_numbers = dict()
        # If self is P1 x Z with use Z.get_all_chern_numbers_times_P1()
        # This applies recursion and should make things much faster
        if len(self.mp.dims) > 1:
            # try to find a P1 factor we can split off
            p1_indices = [i for i in range(len(self.mp.dims)) if self.mp.dims[i]==1]
            if len(p1_indices) > 0:
                for i in p1_indices:
                    if all([d[i] == 0 for d in self.degs]):
                        new_dims = list(self.mp.dims[:])
                        new_degs = list(self.degs)
                        new_dims.pop(i)
                        for d in new_degs:
                            d.pop(i)
                        all_chern_numbers = CompIntersection(
                            MultiProj(new_dims), new_degs).get_all_chern_numbers_times_p1()
                        _set_entry_in_db(self.mp.dims, self.degs, all_chern_numbers)
                        return all_chern_numbers
        for part in all_parts:
            all_chern_numbers[part] = self.get_chern_number(part)
        _set_entry_in_db(self.mp.dims, self.degs, all_chern_numbers)
        return all_chern_numbers

    def get_all_chern_numbers_times_p1(self):
        """Find all of the self times P1. This gives an easy way to speed up computations
        The explicit formula we're leveraging is that for I = (I_1,...,I_k)

        c_I(P1 x M) = 2(c_{I_1 - 1,I2,..,I_K}(M) + ... + c_{I_1,I2,..,I_K-1}(M))
        """
        upper_parts = all_partitions(self.total_dim + 1)
        all_chern_numbers = self.get_all_chern_numbers()
        output = {p: 0 for p in upper_parts}
        for part in upper_parts:
            for i in range(len(part)):
                # copy the partition so it can be modified
                modified_part = list(part[:])
                modified_part[i] -= 1
                modified_part = tuple(sorted([x for x in modified_part if x > 0]))
                try:
                    output[part] += 2 * all_chern_numbers[tuple(sorted(modified_part))]
                except  KeyError as e:
                    LOGGER.info(f"Keyerror for get_all_chern_numbers_times_p1 of {self}")
                    raise e
        return output

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
        for part in self.get_all_chern_numbers().keys():
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
                LOGGER.info(f"Izawa delta polynomial for partition {part} \n {delta_poly}")
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


def cob_to_multiproj(mflds, log=False):
    """For a list of multiprojective spaces of the same dimension, express their rational bordism class
    as a collection of multiprojective spaces. Output is a list of dictionaries of the form
        { partition: coeff }
    """
    # ensure that we have a list of multiprojs all having the same dimensions
    assert isinstance(mflds, list)
    assert len(mflds) > 0
    assert all([isinstance(x,CompIntersection) for x in mflds])
    assert len(set([x.total_dim for x in mflds])) == 1
    mfld_cherns = [x.get_all_chern_numbers() for x in mflds]
    # this orders partitions, so that they can be used as an index
    parts = list(mfld_cherns[0].keys())
    mps = [CompIntersection(MultiProj(p), []) for p in parts]
    mp_cherns = [x.get_all_chern_numbers(log=log) for x in mps]
    # Put these into a matrix with...
    # the jth row corresponding to chern numbers indexed by jth partitions
    # the ith column corresponds to a product of projective spaces indexed by ith partition
    # We know that this matrix will be invertible over the rationals since Omega^{U}
    # tensor Q is a polynomial algebra on the P^i
    output = list()
    m = sympy.Matrix([[mpc[ind] for ind in parts] for mpc in mp_cherns]).T
    for mc in mfld_cherns:
        v = sympy.Matrix([mc[p] for p in parts])
        solution = list(m.LUsolve(v).T)
        solution_dict = {parts[i]: solution[i] for i in range(len(parts))}
        output.append(solution_dict)
    return output


def milnors_cobordant_to_multiproj(dim, log=False):
    """Generate a list of indices of milnor hypersurfaces of a given complex dim which are cobordant
    to products of projective spaces
    """
    milnors = list()
    i = 1
    while 2*i <= dim + 1:
        milnors.append(CompIntersection(MultiProj([i, dim + 1 - i]), [[1,1]]))
        i += 1
    if len(milnors) == 0:
        return list()
    coeffs = cob_to_multiproj(milnors, log=log)
    output = list()
    for i in range(len(milnors)):
        if all([v.denominator == 1 for v in coeffs[i].values()]):
            output.append(milnors[i].mp.dims)
    return output


def get_additive_cob_gens(n, log=False):
    """Get a set of complete intersections which additively spans Omega^{U}_{2n}.
    This will contain lots of redundancies!
    We take all products of projective spaces and milnor hypersurfaces of dim > 2
    Recall the milnor hypersurface of dim = i + j - 1 is the deg=[1,1] hypersurface in (P^i)x(P^j)
    """
    # want to generate a list of all possible multiplicative generators of dim < n
    # I'll just record the dims and degs as I don't need the actual objects yet
    gens=list()
    i = 1
    while i <= n:
        # add all of the projective spaces
        # d is tuple of dims, degs
        gens.append(CompIntersection(MultiProj([i]), []))
        # all of the milnor hypersurfaces except those cobordand to multiprojs
        to_exclude = milnors_cobordant_to_multiproj(i)
        if i >= 3:
            j = 1
            while 2*j <= i + 1:
                mfld = CompIntersection(MultiProj([j, i+1-j]), [[1,1]])
                if mfld not in to_exclude:
                    gens.append(mfld)
                j+=1
        i += 1
    output = list()
    under_construction = [{"last_index": i,"manifold": v} for i, v in enumerate(gens)]
    while len(under_construction) > 0:
        for d in under_construction:
            if d["manifold"].total_dim == n:
                output.append(d["manifold"])
                under_construction.remove(d)
            else:
                for i, v in enumerate(gens):
                    if d["last_index"] <= i and d["manifold"].total_dim + v.total_dim <= n:
                        under_construction.append({
                            "last_index": i,
                            "manifold": d["manifold"].get_product(v)
                        })
                under_construction.remove(d)
    output = sorted(output, key=lambda mfld: mfld.mp.dims)
    assert all([x.total_dim == n for x in output])
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
        LOGGER.info("Starting enumeration of cobordism generators")
        gens = get_additive_cob_gens(n, log=True)
        LOGGER.info(f"Found {len(gens)} generators. Starting chern number computation")
        chern_numbers = dict()
        for g in gens:
            chern_numbers[str(g)] = g.get_all_chern_numbers(log=True)
        LOGGER.info("Finish computing chern numbers. Starting matrix initialization")
        # put these into a matrix with each row corresponding to a chern number other than c_n
        # each column corresponds to a product of projective spaces
        m = []
        for g in gens:
            try:
                m.append([chern_numbers[str(g)][ind] for ind in parts if not ind==_EULER_INDEX])
            except KeyError as e:
                LOGGER.info(f"Keyerror for {g} with cher_numbers: \n {chern_numbers[str(g)]}")
                raise e
        LOGGER.info("Starting nullspace computation")
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
            LOGGER.info("Starting enumeration of cobordism generators")
            upper_gens = get_additive_cob_gens(n+1, log=True)
            LOGGER.info(f"Found {len(upper_gens)} generators. Starting chern number computation")
            chern_numbers = dict()
            for g in upper_gens:
                chern_numbers[str(g)] = g.get_all_chern_numbers(log=True)
            LOGGER.info("Finish computing chern numbers. Starting matrix initialization")
            upper_parts = chern_numbers[str(upper_gens[0])].keys()
            vanishing_parts = [c for c in upper_parts if c[0]==1 and not c == _EULER_INDEX]
            # find which linear combos of them have all chern numbers of the form c1*c1*... equal to zero
            m = [
                [chern_numbers[str(g)][ind] for ind in vanishing_parts]
                for g in upper_gens
            ]
            # the kernel will consist of those Y we are interested in
            LOGGER.info("Starting nullspace computation")
            nullspace = get_z_kernel(m, transpose=True)
            LOGGER.info(f"Found nullspace with {len(nullspace)} elements")
            eulers = [
                sum(v[i] * chern_numbers[str(upper_gens[i])][_EULER_INDEX] for i in range(len(upper_gens)))
                for v in nullspace
            ]
    return abs(sympy.gcd(eulers))
