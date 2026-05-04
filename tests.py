"""
tests.py — Unit tests for the mpci library.

Run with:
    python -m unittest tests.py
    python -m unittest tests.py -v          # verbose output
    python -m unittest tests.TestChernNumbers  # run one test class

Each test class includes mathematical explanations in docstrings so that
the expected values can be verified by a human consulting standard references.
"""

import unittest
import logging
import sympy

import mpci

# Suppress logging during tests
logging.disable(logging.INFO)


# =========================================================================
# Helpers: convenient constructors
# =========================================================================

def CPn(n):
    """CP^n as a CompIntersection (no defining equations in P^n)."""
    return mpci.CompIntersection(mpci.MultiProj([n]), [])


def product(*factors):
    """Product of CompIntersections, e.g. product(CPn(1), CPn(2))."""
    result = factors[0]
    for f in factors[1:]:
        result = result.get_product(f)
    return result


def hypersurface(ambient_dim, degree):
    """Degree-d smooth hypersurface in CP^{ambient_dim}.
    Has complex dimension ambient_dim - 1."""
    return mpci.CompIntersection(mpci.MultiProj([ambient_dim]), [[degree]])


def complete_intersection(ambient_dim, degrees):
    """Complete intersection of hypersurfaces of given degrees in CP^{ambient_dim}.
    Has complex dimension ambient_dim - len(degrees)."""
    return mpci.CompIntersection(mpci.MultiProj([ambient_dim]), [[d] for d in degrees])


def milnor_hypersurface(r, s):
    """The Milnor hypersurface H_{r,s}: the (1,1)-divisor in CP^r x CP^s.
    Has complex dimension r + s - 1."""
    return mpci.CompIntersection(mpci.MultiProj([r, s]), [[1, 1]])


# =========================================================================
# Layer 1: Chern numbers of specific varieties
# =========================================================================

class TestChernNumbersProjectiveSpaces(unittest.TestCase):
    """Test Chern numbers of projective spaces.

    MATHEMATICAL BACKGROUND:
    The tangent bundle of CP^n satisfies TCP^n + O = O(1)^{n+1}, so
    c(TCP^n) = (1+h)^{n+1} where h generates H^2(CP^n).
    This gives c_k(CP^n) = binom(n+1, k) * h^k.

    The top Chern number c_n[CP^n] = integral of c_n over [CP^n].
    Since c_n = binom(n+1, n) * h^n = (n+1) * h^n and integral of h^n
    over [CP^n] equals 1, we get c_n[CP^n] = n+1.

    This equals the Euler characteristic chi(CP^n) = n+1, which is a
    standard check (the Euler characteristic of CP^n has n+1 cells in
    its CW decomposition: one cell in each even dimension 0, 2, ..., 2n).
    """

    def test_cp1_euler_char(self):
        """chi(CP^1) = 2. (CP^1 = S^2, chi = 2-0+0 = 2.)"""
        cn = CPn(1).get_all_chern_numbers()
        self.assertEqual(cn[(1,)], 2)

    def test_cp2_euler_char(self):
        """chi(CP^2) = 3. Also c_{1,1}(CP^2) = binom(3,1)^2 = 9."""
        cn = CPn(2).get_all_chern_numbers()
        self.assertEqual(cn[(2,)], 3)
        self.assertEqual(cn[(1, 1)], 9)

    def test_cp3_euler_char(self):
        """chi(CP^3) = 4. c_3 = binom(4,3) = 4."""
        cn = CPn(3).get_all_chern_numbers()
        self.assertEqual(cn[(3,)], 4)

    def test_cp3_all_chern_numbers(self):
        """All Chern numbers of CP^3.

        c(CP^3) = (1+h)^4 = 1 + 4h + 6h^2 + 4h^3.
        c_1 = 4h, c_2 = 6h^2, c_3 = 4h^3.
        Chern numbers (integrate monomials of h-degree 3 over [CP^3]):
          c_{3}     = 4
          c_{1,2}   = 4 * 6 = 24
          c_{1,1,1} = 4^3 = 64
        """
        cn = CPn(3).get_all_chern_numbers()
        self.assertEqual(cn[(3,)], 4)
        self.assertEqual(cn[(1, 2)], 24)
        self.assertEqual(cn[(1, 1, 1)], 64)

    def test_cp4_euler_char(self):
        """chi(CP^4) = 5."""
        cn = CPn(4).get_all_chern_numbers()
        self.assertEqual(cn[(4,)], 5)


class TestChernNumbersProducts(unittest.TestCase):
    """Test Chern numbers of products of projective spaces.

    MATHEMATICAL BACKGROUND:
    For a product M x N, the total Chern class is c(M x N) = c(M) * c(N)
    by the Whitney sum formula (since T(M x N) = TM + TN).

    Key identity: c_k(M x N) = sum_{i+j=k} c_i(M) * c_j(N).

    For CP^1 x CP^1: c(CP^1 x CP^1) = (1+2a)(1+2b) where a,b generate
    H^2 of each factor. The fundamental class is [CP^1 x CP^1] with
    integral of a*b over it equal to 1.
      c_1 = 2a + 2b, c_2 = 4ab.
      c_{1,1} = (2a+2b)^2 = 4a^2 + 8ab + 4b^2. Only the ab term survives: 8.
      c_{2} = 4ab, integral = 4.
    The Euler characteristic is 4 = chi(CP^1) * chi(CP^1) = 2*2.
    """

    def test_p1_times_p1(self):
        cn = product(CPn(1), CPn(1)).get_all_chern_numbers()
        self.assertEqual(cn[(2,)], 4)
        self.assertEqual(cn[(1, 1)], 8)

    def test_p1_times_p2(self):
        """CP^1 x CP^2: chi = 2*3 = 6.

        c(CP^1 x CP^2) = (1+2a)(1+3b+3b^2) where a in H^2(CP^1), b in H^2(CP^2).
        Integrate over [CP^1 x CP^2] using int(a * b^2) = 1.
          c_3 = c_1(CP^1)*c_2(CP^2) = 2a * 3b^2 => integral = 6.
          c_{1,2} = (2a+3b)(3b^2+2a*3b) = ... easier to just compute.
        """
        cn = product(CPn(1), CPn(2)).get_all_chern_numbers()
        self.assertEqual(cn[(3,)], 6)

    def test_product_euler_chars(self):
        """chi(M x N) = chi(M) * chi(N) for all products.

        This is a consequence of the Kunneth theorem: the Euler characteristic
        is multiplicative under products.
        """
        for (a, b) in [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3)]:
            cn = product(CPn(a), CPn(b)).get_all_chern_numbers()
            top_chern = (a + b,)
            expected_chi = (a + 1) * (b + 1)
            self.assertEqual(cn[top_chern], expected_chi,
                             f"chi(CP^{a} x CP^{b}) should be {expected_chi}")


class TestChernNumbersHypersurfaces(unittest.TestCase):
    """Test Chern numbers of smooth hypersurfaces.

    MATHEMATICAL BACKGROUND:
    A smooth degree-d hypersurface V in CP^{n+1} has complex dimension n.
    Its normal bundle is O(d), so by adjunction:
        c(TV) = c(TCP^{n+1})|_V / c(O(d)) = (1+h)^{n+2} / (1+dh)
    restricted to V, where h is the hyperplane class.
    Integration over [V] gives: int_V h^n = d (the degree).

    The Euler characteristic of V is:
        chi(V) = int_V c_n(TV)
    """

    def test_quadric_surface(self):
        """Degree-2 surface in CP^3 (smooth quadric = CP^1 x CP^1).

        c(TV) = (1+h)^4 / (1+2h) restricted to V.
        (1+h)^4 = 1+4h+6h^2+4h^3, (1+2h)^{-1} = 1-2h+4h^2-8h^3 (mod h^4).
        c(TV) = (1+4h+6h^2+4h^3)(1-2h+4h^2-8h^3)
              = 1 + 2h + 2h^2 + ...
        So c_1 = 2h, c_2 = 2h^2.
        int_V h^2 = deg = 2.
        c_{2} = 2*2 = 4, c_{1,1} = (2h)^2 integrated = 4*2 = 8.
        Same as CP^1 x CP^1 (they are diffeomorphic!).
        """
        V = hypersurface(3, 2)
        cn = V.get_all_chern_numbers()
        self.assertEqual(cn[(2,)], 4)
        self.assertEqual(cn[(1, 1)], 8)

    def test_cubic_surface(self):
        """Degree-3 surface in CP^3 (del Pezzo of degree 3, blowup of CP^2 at 6 pts).

        chi = 9 (since chi(CP^2) = 3, blowing up adds 1 each time, 3 + 6 = 9).
        By Noether's formula: c_2 + c_{1,1} = 12*chi(O_V) = 12 for any surface,
        and for a del Pezzo of degree d, c_{1,1} = K^2 = d.
        So c_{1,1} = 3 (this is K^2 = degree), c_{2} = 9 (= chi).
        """
        V = hypersurface(3, 3)
        cn = V.get_all_chern_numbers()
        self.assertEqual(cn[(2,)], 9)
        self.assertEqual(cn[(1, 1)], 3)

    def test_quartic_surface_is_k3(self):
        """Degree-4 surface in CP^3 is a K3 surface.

        A K3 surface has c_1 = 0 (trivial canonical bundle), chi = 24.
        By adjunction: c_1(V) = c_1(CP^3) - [V] = 4h - 4h = 0. Confirmed.
        c_{1,1} = 0 (since c_1 = 0).
        c_{2} = chi = 24.
        """
        V = hypersurface(3, 4)
        cn = V.get_all_chern_numbers()
        self.assertEqual(cn[(1, 1)], 0)
        self.assertEqual(cn[(2,)], 24)

    def test_quintic_surface(self):
        """Degree-5 surface in CP^3 is a surface of general type.

        c_1(V) = (4-5)h = -h, so c_{1,1} = (-h)^2 int = 1*5 = 5.
        Noether: c_2 = 12*chi(O) - c_{1,1}. chi(O) = 1 + p_g - q.
        For a degree-5 surface: p_g = binom(5-1,3) = 4, q = 0, so chi(O) = 5.
        c_2 = 12*5 - 5 = 55.
        """
        V = hypersurface(3, 5)
        cn = V.get_all_chern_numbers()
        self.assertEqual(cn[(1, 1)], 5)
        self.assertEqual(cn[(2,)], 55)

    def test_noether_formula(self):
        """Noether's formula: (c_{1,1} + c_2) / 12 = chi(O_V) is a positive integer
        for any smooth projective surface.

        c_{1,1} + c_2 = 12 * chi(O_V).
        For degree-d surfaces in CP^3: chi(O) = binom(d-1, 3) + 1 = (d-1)(d-2)(d-3)/6 + 1.
        """
        for d in range(2, 8):
            V = hypersurface(3, d)
            cn = V.get_all_chern_numbers()
            noether = cn[(1, 1)] + cn[(2,)]
            self.assertEqual(noether % 12, 0,
                             f"Noether's formula fails for degree-{d} surface")
            chi_O = (d - 1) * (d - 2) * (d - 3) // 6 + 1
            self.assertEqual(noether // 12, chi_O,
                             f"chi(O) wrong for degree-{d} surface")


class TestChernNumbersK3(unittest.TestCase):
    """Test that various K3 surface constructions give the same Chern numbers.

    MATHEMATICAL BACKGROUND:
    A K3 surface is a simply connected compact complex surface with
    trivial canonical bundle. All K3 surfaces are diffeomorphic, so
    they have the same Chern numbers: c_{1,1} = 0, c_2 = 24.

    Several realizations as complete intersections:
      - Quartic in CP^3: degree [4]
      - Intersection of a quadric and a cubic in CP^4: degrees [2, 3]
      - Intersection of three quadrics in CP^5: degrees [2, 2, 2]
      - Double cover of CP^2 branched over a sextic (tested elsewhere)

    In each case, c_1 = 0 follows from adjunction:
    c_1(V) = (n+2 - sum(degrees)) * h, which equals 0 in all cases above.
    """

    def _check_k3(self, V, description):
        cn = V.get_all_chern_numbers()
        self.assertEqual(cn[(1, 1)], 0, f"c_{{1,1}} != 0 for {description}")
        self.assertEqual(cn[(2,)], 24, f"c_2 != 24 for {description}")

    def test_quartic_in_cp3(self):
        self._check_k3(hypersurface(3, 4), "quartic in CP^3")

    def test_quadric_cubic_in_cp4(self):
        self._check_k3(complete_intersection(4, [2, 3]),
                        "(2,3) CI in CP^4")

    def test_three_quadrics_in_cp5(self):
        self._check_k3(complete_intersection(5, [2, 2, 2]),
                        "(2,2,2) CI in CP^5")


# =========================================================================
# Layer 2: Product recursion formulas
# =========================================================================

class TestProductFormula(unittest.TestCase):
    """Test get_all_chern_numbers_of_product against direct computation.

    MATHEMATICAL BACKGROUND:
    For M^{2a} x N^{2b}, the Chern number c_I(M x N) with |I| = a + b is:

        c_I(M x N) = sum         c_{I_M}[M] * c_{I_N}[N]
                      splittings

    where the sum is over all ways to split each entry i_ell = p_ell + q_ell
    with sum(p_ell) = a and sum(q_ell) = b. The partition I_M is obtained
    from (p_1,...,p_r) by dropping zeros and sorting, likewise for I_N.

    We verify by comparing against the direct computation (building the product
    as a CompIntersection in the ambient product of projective spaces).
    """

    def _compare_product(self, M, N, label):
        """Compute Chern numbers of M x N two ways and compare."""
        direct = product(M, N).get_all_chern_numbers()
        via_formula = M.get_all_chern_numbers_of_product(N)
        self.assertEqual(direct, via_formula, f"Product formula mismatch for {label}")

    def test_cp1_times_cp1(self):
        self._compare_product(CPn(1), CPn(1), "CP^1 x CP^1")

    def test_cp1_times_cp2(self):
        self._compare_product(CPn(1), CPn(2), "CP^1 x CP^2")

    def test_cp2_times_cp2(self):
        self._compare_product(CPn(2), CPn(2), "CP^2 x CP^2")

    def test_cp1_times_cp3(self):
        self._compare_product(CPn(1), CPn(3), "CP^1 x CP^3")

    def test_cp2_times_cp3(self):
        self._compare_product(CPn(2), CPn(3), "CP^2 x CP^3")

    def test_cp3_times_cp3(self):
        self._compare_product(CPn(3), CPn(3), "CP^3 x CP^3")

    def test_cp1_times_h22(self):
        self._compare_product(CPn(1), milnor_hypersurface(2, 2), "CP^1 x H_{2,2}")

    def test_cp2_times_h22(self):
        self._compare_product(CPn(2), milnor_hypersurface(2, 2), "CP^2 x H_{2,2}")

    def test_cp3_times_h22(self):
        self._compare_product(CPn(3), milnor_hypersurface(2, 2), "CP^3 x H_{2,2}")

    def test_h22_times_h22(self):
        """H_{2,2} x H_{2,2}: a product of two Milnor hypersurfaces.

        This is the type of computation that previously required expensive
        polynomial arithmetic in an 4-fold multi-projective ambient space.
        """
        h22 = milnor_hypersurface(2, 2)
        self._compare_product(h22, h22, "H_{2,2} x H_{2,2}")

    def test_commutativity(self):
        """c_I(M x N) = c_I(N x M) for all I.

        The product formula is manifestly symmetric in M and N (just swap
        p_ell and q_ell in every splitting).
        """
        M, N = CPn(2), milnor_hypersurface(2, 2)
        cn_MN = M.get_all_chern_numbers_of_product(N)
        cn_NM = N.get_all_chern_numbers_of_product(M)
        self.assertEqual(cn_MN, cn_NM)


class TestProductSplitting(unittest.TestCase):
    """Test that _try_split_as_product correctly detects product structures.

    MATHEMATICAL BACKGROUND:
    A complete intersection in CP^{d_1} x ... x CP^{d_m} is a product X x Y
    if the ambient factors can be partitioned into two groups such that no
    defining equation involves factors from both groups. This is detected by
    finding connected components in the factor-equation incidence graph.
    """

    def test_bare_product_splits(self):
        """CP^2 x CP^3 (no equations) should split."""
        V = mpci.CompIntersection(mpci.MultiProj([2, 3]), [])
        split = V._try_split_as_product()
        self.assertIsNotNone(split)

    def test_milnor_does_not_split(self):
        """H_{2,2} = deg-(1,1) in CP^2 x CP^2 should NOT split.

        The equation links both ambient factors.
        """
        V = milnor_hypersurface(2, 2)
        split = V._try_split_as_product()
        self.assertIsNone(split)

    def test_single_factor_does_not_split(self):
        """CP^3 (single ambient factor) cannot split."""
        V = CPn(3)
        split = V._try_split_as_product()
        self.assertIsNone(split)

    def test_milnor_times_p1_splits(self):
        """CP^1 x H_{2,2} in CP^1 x CP^2 x CP^2 with degs [[0,1,1]] should split.

        The CP^1 factor (index 0) has degree 0 in the equation, so it's
        disconnected from the other factors.
        """
        V = mpci.CompIntersection(mpci.MultiProj([1, 2, 2]), [[0, 1, 1]])
        split = V._try_split_as_product()
        self.assertIsNotNone(split)

    def test_three_factor_product_splits(self):
        """CP^1 x CP^2 x CP^3 should split (into two groups on first call).

        The recursion will split the remaining group further.
        """
        V = mpci.CompIntersection(mpci.MultiProj([1, 2, 3]), [])
        split = V._try_split_as_product()
        self.assertIsNotNone(split)

    def test_split_gives_correct_chern_numbers(self):
        """Verify that splitting + product formula gives the same result
        as direct computation for CP^1 x H_{2,2}."""
        V = mpci.CompIntersection(mpci.MultiProj([1, 2, 2]), [[0, 1, 1]])
        # Via auto-splitting (get_all_chern_numbers will detect the split)
        cn_split = V.get_all_chern_numbers()
        # Via direct construction as a product
        cn_direct = product(CPn(1), milnor_hypersurface(2, 2)).get_all_chern_numbers()
        self.assertEqual(cn_split, cn_direct)

    def test_h22_times_h22_auto(self):
        """H_{2,2} x H_{2,2} constructed as CI in (CP^2)^4 should auto-split.

        This is the computation that was causing memory issues. With product
        splitting, it decomposes into two H_{2,2} factors.
        """
        V = mpci.CompIntersection(
            mpci.MultiProj([2, 2, 2, 2]),
            [[1, 1, 0, 0], [0, 0, 1, 1]])
        cn = V.get_all_chern_numbers()
        # Verify: chi(H_{2,2} x H_{2,2}) = chi(H_{2,2})^2
        h22 = milnor_hypersurface(2, 2)
        chi_h22 = h22.get_all_chern_numbers()[(3,)]
        self.assertEqual(cn[(6,)], chi_h22 ** 2)


# =========================================================================
# Layer 3: s-number and cobordism generators
# =========================================================================

class TestSNumber(unittest.TestCase):
    """Test the Milnor s-number computation.

    MATHEMATICAL BACKGROUND:
    The s-number s_n[M] is defined by expressing the n-th power sum
    symmetric function in terms of elementary symmetric functions (= Chern
    classes) via Newton's identities, then evaluating as a Chern number.

    Key values:
      s_n(CP^n) = n+1                    (standard result)
      s_n(H_{r,s}) = -binom(r+s, r)      (for the Milnor hypersurface)
      s_n(M x N) = 0  when dim M, dim N > 0  (s_n detects indecomposables)
    """

    def test_s_number_cpn(self):
        """s_n(CP^n) = n+1 for n = 1, ..., 6."""
        for n in range(1, 7):
            self.assertEqual(CPn(n).get_s_number(), n + 1)

    def test_s_number_milnor(self):
        """s_n(H_{r,s}) = -binom(r+s, r) for several Milnor hypersurfaces."""
        cases = [(2, 2), (2, 3), (2, 4), (3, 3), (4, 4)]
        for r, s in cases:
            h = milnor_hypersurface(r, s)
            expected = -int(sympy.binomial(r + s, r))
            self.assertEqual(h.get_s_number(), expected,
                             f"s-number wrong for H_{{{r},{s}}}")

    def test_s_number_products_vanish(self):
        """s_n(M x N) = 0 for products of positive-dimensional manifolds.

        This is because the power sum p_n(x_1,...,x_a, y_1,...,y_b)
        = p_n(x) + p_n(y), and when we evaluate on a product [M x N],
        the x-part lives in H^{2n}(M) with n > dim M, so it vanishes,
        and similarly for the y-part when n > dim N.
        More precisely: s_n of a product is 0 when n exceeds the dimension
        of each factor, which holds when n = dim(M) + dim(N) > max(dim M, dim N).
        """
        products_to_test = [
            (CPn(1), CPn(1)),
            (CPn(1), CPn(2)),
            (CPn(2), CPn(2)),
            (CPn(1), milnor_hypersurface(2, 2)),
        ]
        for M, N in products_to_test:
            P = product(M, N)
            self.assertEqual(P.get_s_number(), 0,
                             f"s-number should vanish on {M} x {N}")


class TestNewtonCoeffs(unittest.TestCase):
    """Test the Newton polynomial coefficients directly.

    MATHEMATICAL BACKGROUND:
    Newton's identities give:
      s_1 = c_1
      s_2 = c_1^2 - 2*c_2
      s_3 = c_1^3 - 3*c_1*c_2 + 3*c_3

    The coefficients dict maps partitions to integers. For example,
    for n=2: {(1,1): 1, (2,): -2} means s_2 = 1*c_{1,1} - 2*c_2.
    """

    def test_newton_n1(self):
        coeffs = mpci.get_newton_coeffs(1)
        self.assertEqual(coeffs, {(1,): 1})

    def test_newton_n2(self):
        coeffs = mpci.get_newton_coeffs(2)
        self.assertEqual(coeffs, {(1, 1): 1, (2,): -2})

    def test_newton_n3(self):
        coeffs = mpci.get_newton_coeffs(3)
        self.assertEqual(coeffs, {(1, 1, 1): 1, (1, 2): -3, (3,): 3})


class TestPolynomialGenerators(unittest.TestCase):
    """Test that polynomial generators satisfy Milnor's criterion.

    MATHEMATICAL BACKGROUND:
    Milnor's theorem: {x_1, x_2, ...} generates Omega^U_* as a polynomial
    ring over Z if and only if for each n:
      - s_n(x_n) = +/- p  when n+1 = p^k (a prime power)
      - s_n(x_n) = +/- 1  when n+1 is not a prime power

    The sign doesn't matter; only |s_n(x_n)| is constrained.
    """

    def test_milnor_criterion(self):
        """Check the Milnor criterion for n = 1 to 11."""
        for n in range(1, 12):
            gen = mpci.get_cobordism_polynomial_generator(n)
            s_val = sum(c * m.get_s_number() for c, m in gen)
            is_pp, p, k = mpci._is_prime_power(n + 1)
            if is_pp:
                self.assertEqual(abs(s_val), p,
                                 f"n={n}: |s_n| should be {p}, got {abs(s_val)}")
            else:
                self.assertEqual(abs(s_val), 1,
                                 f"n={n}: |s_n| should be 1, got {abs(s_val)}")

    def test_generator_dimensions(self):
        """Each polynomial generator should have the correct complex dimension."""
        for n in range(1, 12):
            gen = mpci.get_cobordism_polynomial_generator(n)
            for coeff, mfld in gen:
                self.assertEqual(mfld.total_dim, n,
                                 f"Generator for dim {n} contains manifold of dim {mfld.total_dim}")


# =========================================================================
# Layer 4: M_{2n} computation
# =========================================================================

class TestEulerOnly(unittest.TestCase):
    """Test the computation of M_{2n}.

    MATHEMATICAL BACKGROUND:
    M_{2n} is the GCD of c_{(n)}(X) over all X in Omega^U_{2n} satisfying
    c_I(X) = 0 for all partitions I of n with I != (n).

    Known values (verified by two independent code paths and hand checks):
      M_4  = 12   (n=2)
      M_6  = 2    (n=3)
      M_8  = 720  (n=4)
      M_10 = 24   (n=5)
      M_12 = 30240 (n=6)

    These fit the conjectured formulas:
      Even n: M_{2n} = denom(B_n / n!)  where B_n is the n-th Bernoulli number.
      Odd n >= 3: M_{2n} = (n-1)!
    """

    def test_M4(self):
        self.assertEqual(mpci.get_euler_only_v2(2), 12)

    def test_M6(self):
        self.assertEqual(mpci.get_euler_only_v2(3), 2)

    def test_M8(self):
        self.assertEqual(mpci.get_euler_only_v2(4), 720)

    def test_M10(self):
        self.assertEqual(mpci.get_euler_only_v2(5), 24)

    def test_M12(self):
        self.assertEqual(mpci.get_euler_only_v2(6), 30240)

    def test_even_formula(self):
        """For even n in {2, 4, 6, 8}: M_{2n} = denom(B_n / n!).

        B_2 = 1/6, B_4 = -1/30, B_6 = 1/42, B_8 = -1/30.
        denom(B_2/2!) = 12, denom(B_4/4!) = 720, denom(B_6/6!) = 30240.
        """
        for n in [2, 4, 6]:
            val = mpci.get_euler_only_v2(n)
            bn = sympy.bernoulli(n)
            expected = sympy.Rational(bn, sympy.factorial(n)).denominator
            self.assertEqual(val, expected,
                             f"Even formula fails for n={n}")

    def test_odd_formula(self):
        """For odd n in {3, 5}: M_{2n} = (n-1)!."""
        for n in [3, 5]:
            val = mpci.get_euler_only_v2(n)
            expected = int(sympy.factorial(n - 1))
            self.assertEqual(val, expected,
                             f"Odd formula fails for n={n}")


# =========================================================================
# Edge cases and consistency checks
# =========================================================================

class TestConsistency(unittest.TestCase):
    """Cross-checks and consistency tests."""

    def test_milnor_h12_cobordant_to_p1p1(self):
        """H_{1,2} (the (1,1)-divisor in CP^1 x CP^2) is cobordant to CP^1 x CP^1.

        They have the same Chern numbers, which means they represent the same
        class in Omega^U_4 (by the Hattori-Stong theorem: Chern numbers
        determine the cobordism class).
        """
        h12 = milnor_hypersurface(1, 2)
        p1p1 = product(CPn(1), CPn(1))
        self.assertEqual(h12.get_all_chern_numbers(), p1p1.get_all_chern_numbers())

    def test_hirzebruch_signature_formula(self):
        """For a 4-manifold: signature = (c_{1,1} - 2*c_2) / 3.

        This is Hirzebruch's signature theorem for complex surfaces:
        sigma(V) = (1/3)(c_1^2 - 2*c_2) = (c_{1,1} - 2*c_2) / 3.

        The signature must be an integer for any closed oriented 4-manifold.
        We check this for several surfaces.
        """
        surfaces = [
            CPn(2),
            product(CPn(1), CPn(1)),
            hypersurface(3, 2),
            hypersurface(3, 3),
            hypersurface(3, 4),
            hypersurface(3, 5),
        ]
        for V in surfaces:
            cn = V.get_all_chern_numbers()
            numerator = cn[(1, 1)] - 2 * cn[(2,)]
            self.assertEqual(numerator % 3, 0,
                             f"Hirzebruch formula fails for {V}: "
                             f"c_{{1,1}}={cn[(1,1)]}, c_2={cn[(2,)]}")

    def test_euler_char_of_hypersurfaces(self):
        """The Euler characteristic of a degree-d hypersurface in CP^{n+1} is:

        chi(V) = ((1-d)^{n+2} - 1) / d  +  (n+2)

        This follows from integrating the top Chern class using the
        adjunction formula c(TV) = (1+h)^{n+2}/(1+dh), with int_V h^n = d.

        We verify this for small cases. Note: this formula is equivalent to
        chi(V) = sum_{k=0}^{n} (-1)^k * binom(n+2, k+1) * d^k * (d-1)^0...
        actually the cleanest closed form is the one above.
        """
        def expected_chi(n, d):
            """chi of degree-d hypersurface in CP^{n+1}, complex dim n."""
            return ((1 - d) ** (n + 2) - 1) // d + (n + 2)

        for n in range(1, 5):
            for d in range(2, 6):
                V = hypersurface(n + 1, d)
                cn = V.get_all_chern_numbers()
                chi = cn[(n,)]
                self.assertEqual(chi, expected_chi(n, d),
                                 f"chi wrong for deg-{d} in CP^{n+1}")


if __name__ == '__main__':
    unittest.main()