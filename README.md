# MPCI: Multi-projective Complete Intersections

A python package for studying complete intersections in products of projective spaces.
The goal is to carry out some quick computations and search for interesting examples.
There are many existing packages which deal with varieties in a single projective space,
and looking at varieties in products yields more examples.

Currently, the package can compute Chern numbers of these varieties.
It can also compute Chern numbers of branched covers
when the branch locus is another complete intersection.

There is also some functionality for studying numerical properties
of weighted projective spaces and the complex cobordism ring.
See the examples below.

Users are warned that the computational results have not been rigorously tested,
only checked against examples which I've worked out by hand.
If you spot any errors, or would be interested in adding some functionality to the code base 
please contact me at my gmail address.

# Basic usage

To use the package, launch python from the root of this directory.
You'll need to have the `sympy` python package installed or use a virtual environment.
See the installation details below.
For documentation of the functions, you'll have to look at the code and follow these examples for now :)
```
$ python
>>> # Import the package
>>> import mpci
>>> # First make a P1 x P2 and a P2 and look at their total Chern classes.
>>> p21 = mpci.MultiProj([1,2])
>>> p12.chern
Poly(6*x_0*x_1**2 + 6*x_0*x_1 + 2*x_0 + 3*x_1**2 + 3*x_1 + 1, x_0, x_1, domain='ZZ')
>>> p2 = mpci.MultiProj([2])
>>> p2.chern
Poly(3*x_0**2 + 3*x_0 + 1, x_0, domain='ZZ')
>>> # Then get all of the Chern numbers of a divisor of degree (1,2)
>>> # This is a degree=5 del Pezzo, which is P2 blown up at 4 points.
>>> mpci.CompIntersection(p21, [[1,2]]).get_all_chern_numbers()
{(1, 1): 5, (2,): 7}
>>> # Here's a del Pezzo double plane we find by double branching P2 over a deg=4 hypersurface 
>>> mpci.CompIntersection(p2, []).get_chern_numbers_branched([4], 2)
{(1, 1): 2, (2,): 10}
>>> # Compute chern numbers of a divisor of degree (2,1) inside of P1 x P2
>>> mpci.CompIntersection(p12, [[2, 1]]).get_all_chern_numbers()
{(1, 1): 8, (2,): 4}
>>> # Compute chern numbers of P2 branched over a degree 2 curve. We'll get the same result as above.
>>> mpci.CompIntersection(p2, []).get_chern_numbers_branched([2], 2)
{(1, 1): 8, (2,): 4}
>>> # Here we make a product of three P1s and a P2 x P2
>>> p111 = mpci.MultiProj([1,1,1])
>>> p22 = mpci.MultiProj([2,2])
>>> # We find degree=6 del Pezzo surfaces in each
>>> mpci.CompIntersection(p111, [[1, 1, 1]]).get_all_chern_numbers()
{(1, 1): 6, (2,): 6}
>>> mpci.CompIntersection(p22, [[1, 1], [1,1]]).get_all_chern_numbers()
{(1, 1): 6, (2,): 6}
>>> # View the second as a subvariety of the deg=(1,1) hypersurface in P2 x P2.
>>> # Now let's increase the (multi)degrees of these complete intersections...
>>> # First multiply by 2, yielding K3 surfaces
>>> mpci.CompIntersection(p111, [[2,2,2]]).get_all_chern_numbers()
{(1, 1): 0, (2,): 24}
>>> mpci.CompIntersection(p22, [[1, 1], [2,2]]).get_all_chern_numbers()
{(1, 1): 0, (2,): 24}
>>> # Now we'll multiply by 3
>>> mpci.CompIntersection(p111, [[3,3,3]]).get_all_chern_numbers()
{(1, 1): 18, (2,): 90}
>>> mpci.CompIntersection(p22, [[1, 1], [3,3]]).get_all_chern_numbers()
{(1, 1): 18, (2,): 90}
```

Additional functionality computes the ring structure of weighted projective spaces 
following Kawasaki "Cohomology of twisted projective lens spaces". 
According to that paper, the additive cohomology has one generator `g_k` 
of degree 2k for each k. We warn that this code is neither highly optimized, nor tested.
```
>>> # create a dim=2 WPS with weight vector [1,1,3]
>>> wps = mpci.WPS([1,1,3]).coefs
[3]
```
Here the output means that `g_1 * g_1 = 3 * g_2`.
```
>>> mpci.WPS([1,2,3,4,5]).coefs
[30, 60, 60]
```
The output here means that `g_1 * g_1 = 30 * g_2`, `g_1 * g_2 = 60 * g_3`, and `g_1 * g_3 = 60 * g_4`.

# Installation

The project only requires `sympy`. To install it in a virtual environment...
```
$ python3 -m pipenv install Pipfile
```
You only need to run the above upon first installation and might want to change the python version in Pipfile
to ensure that this can run on your computer.
Then to run the virtual environment...
```
$ python3 -m pipenv shell
```
