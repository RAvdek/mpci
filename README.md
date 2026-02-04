# MPCI: Multi-projective Complete Intersections

A python package for studying complete intersections in products of projective spaces and complex cobordism rings.
The goal is to carry out some quick computations and search for interesting examples.
There are many existing packages which deal with varieties in a single projective space,
and looking at varieties in products yields more examples.

The package can compute Chern numbers of these varieties and describe varieties spanning `Omega^{U}_{2n}`
with redundancies. It can also compute Chern numbers of branched covers when the branch locus
is another complete intersection. There is also some functionality for studying numerical properties
of weighted projective spaces. See the examples below.

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

## Memory

Computing Chern numbers and partitions over and over can be quite expensive. 
The following functions will allow you to store your computations in a local database file and retrieve them in future
sessions.
```
>>> mpci.save_memory() # saves some current computations in a local file
>>> mpci.load_memory() # self-explanatory 
```
Note that these work with a single file! If you use `save_memory`, the information from your current session
will overwrite any previous computations. Likewise, `load_memory` will overwrite your current computations.

## Brieskorn-Pham Milnor fibers

Here we compute the dimension, rank of `H_2` (AKA Milnor number), Euler char, and signature 
of the `[2,3,5]` Brieskorn-Pham Milnor fiber, whose link is the Poincaré homology sphere:   
```
>>> f = mpci.BrieskornFiber([2,3,5])
>>> f.dim
2
>>> f.get_rank_middle_hom()
8
>>> f.get_chi()
9
>>> f.get_sigma()
-6
```
This functionality is new and needs to be tested against examples.

## Cobordism rings

Recall that multiprojective spaces generate `Omega^{U}_{*}` tensor Q and that milnor hypersurfaces
generate the complex cobordism ring `Omega^{U}_{*}`. Cf. this page on
[the manifold atlas][http://www.map.mpim-bonn.mpg.de/Complex_bordism#Milnor_hypersurfaces].

The function `mpci.cob_to_multiproj(mfld)` takes a list of complete intersections and expresses
each as a rational linear combination of multiprojectives. This computation...
```
>>> mpci.cob_to_multiproj([mpci.CompIntersection(mpci.MultiProj([1,2]),[[1,1]])])
[{(1, 1): 1, (2,): 0}]
```
shows that the index `[1,2]` milnor hypersurface (which is a blowup of `P2` at a point) is cobordant to `P1 x P1`. 
Next, we use the function to study Milnor hypersurfaces of complex dimension 3: 
```
>>> x = mpci.CompIntersection(mpci.MultiProj([1,3]),[[1,1]])
>>> y = mpci.CompIntersection(mpci.MultiProj([2,2]),[[1,1]])
>>> mpci.cob_to_multiproj([x,y])
[{(1, 2): 1, (3,): 0, (1, 1, 1): 0}, {(1, 2): 4, (3,): -3/2, (1, 1, 1): -3/2}]
```
The first line of the above output says that the `[1,3]` Milnor hypersurface is cobordant to `P^1 x P^2`
but that the `[2,2]` hypersurface cannot be expressed as an integral linear combination of multiprojectives. 
Next, we see that all of the milnor hypersurfaces of complex dimension 4 can be expresses
as integer linear combinations of projective spaces
```
>>> x = mpci.CompIntersection(mpci.MultiProj([1,4]),[[1,1]])
>>> y = mpci.CompIntersection(mpci.MultiProj([2,3]),[[1,1]])
>>> mpci.cob_to_multiproj([x,y])
[{(4,): 0, (1, 1, 2): 0, (2, 2): 0, (1, 3): 1, (1, 1, 1, 1): 0},
 {(4,): -2, (1, 1, 2): -6, (2, 2): 3, (1, 3): 4, (1, 1, 1, 1): 2}]
```
In particular, first line of the above output says that the `[1,4]` Milnor hypersurface is cobordant to `P^1 x P^3`.
Now let's look at the complex dimension 5 case:
```
>>> x = mpci.CompIntersection(mpci.MultiProj([1,5]),[[1,1]])
>>> y = mpci.CompIntersection(mpci.MultiProj([2,4]),[[1,1]])
>>> z = mpci.CompIntersection(mpci.MultiProj([3,3]),[[1,1]])
>>> mpci.cob_to_multiproj([x,y,z])
[{(5,): 0, (1, 1, 3): 0, (1, 4): 1, (2, 3): 0, (1, 2, 2): 0, (1, 1, 1, 2): 0, (1, 1, 1, 1, 1): 0},
 {(5,): -5/2, (1, 1, 3): -15/2, (1, 4): 5, (2, 3): 6, (1, 2, 2): -15/2, (1, 1, 1, 2): 10, (1, 1, 1, 1, 1): -5/2},
 {(5,): -10/3, (1, 1, 3): -12, (1, 4): 6, (2, 3): 11, (1, 2, 2): -16, (1, 1, 1, 2): 64/3, (1, 1, 1, 1, 1): -6}]
```
We see that the `[1,5]` milnor hypersurface is bordant to `P^1 x P^4` but that the `[2,4]` and `[3,3]`
milnor hypersufaces cannot be expressed as integer linear combinations of multiprojective spaces.

The function `milnors_cobordant_to_multiproj` enumerates the indices of Milnor hypersurfaces which are integrally 
cobordant to products of projective spaces:
```
>>> mpci.milnors_cobordant_to_multiproj(2)
[[1, 2]]
>>> mpci.milnors_cobordant_to_multiproj(3)
[[1, 3]]
>>> mpci.milnors_cobordant_to_multiproj(4)
[[1, 4], [2, 3]]
>>> mpci.milnors_cobordant_to_multiproj(5)
[[1, 5]]
```
The results are validated by the preceding computations with `cob_to_multiproj`.

The function `mpci.get_additive_cob_gens(n)` provides a list
of complete intersection varieties spanning `Omega^{U}_{2n}`, which in general will contain redundancies.
In practice, we have found that it is computationally expensive to compute chern numbers of products of milnor
hypersurfaces whose total dimension is large. This is probably due to sympy's polynomial methods not being 
computationally optimized. Algorithmically, it's also a bit tricky to list a minimal set of generators.
So we opt to have lists of generators with redundancies.
Below we see that `P1 x P1` and `P2` additively generate `Omega^{U}_{4}`:
```
>>> for g in mpci.get_additive_cob_gens(2):
...     print(g)
... 
CompIntersection in P^((1, 1)) of deg=[]
CompIntersection in P^((2,)) of deg=[]
```
Next we see that `P1 x P1`, `P3`, `P1 x P1 x P1`, and the `[2,2]` milnor hypersurface span `Omega^{U}_{6}`.
This is not a minimal set of generators!
```
>>> for g in mpci.get_additive_cob_gens(3):
...     print(g)
... 
CompIntersection in P^((1, 2)) of deg=[]
CompIntersection in P^((3,)) of deg=[]
CompIntersection in P^((1, 1, 1)) of deg=[]
CompIntersection in P^([2, 2]) of deg=[[1, 1]]
```
Here are the complex dim 4 and 5 cases:
```
>>> for g in mpci.get_additive_cob_gens(4):
...     print(g)
... 
CompIntersection in P^((4,)) of deg=[]
CompIntersection in P^((1, 1, 2)) of deg=[]
CompIntersection in P^((2, 2)) of deg=[]
CompIntersection in P^((1, 3)) of deg=[]
CompIntersection in P^((1, 1, 1, 1)) of deg=[]
CompIntersection in P^([2, 2, 1]) of deg=[[1, 1, 0]]
```
In the case n=4, we have `P^4`, the product `P^1 x P^1 x P2`, the product `P^2 x P^2`, the product `P^1 x P^3`, 
the product `(P^1)^4`, and the product of a (2,2) milnor 
hypersurface with `P^1`. Here is the complex dim 5 case:
```
>>> for g in mpci.get_additive_cob_gens(5):
...     print(g)
... 
CompIntersection in P^((5,)) of deg=[]
CompIntersection in P^((1, 1, 3)) of deg=[]
CompIntersection in P^((1, 4)) of deg=[]
CompIntersection in P^((2, 3)) of deg=[]
CompIntersection in P^((1, 2, 2)) of deg=[]
CompIntersection in P^((1, 1, 1, 2)) of deg=[]
CompIntersection in P^((1, 1, 1, 1, 1)) of deg=[]
CompIntersection in P^([2, 2, 1, 1]) of deg=[[1, 1, 0, 0]]
CompIntersection in P^([2, 2, 2]) of deg=[[1, 1, 0]]
CompIntersection in P^([2, 4]) of deg=[[1, 1]]
CompIntersection in P^([3, 3]) of deg=[[1, 1]]
```


The function `get_euler_only(n, su=False)` returns the GCD of the top chern number, `c_n(X)`, 
over all X in `Omega^{U}_{2n}` such that all other Chern numbers are zero.
If `su==True`, restrict to image of `Omega^{SU}_{2n}` in `Omega^{U}_{2n}` instead.
In the latter case, we apply Theorem 5.11 of [this article][https://arxiv.org/pdf/1903.07178] to carry out 
the computation. The output is independent of `su` output unless `2n = 4 mod 8`. See the cases `n=2,6` below. 
This computation is very heavy for `n > 7`, so logging is included to track progress of the computation. 
```
>>> mpci.get_euler_only(2)
-12
>>> mpci.get_euler_only(2, su=True) # Euler char of K3 surface
24
>>> mpci.get_euler_only(3, su=True) # Euler char of S6
2
>>> mpci.get_euler_only(4, su=True)
720
>>> mpci.get_euler_only(5, su=True)
24
>>> mpci.get_euler_only(6)
30240
>>> mpci.get_euler_only(6, su=True)
1391040
>>> mpci.get_euler_only(7, su=True)
1440
>>> mpci.get_euler_only(8, su=True) # currently 3 secs (down from 22 secs)
1209600
>>> mpci.get_euler_only(9, su=True) # this took 18 secs (down from 2:34, down from 4 minutes) on my laptop...
80640
>>> mpci.get_euler_only(10, su=True) # this took 7 minutes
95800320
>>> mpci.get_euler_only(10) # without the SU condition, we get half :) Took about 1:10
47900160
```
This code is not very well optimized :)

## Weighted projective spaces

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
