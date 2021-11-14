# jackpy

<!-- badges: start -->
[![Documentation status](https://readthedocs.org/projects/jackpy/badge/)](http://jackpy.readthedocs.io)
<!-- badges: end -->

Jack polynomials with Python.

```
>>> from gmpy2 import mpq
>>> from sympy.combinatorics.partitions import IntegerPartition
>>> from jackpy.jack import JackPol
>>>
>>> poly = JackPol(3, IntegerPartition([2, 1]), alpha = mpq(3, 2))
>>> print(poly)
Poly(7/2*x_0**2*x_1 + 7/2*x_0**2*x_2 + 7/2*x_0*x_1**2 + 6*x_0*x_1*x_2
 + 7/2*x_0*x_2**2 + 7/2*x_1**2*x_2 + 7/2*x_1*x_2**2, x_0, x_1, x_2, domain='QQ')
```