# jackpolynomials

<!-- badges: start -->
[![Documentation status](https://readthedocs.org/projects/jackpy/badge/)](http://jackpy.readthedocs.io)
<!-- badges: end -->

Jack polynomials with Python.

```python
>>> from gmpy2 import mpq
>>> from sympy.combinatorics.partitions import IntegerPartition
>>> from jackpy.jack import JackPol
>>>
>>> poly = JackPol(3, IntegerPartition([2, 1]), alpha = mpq(3, 2))
>>> print(poly)
Poly(7/2*x_1**2*x_2 + 7/2*x_1**2*x_3 + 7/2*x_1*x_2**2 + 6*x_1*x_2*x_3 
 + 7/2*x_1*x_3**2 + 7/2*x_2**2*x_3 + 7/2*x_2*x_3**2, x_1, x_2, x_3, domain='QQ')
```

As of the development version `0.1.0.9000`, it is possible to get Jack polynomials with a symbolic 
Jack parameter. These polynomials have domain `'QQ(alpha)'`, where `alpha` is
the Jack parameter.

```python
>>> from sympy.combinatorics.partitions import IntegerPartition
>>> from jackpy.jack import JackSymbolicPol
>>>
>>> poly = JackSymbolicPol(3, IntegerPartition([2, 2]))
>>> print(poly)
Poly((2*alpha**2 + 6*alpha + 4)*x_1**2*x_2**2 
 + (4*alpha + 8)*x_1**2*x_2*x_3 
 + (2*alpha**2 + 6*alpha + 4)*x_1**2*x_3**2 
 + (4*alpha + 8)*x_1*x_2**2*x_3 
 + (4*alpha + 8)*x_1*x_2*x_3**2 
 + (2*alpha**2 + 6*alpha + 4)*x_2**2*x_3**2, 
 x_1, x_2, x_3, domain='QQ(alpha)')
```
