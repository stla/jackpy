# jackpolynomials

<!-- badges: start -->
[![Documentation status](https://readthedocs.org/projects/jackpy/badge/)](http://jackpy.readthedocs.io)
<!-- badges: end -->

***Jack, zonal, and Schur polynomials with Python.***

Schur polynomials have applications in combinatorics and zonal polynomials have
applications in multivariate statistics. They are particular cases of
[Jack polynomials](https://en.wikipedia.org/wiki/Jack_function). This package
allows to get these polynomials.

___

The Jack polynomial in three variables of the integer partition $(2,1)$, with 
parameter $3/2$:

```python
>>> from gmpy2 import mpq
>>> from jackpy.jack import JackPol
>>>
>>> poly = JackPol(3, [2, 1], alpha = mpq(3, 2))
>>> print(poly)
Poly(7/2*x_1**2*x_2 + 7/2*x_1**2*x_3 + 7/2*x_1*x_2**2 + 6*x_1*x_2*x_3 
 + 7/2*x_1*x_3**2 + 7/2*x_2**2*x_3 + 7/2*x_2*x_3**2, x_1, x_2, x_3, domain='QQ')
```

By default, `JackPol` returns the Jack $J$-polynomial. It can also return 
the Jack $C$-polynomial, the Jack $P$-polynomial, and the Jack $Q$-polynomial,
by using the `which` argument.

The expression of a Jack polynomial can be long. Since a Jack polynomial is 
symmetric, it is possible to write it as a linear combination of some 
monomial symmetric polynomials, and this gives a shorter expression.

```python
>>> from jackpy.jack import JackPol
>>> from jackpy.monomial_symmetric_polynomials import msp_combination_expr
>>> poly = JackPol(4, [2, 2], alpha = 5)
>>> print(poly)
Poly(84*x_1**2*x_2**2 + 28*x_1**2*x_2*x_3 + 28*x_1**2*x_2*x_4
 + 84*x_1**2*x_3**2 + 28*x_1**2*x_3*x_4 + 84*x_1**2*x_4**2
 + 28*x_1*x_2**2*x_3 + 28*x_1*x_2**2*x_4 + 28*x_1*x_2*x_3**2
 + 24*x_1*x_2*x_3*x_4 + 28*x_1*x_2*x_4**2 + 28*x_1*x_3**2*x_4 
 + 28*x_1*x_3*x_4**2 + 84*x_2**2*x_3**2 + 28*x_2**2*x_3*x_4 
 + 84*x_2**2*x_4**2 + 28*x_2*x_3**2*x_4 + 28*x_2*x_3*x_4**2 
 + 84*x_3**2*x_4**2, x_1, x_2, x_3, x_4, domain='QQ')
>>> expr = msp_combination_expr(poly)
>>> print(expr)
24*M[1;1;1;1] + 28*M[2;1;1] + 84*M[2;2]
```

Here the symbol `M[1;1;1;1]` represents the monomial symmetric polynomial of 
the integer partition $(1,1,1,1)$.

As of the development version `0.1.0.9000`, it is possible to get Jack polynomials with a symbolic 
Jack parameter. These polynomials have domain `'QQ(alpha)'`, where `alpha` is
the Jack parameter.

```python
>>> from jackpy.jack import JackSymbolicPol
>>>
>>> poly = JackSymbolicPol(3, [2, 2])
>>> print(poly)
Poly((2*alpha**2 + 6*alpha + 4)*x_1**2*x_2**2 
 + (4*alpha + 8)*x_1**2*x_2*x_3 
 + (2*alpha**2 + 6*alpha + 4)*x_1**2*x_3**2 
 + (4*alpha + 8)*x_1*x_2**2*x_3 
 + (4*alpha + 8)*x_1*x_2*x_3**2 
 + (2*alpha**2 + 6*alpha + 4)*x_2**2*x_3**2, 
 x_1, x_2, x_3, domain='QQ(alpha)')
```
