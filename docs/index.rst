.. Jack polynomials documentation master file, created by
   sphinx-quickstart on Sun Nov 14 21:14:22 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Jack polynomials's documentation!
============================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. role:: raw-html(raw)
    :format: html

:raw-html:`<br />`
	
**jackpy.jack.JackPol(n, kappa, alpha)**

   Jack polynomial of an integer partition.

   :Parameters:
      *  **n** (*int*) – Positive integer, the number of variables of
         the polynomial.

      *  **kappa** (*IntegerPartition*) – An integer partition
         obtained with *sympy.combinatorics.partitions*.

      *  **alpha** (*number*) – A positive number, the parameter of
         the Jack polynomial.

   :Returns:
      The Jack polynomial of *kappa* with parameter *alpha*. The type
      of its coefficients depend on the type of *alpha*.

   :Return type:
      Poly

   -[ Examples ]-

   >>> from gmpy2 import mpq
   >>> from sympy.combinatorics.partitions import IntegerPartition
   >>> from jackpy.jack import JackPol
   >>>
   >>> poly = JackPol(3, IntegerPartition([2, 1]), alpha = mpq(3, 2))
   >>> print(poly)
   Poly(7/2*x_0**2*x_1 + 7/2*x_0**2*x_2 + 7/2*x_0*x_1**2 + 6*x_0*x_1*x_2
   + 7/2*x_0*x_2**2 + 7/2*x_1**2*x_2 + 7/2*x_1*x_2**2, x_0, x_1, x_2, domain='QQ')

**jackpy.jack.SchurPol(n, kappa)**

   Schur polynomial of an integer partition.

   :Parameters:
      *  **n** (*int*) – Positive integer, the number of variables of
         the polynomial.

      *  **kappa** (*IntegerPartition*) – An integer partition
         obtained with *sympy.combinatorics.partitions*.

   :Returns:
      The Schur polynomial of *kappa*. It has integer coefficents.

   :Return type:
      Poly

   -[ Examples ]-

   >>> from sympy.combinatorics.partitions import IntegerPartition
   >>> from jackpy.jack import SchurPol
   >>> p = SchurPol(2, IntegerPartition([2,1]))
   >>> print(p)
   Poly(x_1**2*x_0 + x_1*x_0**2, x_1, x_0, domain='ZZ')
   >>> y = p.eval({x_0: 1, x_1: 1})
   >>> print(y)
   2

**jackpy.jack.ZonalPol(m, kappa)**

   Zonal polynomial of an integer partition. Up to a normalization,
   this is the Jack polynomial of the integer partition with paramater
   *alpha = 2*.

   :Parameters:
      *  **m** (*int*) – Positive integer, the number of variables of
         the polynomial.

      *  **kappa** (*IntegerPartition*) – An integer partition
         obtained with *sympy.combinatorics.partitions*.

   :Returns:
      The zonal polynomial of *kappa*. It has rational coefficients.

   :Return type:
      Poly

**jackpy.jack.ZonalQPol(m, kappa)**

   Quaternionic zonal polynomial of an integer partition. Up to a
   normalization, this is the Jack polynomial of the integer partition
   with paramater *alpha = 1/2*.

   :Parameters:
      *  **m** (*int*) – Positive integer, the number of variables of
         the polynomial.

      *  **kappa** (*IntegerPartition*) – An integer partition
         obtained with *sympy.combinatorics.partitions*.

   :Returns:
      The quaternionic zonal polynomial of *kappa*. It has rational
      coefficients.

   :Return type:
      Poly

:raw-html:`<br />`
	
**jackpy.msf.msf_poly(m, kappa)**

   Monomial symmetric polynomial. This function is not highly related
   to the package; it has been introduced to facilitate some unit
   tests.

   :Parameters:
      *  **m** (*int*) â€“ Positive integer, the number of variables of
         the polynomial.

      *  **kappa** (*IntegerPartition*) â€“ An integer partition
         obtained with *sympy.combinatorics.partitions*.

   :Returns:
      The monomial symmetric polynomial corresponding to *kappa*.

   :Return type:
      Poly
