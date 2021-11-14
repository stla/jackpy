# -*- coding: utf-8 -*-
from sympy.combinatorics.partitions import IntegerPartition
from jackpy.jack import ZonalPol
from sympy import symbols, Poly


def test_zonalpol():
    n = 4
    x_0, x_1, x_2, x_3 = symbols("x_0, x_1, x_2, x_3")
    expected = (
        Poly(x_0, x_0) + Poly(x_1, x_1) + Poly(x_2, x_2) + Poly(x_3, x_3) 
    )**(n - 1)
    obtained = (
        ZonalPol(n, IntegerPartition([3]))
        + ZonalPol(n, IntegerPartition([2,1]))
        + ZonalPol(n, IntegerPartition([1,1,1]))
    )
    assert obtained == expected.set_domain('QQ')

