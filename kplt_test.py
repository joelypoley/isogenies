from __future__ import print_function

import unittest

from sage.all import *
from kplt import prime_norm_representative
from kplt import element_of_norm


class IdealsTest(unittest.TestCase):

    def test_prime_norm_representative(self):
        B = QuaternionAlgebra(59)
        max_order = B.maximal_order()
        I = max_order.left_ideal([2 * a for a in max_order.basis()])
        J = prime_norm_representative(I, max_order)
        self.assertTrue(J.norm() in Primes() and J.left_order() == max_order)

    def test_element_of_norm(self):
        B = QuaternionAlgebra(59)
        max_order = B.maximal_order()
        I = max_order.left_ideal([2 * a for a in max_order.basis()])
        M = 200001
        gamma = element_of_norm(M, max_order)
        self.assertTrue(gamma.reduced_norm() == M)

    def test_element_of_norm_large(self):
        B = QuaternionAlgebra(1019)
        max_order = B.maximal_order()
        I = max_order.left_ideal([2 * a for a in max_order.basis()])
        M = 200000
        gamma = element_of_norm(M, max_order)
        self.assertTrue(gamma.reduced_norm() == M)


if __name__ == '__main__':
    unittest.main()