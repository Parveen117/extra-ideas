"""Focused exact checks of R2 identities and their hypothesis boundaries."""

import unittest
from fractions import Fraction as Q

from aghora_return import (
    add, cayley_flow, compressed_return, conjugate, identity, inverse, matrix,
    mul, nilpotent_flow, orthogonal_example, parity_projectors,
    projection_for_scalar, return_operator, scale, sub, transpose,
    validate_source_relations,
)


class AghoraReturnTests(unittest.TestCase):
    def setUp(self):
        self.a = matrix(((1, 0), (0, -1)))
        self.g = matrix(((0, 1), (0, 0)))
        self.unit = identity(2)

    def returned(self, t):
        return return_operator(self.a, self.g, nilpotent_flow(self.g, t))

    def test_source_relations_allow_a_nonzero_generator(self):
        validate_source_relations(self.a, self.g)
        self.assertNotEqual(self.g, scale(self.unit, 0))
        self.assertEqual(mul(self.g, self.g), scale(self.unit, 0))

    def test_same_parameter_double_return_is_identity(self):
        for t in (Q(0), Q(1, 3), Q(-2), Q(7, 5)):
            with self.subTest(t=t):
                returned = self.returned(t)
                self.assertEqual(returned, matrix(((1, t), (0, -1))))
                self.assertEqual(mul(returned, returned), self.unit)

    def test_unequal_returns_leave_the_predicted_residual_flow(self):
        first, second = Q(5, 3), Q(2, 3)
        self.assertEqual(mul(self.returned(second), self.returned(first)),
                         nilpotent_flow(self.g, first - second))
        self.assertNotEqual(mul(self.returned(second), self.returned(first)), self.unit)

    def test_return_is_similar_to_aghora(self):
        t = Q(7, 3)
        change = nilpotent_flow(self.g, -t / 2)
        self.assertEqual(conjugate(self.a, change), self.returned(t))

    def test_parity_projectors_select_both_exact_signs(self):
        returned = self.returned(Q(2, 3))
        plus, minus = parity_projectors(returned)
        self.assertEqual(add(plus, minus), self.unit)
        self.assertEqual(mul(plus, minus), scale(self.unit, 0))
        self.assertEqual(mul(plus, plus), plus)
        self.assertEqual(mul(minus, minus), minus)
        self.assertEqual(mul(returned, plus), plus)
        self.assertEqual(mul(returned, minus), scale(minus, -1))
        self.assertNotEqual(plus, scale(self.unit, 0))
        self.assertNotEqual(minus, scale(self.unit, 0))

    def test_one_dimensional_source_relation_forces_zero_generator(self):
        for sign in (1, -1):
            validate_source_relations(matrix(((sign,),)), matrix(((0,),)))
            with self.assertRaises(ValueError):
                validate_source_relations(matrix(((sign,),)), matrix(((2,),)))

    def test_commuting_bond_does_not_choose_a_global_sign(self):
        returned = self.returned(1)
        compressed, defect, _ = compressed_return(returned, self.unit)
        self.assertEqual(compressed, returned)
        self.assertEqual(defect, scale(self.unit, 0))
        self.assertNotEqual(compressed, self.unit)
        self.assertNotEqual(compressed, scale(self.unit, -1))

    def test_arbitrary_algebraic_cut_coefficients(self):
        for value in (Q(-2), Q(-1), Q(0), Q(1, 3), Q(1), Q(3)):
            with self.subTest(value=value):
                bond = projection_for_scalar(value)
                self.assertEqual(mul(bond, bond), bond)
                compressed, defect, excursion = compressed_return(self.a, bond)
                self.assertEqual(compressed, scale(bond, value))
                self.assertEqual(defect, scale(bond, 1 - value * value))
                self.assertEqual(defect, excursion)

    def test_defect_identity_with_nontrivial_flow_and_bond(self):
        returned = self.returned(Q(5, 2))
        bond = projection_for_scalar(Q(2, 7))
        _, defect, excursion = compressed_return(returned, bond)
        self.assertEqual(defect, excursion)
        self.assertNotEqual(defect, scale(self.unit, 0))

    def test_similarity_transports_the_return_and_the_cut_together(self):
        flow = nilpotent_flow(self.g, Q(2, 3))
        returned = return_operator(self.a, self.g, flow)
        bond = projection_for_scalar(Q(1, 3))
        change = matrix(((2, 1), (1, 1)))
        transformed = return_operator(conjugate(self.a, change),
                                      conjugate(self.g, change),
                                      conjugate(flow, change))
        self.assertEqual(transformed, conjugate(returned, change))
        original_data = compressed_return(returned, bond)
        transformed_data = compressed_return(transformed, conjugate(bond, change))
        self.assertEqual(transformed_data, tuple(conjugate(x, change) for x in original_data))

    def test_orthogonal_example_has_continuous_compressed_response(self):
        _, _, flow, returned, bond = orthogonal_example()
        self.assertEqual(mul(transpose(flow), flow), self.unit)
        self.assertEqual(transpose(returned), returned)
        self.assertEqual(mul(returned, returned), self.unit)
        self.assertEqual(returned, matrix(((Q(4, 5), Q(3, 5)), (Q(3, 5), Q(-4, 5)))))
        compressed, defect, excursion = compressed_return(returned, bond)
        self.assertEqual(compressed, scale(bond, Q(4, 5)))
        self.assertEqual(defect, scale(bond, Q(9, 25)))
        self.assertEqual(defect, excursion)

    def test_positive_defect_is_the_omitted_channel_gram_matrix(self):
        _, _, _, returned, bond = orthogonal_example()
        omitted = mul(mul(sub(self.unit, bond), returned), bond)
        gram = mul(transpose(omitted), omitted)
        self.assertEqual(gram, compressed_return(returned, bond)[1])
        self.assertEqual(Q(4, 5) ** 2 + Q(3, 5) ** 2, Q(1))

    def test_orthogonal_bond_alone_does_not_imply_a_positive_defect(self):
        returned = self.returned(4)
        bond = scale(matrix(((1, 1), (1, 1))), Q(1, 2))
        self.assertEqual(transpose(bond), bond)
        self.assertNotEqual(transpose(returned), returned)
        compressed, defect, _ = compressed_return(returned, bond)
        self.assertEqual(compressed, scale(bond, 2))
        self.assertEqual(defect, scale(bond, -3))

    def test_generator_scale_remains_free_despite_fixed_return_signs(self):
        a, g, _, first, bond = orthogonal_example()
        scaled = scale(g, 2)
        second = return_operator(a, scaled, cayley_flow(scaled, Q(1, 3)))
        self.assertEqual(mul(first, first), self.unit)
        self.assertEqual(mul(second, second), self.unit)
        self.assertEqual(compressed_return(first, bond)[0], scale(bond, Q(4, 5)))
        self.assertEqual(compressed_return(second, bond)[0], scale(bond, Q(5, 13)))

    def test_an_unrelated_flow_is_not_licensed_by_aghora(self):
        with self.assertRaises(ValueError):
            return_operator(self.a, self.g, scale(self.unit, 2))
        with self.assertRaises(ValueError):
            validate_source_relations(self.unit, self.g)
        with self.assertRaises(ValueError):
            compressed_return(self.returned(1), scale(self.unit, 2))
        with self.assertRaises(ValueError):
            inverse(scale(self.unit, 0))
        with self.assertRaises(TypeError):
            matrix(((0.5,),))


if __name__ == '__main__':
    unittest.main()
