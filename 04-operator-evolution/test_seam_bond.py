"""Exact theorem instances and counterexamples for the R3 hypotheses."""

from fractions import Fraction as Q
import unittest

from aghora_return import (
    add, compressed_return, conjugate, identity, inverse, matrix, mul,
    nilpotent_flow, return_operator, scale, sub, transpose,
)
from seam_bond import (
    AGHORA, apply, compatible_generator, metric_adjoint,
    projector_from_line_and_involution, rational_seam_return,
    seam_constraint, seam_geometry, weighted_norm2,
)


I = identity(2)
ZERO = scale(I, 0)


class SeamBondTests(unittest.TestCase):
    def test_typed_constraint_and_retained_removed_lines(self):
        for gamma in (Q(2), Q(-3, 5), Q(1, 7)):
            with self.subTest(gamma=gamma):
                data = seam_geometry(gamma)
                bond, u, v = data['bond'], data['seam_vector'], data['mirror_vector']
                self.assertEqual(seam_constraint(gamma, u), 0)
                self.assertEqual(seam_constraint(gamma, v), 2 * gamma)
                self.assertEqual(apply(bond, u), u)
                self.assertEqual(apply(bond, v), (0, 0))
                self.assertEqual(mul(bond, bond), bond)
                self.assertEqual(bond, matrix(((Q(1, 2), gamma / 2),
                                                (1 / (2 * gamma), Q(1, 2)))))

    def test_bare_seam_does_not_choose_complement(self):
        data = seam_geometry(2)
        other = matrix(((0, 2), (0, 1)))
        self.assertEqual(mul(other, other), other)
        self.assertEqual(apply(other, data['seam_vector']), data['seam_vector'])
        self.assertEqual(seam_constraint(2, apply(other, (5, 7))), 0)
        self.assertNotEqual(other, data['bond'])
        self.assertNotEqual(mul(mul(AGHORA, other), AGHORA), sub(I, other))

    def test_general_coordinate_covariance_and_line_normalization(self):
        data = seam_geometry(Q(-2, 3))
        change = matrix(((2, 1), (1, 1)))
        new_a, new_u = conjugate(AGHORA, change), apply(change, data['seam_vector'])
        expected = conjugate(data['bond'], change)
        self.assertEqual(projector_from_line_and_involution(new_a, new_u), expected)
        self.assertEqual(projector_from_line_and_involution(
            new_a, tuple(-7 * x for x in new_u)), expected)
        # Images of the two independent basis vectors determine B uniquely.
        self.assertEqual(mul(expected, conjugate(expected, new_a)), ZERO)

    def test_degenerate_seams_and_nonexact_inputs_are_rejected(self):
        for vector in ((0, 0), (0, 1), (1, 0)):
            with self.assertRaises(ValueError):
                projector_from_line_and_involution(AGHORA, vector)
        with self.assertRaises(ValueError):
            seam_geometry(0)
        with self.assertRaises(ValueError):
            projector_from_line_and_involution(matrix(((1, 1), (0, 1))), (2, 1))
        with self.assertRaises(TypeError):
            seam_geometry(0.5)
        with self.assertRaises(TypeError):
            rational_seam_return(2, 3, 0.1)

    def test_exchange_law_and_seam_reflection(self):
        bond = seam_geometry(Q(-7, 3))['bond']
        reflection = sub(scale(bond, 2), I)
        self.assertEqual(mul(mul(AGHORA, bond), AGHORA), sub(I, bond))
        self.assertEqual(mul(reflection, reflection), I)
        self.assertEqual(mul(mul(AGHORA, reflection), AGHORA), scale(reflection, -1))

    def test_commutator_constructs_complex_structure(self):
        data = seam_geometry(2)
        bond, structure = data['bond'], data['complex_structure']
        self.assertEqual(structure, matrix(((0, 2), (Q(-1, 2), 0))))
        self.assertEqual(structure, sub(mul(AGHORA, bond), mul(bond, AGHORA)))
        self.assertEqual(mul(structure, structure), scale(I, -1))
        self.assertEqual(apply(structure, data['seam_vector']), data['mirror_vector'])
        self.assertEqual(apply(structure, data['mirror_vector']),
                         tuple(-x for x in data['seam_vector']))

    def test_complex_operator_multiplication(self):
        structure = seam_geometry(Q(4, 9))['complex_structure']
        def embed(real, imag):
            return add(scale(I, real), scale(structure, imag))
        a, b, c, d = Q(2, 3), Q(-4, 7), Q(5, 2), Q(8, 9)
        self.assertEqual(mul(embed(a, b), embed(c, d)), embed(a*c - b*d, a*d + b*c))
        self.assertEqual(mul(embed(a, b), embed(a/(a*a+b*b), -b/(a*a+b*b))), I)

    def test_relative_metric_compatibility(self):
        data = seam_geometry(Q(-5, 2))
        h, bond, structure = data['metric'], data['bond'], data['complex_structure']
        self.assertEqual(metric_adjoint(AGHORA, h), AGHORA)
        self.assertEqual(metric_adjoint(bond, h), bond)
        self.assertEqual(metric_adjoint(structure, h), scale(structure, -1))
        self.assertEqual(mul(mul(transpose(data['basis']), h), data['basis']), scale(I, 2))
        self.assertEqual(metric_adjoint(bond, scale(h, 7)), bond)

    def test_diagonal_reference_covariance_including_metric(self):
        gamma, a, b = Q(-3, 5), Q(7, 2), Q(-4, 3)
        old, new = seam_geometry(gamma), seam_geometry(a*gamma/b)
        change = matrix(((a, 0), (0, b)))
        for name in ('bond', 'complex_structure'):
            self.assertEqual(new[name], conjugate(old[name], change))
        transported_metric = mul(mul(transpose(inverse(change)), old['metric']), inverse(change))
        self.assertEqual(new['metric'], scale(transported_metric, b*b))
        self.assertEqual(rational_seam_return(gamma, 3, Q(2, 9))['mu'],
                         rational_seam_return(a*gamma/b, 3, Q(2, 9))['mu'])

    def test_odd_compatible_generator_and_an_incompatible_one(self):
        data = seam_geometry(2)
        generator = compatible_generator(2, 3)
        self.assertEqual(generator, matrix(((0, 6), (Q(-3, 2), 0))))
        self.assertEqual(mul(mul(AGHORA, generator), AGHORA), scale(generator, -1))
        self.assertEqual(metric_adjoint(generator, data['metric']), scale(generator, -1))
        other = matrix(((0, 6), (1, 0)))
        self.assertEqual(mul(mul(AGHORA, other), AGHORA), scale(other, -1))
        self.assertNotEqual(metric_adjoint(other, data['metric']), scale(other, -1))

    def test_chosen_rational_example(self):
        data = rational_seam_return(2, 3, Q(1, 9))
        self.assertEqual(data['returned'], matrix(((Q(4, 5), Q(6, 5)),
                                                  (Q(3, 10), Q(-4, 5)))))
        seam_return = mul(mul(inverse(data['basis']), data['returned']), data['basis'])
        self.assertEqual(seam_return, matrix(((Q(3, 5), Q(4, 5)),
                                              (Q(4, 5), Q(-3, 5)))))
        self.assertEqual(data['compressed'], scale(data['bond'], Q(3, 5)))
        self.assertEqual(data['defect'], scale(data['bond'], Q(16, 25)))
        self.assertEqual(data['defect'], data['excursion'])

    def test_cayley_response_formula_for_exact_cases(self):
        for gamma, omega, parameter in ((2, 0, Q(5)), (-3, -2, Q(1, 5)),
                                         (Q(2, 7), 4, Q(-3, 8)), (1, 3, Q(2))):
            with self.subTest(gamma=gamma, omega=omega, parameter=parameter):
                data = rational_seam_return(gamma, omega, parameter)
                z = omega * parameter
                self.assertEqual(data['mu'], 2*z/(1+z*z))
                self.assertEqual(data['defect_coefficient'], (1-z*z)**2/(1+z*z)**2)
                self.assertEqual(mul(data['returned'], data['returned']), I)
                self.assertEqual(data['compressed'], scale(data['bond'], data['mu']))
                self.assertLessEqual(abs(data['mu']), 1)

    def test_weighted_norm_balance_and_positive_defect(self):
        data = rational_seam_return(2, 3, Q(1, 9))
        h, bond, returned, u = data['metric'], data['bond'], data['returned'], data['seam_vector']
        leakage = mul(mul(sub(I, bond), returned), bond)
        self.assertEqual(data['defect'], mul(metric_adjoint(leakage, h), leakage))
        norm = weighted_norm2(u, h)
        self.assertEqual(weighted_norm2(apply(data['compressed'], u), h)/norm, Q(9, 25))
        self.assertEqual(weighted_norm2(apply(leakage, u), h)/norm, Q(16, 25))
        for state in ((2, 1), (5, -7), (0, 3)):
            self.assertEqual(weighted_norm2(apply(returned, state), h), weighted_norm2(state, h))

    def test_both_return_signs_are_available_at_zero_leakage(self):
        for sign in (-1, 1):
            data = rational_seam_return(2, 3, Q(sign, 3))
            self.assertEqual(data['mu'], sign)
            self.assertEqual(data['defect'], ZERO)
            self.assertEqual(mul(data['returned'], data['bond']), scale(data['bond'], sign))
            self.assertEqual(data['returned'], scale(sub(scale(data['bond'], 2), I), sign))

    def test_free_rate_changes_observed_response(self):
        first = rational_seam_return(2, 1, Q(1, 3))
        second = rational_seam_return(2, 3, Q(1, 3))
        self.assertEqual(first['bond'], second['bond'])
        self.assertEqual(first['metric'], second['metric'])
        self.assertEqual((first['mu'], second['mu']), (Q(3, 5), Q(1)))

    def test_scalar_closure_can_hide_leakage_without_metric_compatibility(self):
        data = seam_geometry(1)
        generator = matrix(((0, 1), (0, 0)))
        returned = return_operator(AGHORA, generator, nilpotent_flow(generator, 2))
        compressed, defect, _ = compressed_return(returned, data['bond'])
        self.assertEqual(compressed, data['bond'])
        self.assertEqual(defect, ZERO)
        self.assertNotEqual(mul(returned, data['bond']), data['bond'])
        self.assertNotEqual(mul(mul(sub(I, data['bond']), returned), data['bond']), ZERO)
        self.assertNotEqual(metric_adjoint(generator, data['metric']), scale(generator, -1))

    def test_incompatible_dynamics_can_exceed_unit_response(self):
        data = seam_geometry(1)
        generator = matrix(((0, 1), (0, 0)))
        returned = return_operator(AGHORA, generator, nilpotent_flow(generator, 4))
        compressed, defect, _ = compressed_return(returned, data['bond'])
        self.assertEqual(compressed, scale(data['bond'], 2))
        self.assertEqual(defect, scale(data['bond'], -3))

    def test_retaining_opposite_seam_reverses_complex_orientation(self):
        data = seam_geometry(2)
        opposite = sub(I, data['bond'])
        opposite_k = sub(mul(AGHORA, opposite), mul(opposite, AGHORA))
        self.assertEqual(opposite, seam_geometry(-2)['bond'])
        self.assertEqual(opposite_k, scale(data['complex_structure'], -1))
        self.assertEqual(mul(opposite_k, opposite_k), scale(I, -1))


if __name__ == '__main__':
    unittest.main()
