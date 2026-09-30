"""R7: typed finite tensors and transport on admitted EMK carriers.

The scalar field of this executable sector is Q. KIR is an operator
algebra, and integer sheets are a separate path ledger, not scalar entries.
No manifold, clock, metric, or Hilbert space is required by the core.
General proofs and the smooth-sector extension are in EMK_TENSOR_CALCULUS_R7.md.
"""

from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import permutations, product
from math import factorial, prod

from aghora_return import (
    add, identity, inverse, matrix, mul, scalar, scale, sub, transpose,
)


I = identity(2)
R = matrix(((0, -1), (1, 0)))
K = matrix(((0, 1), (1, 0)))
RK = mul(R, K)
EPSILON = scale(R, -1)


def commutator(a, b):
    return sub(mul(a, b), mul(b, a))


@dataclass(frozen=True)
class Slot:
    fiber: str
    variance: int
    dimension: int = 2

    def __post_init__(self):
        if not isinstance(self.fiber, str) or not self.fiber:
            raise ValueError('Each slot needs a named carrier')
        if type(self.variance) is not int or self.variance not in (-1, 1):
            raise ValueError('Variance is +1 for a vector, -1 for a dual')
        if type(self.dimension) is not int or self.dimension < 1:
            raise ValueError('Carrier dimension must be a positive integer')


@dataclass(frozen=True)
class Tensor:
    slots: tuple
    data: tuple

    def __post_init__(self):
        slots = tuple(self.slots)
        if not all(isinstance(slot, Slot) for slot in slots):
            raise TypeError('Tensor slots must be Slot objects')
        dimensions = {}
        for slot in slots:
            if slot.fiber in dimensions and dimensions[slot.fiber] != slot.dimension:
                raise ValueError('A named carrier has one dimension')
            dimensions[slot.fiber] = slot.dimension
        data = tuple(scalar(x) for x in self.data)
        if len(data) != prod(slot.dimension for slot in slots):
            raise ValueError('Component count differs from the tensor type')
        object.__setattr__(self, 'slots', slots)
        object.__setattr__(self, 'data', data)

    def indices(self):
        return product(*(range(slot.dimension) for slot in self.slots))

    def at(self, indices):
        indices = tuple(indices)
        if len(indices) != len(self.slots):
            raise ValueError('Wrong number of indices')
        offset = 0
        for index, slot in zip(indices, self.slots):
            if type(index) is not int or not 0 <= index < slot.dimension:
                raise IndexError('Tensor index out of range')
            offset = offset * slot.dimension + index
        return self.data[offset]


def tensor_add(a, b):
    if a.slots != b.slots:
        raise ValueError('Only tensors of the same type can be added')
    return Tensor(a.slots, tuple(x + y for x, y in zip(a.data, b.data)))


def tensor_scale(t, factor):
    factor = scalar(factor)
    return Tensor(t.slots, tuple(factor * x for x in t.data))


def tensor_sub(a, b):
    return tensor_add(a, tensor_scale(b, -1))


def tensor_product(a, b):
    return Tensor(a.slots + b.slots, tuple(x * y for x in a.data for y in b.data))


def tensor_from_matrix(value, row_slot, column_slot):
    value = matrix(value)
    if row_slot.dimension != len(value) or column_slot.dimension != len(value):
        raise ValueError('Matrix dimensions differ from tensor slots')
    return Tensor((row_slot, column_slot), tuple(x for row in value for x in row))


def tensor_matrix(t):
    if len(t.slots) != 2 or t.slots[0].dimension != t.slots[1].dimension:
        raise ValueError('A square rank-two tensor is required')
    n = t.slots[0].dimension
    return matrix([t.data[i*n:(i+1)*n] for i in range(n)])


def apply_slot(t, axis, operator):
    """Apply a linear map to one slot; projection maps need not be invertible."""
    if type(axis) is not int or not 0 <= axis < len(t.slots):
        raise IndexError('Tensor axis out of range')
    operator = matrix(operator)
    n = t.slots[axis].dimension
    if len(operator) != n:
        raise ValueError('Operator dimension differs from its slot')
    values = []
    for indices in t.indices():
        values.append(sum((operator[indices[axis]][j] * t.at(
            indices[:axis] + (j,) + indices[axis+1:]) for j in range(n)), Q(0)))
    return Tensor(t.slots, tuple(values))


def transport(t, operators):
    """rho(U): vectors use U; duals use U^-T, independently by carrier."""
    out = t
    used = {}
    for axis, slot in enumerate(t.slots):
        key = (slot.fiber, slot.variance)
        if key not in used:
            if slot.fiber not in operators:
                raise ValueError('Missing transport for carrier ' + slot.fiber)
            value = matrix(operators[slot.fiber])
            # Require an isomorphism even if this tensor has only upper slots.
            inv = inverse(value)
            used[key] = value if slot.variance == 1 else transpose(inv)
        out = apply_slot(out, axis, used[key])
    return out


def generator_action(t, generators):
    """d rho(A): sum A on vector slots and -A^T on dual slots."""
    out = tensor_scale(t, 0)
    for axis, slot in enumerate(t.slots):
        if slot.fiber not in generators:
            raise ValueError('Missing generator for carrier ' + slot.fiber)
        a = matrix(generators[slot.fiber])
        action = a if slot.variance == 1 else scale(transpose(a), -1)
        out = tensor_add(out, apply_slot(t, axis, action))
    return out


def induced_matrix(slots, operators, *, infinitesimal=False):
    slots = tuple(slots)
    size = prod(s.dimension for s in slots)
    columns = []
    action = generator_action if infinitesimal else transport
    for column in range(size):
        basis = Tensor(slots, tuple(Q(i == column) for i in range(size)))
        columns.append(action(basis, operators).data)
    return matrix(zip(*columns))


def contract(t, first, second):
    if (type(first) is not int or type(second) is not int or first == second
            or not 0 <= first < len(t.slots) or not 0 <= second < len(t.slots)):
        raise IndexError('Two distinct tensor axes are required')
    a, b = t.slots[first], t.slots[second]
    if a.fiber != b.fiber or a.dimension != b.dimension or a.variance != -b.variance:
        raise ValueError('Contraction needs a carrier and its own dual')
    keep = tuple(i for i in range(len(t.slots)) if i not in (first, second))
    slots = tuple(t.slots[i] for i in keep)
    values = []
    for reduced in product(*(range(s.dimension) for s in slots)):
        indices = [0] * len(t.slots)
        for axis, index in zip(keep, reduced):
            indices[axis] = index
        total = Q(0)
        for j in range(a.dimension):
            indices[first] = indices[second] = j
            total += t.at(indices)
        values.append(total)
    return Tensor(slots, tuple(values))


def permute(t, order):
    order = tuple(order)
    if any(type(i) is not int for i in order) or sorted(order) != list(range(len(t.slots))):
        raise ValueError('A permutation of all tensor axes is required')
    slots = tuple(t.slots[i] for i in order)
    values = []
    for indices in product(*(range(s.dimension) for s in slots)):
        old = [0] * len(order)
        for new_axis, old_axis in enumerate(order):
            old[old_axis] = indices[new_axis]
        values.append(t.at(old))
    return Tensor(slots, tuple(values))


def symmetrize(t, axes=None, *, alternating=False):
    axes = tuple(range(len(t.slots))) if axes is None else tuple(axes)
    if (any(type(i) is not int or not 0 <= i < len(t.slots) for i in axes)
            or len(set(axes)) != len(axes)):
        raise IndexError('Symmetry axes must be distinct valid axes')
    if axes and any(t.slots[i] != t.slots[axes[0]] for i in axes):
        raise ValueError('Symmetrization needs identical slot types')
    out = tensor_scale(t, 0)
    for rearrangement in permutations(range(len(axes))):
        order = list(range(len(t.slots)))
        for position, old_position in enumerate(rearrangement):
            order[axes[position]] = axes[old_position]
        parity = sum(rearrangement[i] > rearrangement[j]
                     for i in range(len(axes)) for j in range(i+1, len(axes)))
        sign = -1 if alternating and parity % 2 else 1
        out = tensor_add(out, tensor_scale(permute(t, order), sign))
    return tensor_scale(out, Q(1, factorial(len(axes))))


def wedge(a, b):
    slots = a.slots + b.slots
    if slots and (slots[0].variance != -1 or any(s != slots[0] for s in slots)):
        raise ValueError('Forms must use one identical dual carrier')
    if a != symmetrize(a, alternating=True) or b != symmetrize(b, alternating=True):
        raise ValueError('Inputs to wedge must be alternating forms')
    factor = Q(factorial(len(slots)), factorial(len(a.slots))*factorial(len(b.slots)))
    return tensor_scale(symmetrize(tensor_product(a, b), alternating=True), factor)


def change_variance(t, axis, metric):
    if type(axis) is not int or not 0 <= axis < len(t.slots):
        raise IndexError('Tensor axis out of range')
    metric = matrix(metric)
    if transpose(metric) != metric:
        raise ValueError('Index raising/lowering requires a symmetric metric')
    inv = inverse(metric)
    old = t.slots[axis]
    out = apply_slot(t, axis, metric if old.variance == 1 else inv)
    slots = list(t.slots)
    slots[axis] = Slot(old.fiber, -old.variance, old.dimension)
    return Tensor(tuple(slots), out.data)


def increment(source, target, operators, ledger=None):
    """Backward transported difference rho(U)T_source - T_target - L."""
    out = tensor_sub(transport(source, operators), target)
    return out if ledger is None else tensor_sub(out, ledger)


def covariant_derivative(value, partial, generators):
    return tensor_add(partial, generator_action(value, generators))


def curvature(a_i, a_j, partial_i_a_j=None, partial_j_a_i=None,
              bracket_connection=None):
    """F_ij = d_i A_j - d_j A_i + [A_i,A_j] - A_[i,j].

    Derivative values must be supplied for varying coefficients. Omitting
    them declares constant coefficients; it does not estimate a derivative.
    """
    out = commutator(matrix(a_i), matrix(a_j))
    if partial_i_a_j is not None:
        out = add(out, matrix(partial_i_a_j))
    if partial_j_a_i is not None:
        out = sub(out, matrix(partial_j_a_i))
    if bracket_connection is not None:
        out = sub(out, matrix(bracket_connection))
    return out


def triangle_curvature(u_01, u_12, u_02):
    """F_012: V_0 -> V_2, direct minus ordered two-edge transport."""
    return sub(matrix(u_02), mul(matrix(u_12), matrix(u_01)))


def tetrahedron_bianchi(edges):
    """F_123 U_01 - F_023 + F_013 - U_23 F_012, exactly zero."""
    u01, u12, u23 = (matrix(edges[i]) for i in ('01', '12', '23'))
    u02, u13, u03 = (matrix(edges[i]) for i in ('02', '13', '03'))
    f012 = triangle_curvature(u01, u12, u02)
    f123 = triangle_curvature(u12, u23, u13)
    f023 = triangle_curvature(u02, u23, u03)
    f013 = triangle_curvature(u01, u13, u03)
    return sub(add(sub(mul(f123, u01), f023), f013), mul(u23, f012))


def metric_defect(u, source_metric, target_metric):
    """Pull target metric back: U^T H_target U - H_source."""
    u = matrix(u)
    return sub(mul(mul(transpose(u), matrix(target_metric)), u), matrix(source_metric))


def nullspace(rows, columns):
    """Exact kernel of a possibly rectangular linear system."""
    work = [[scalar(x) for x in row] for row in rows]
    if any(len(row) != columns for row in work):
        raise ValueError('Linear system has inconsistent dimensions')
    pivots = []
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        divisor = work[rank][column]
        work[rank] = [x / divisor for x in work[rank]]
        for row in range(len(work)):
            if row != rank:
                factor = work[row][column]
                work[row] = [x - factor*y for x, y in zip(work[row], work[rank])]
        pivots.append(column)
        rank += 1
    basis = []
    for free in (i for i in range(columns) if i not in pivots):
        vector = [Q(0)] * columns
        vector[free] = Q(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -work[row][free]
        basis.append(tuple(vector))
    return tuple(basis)


def invariant_symmetric_forms(generators):
    """All fixed symmetric H with A^T H + H A = 0 on the two-mode carrier."""
    basis = (matrix(((1, 0), (0, 0))), K, matrix(((0, 0), (0, 1))))
    rows = []
    for a in generators:
        a = matrix(a)
        if len(a) != 2:
            raise ValueError('This invariant-form solver is for the two-mode carrier')
        images = [add(mul(transpose(a), h), mul(h, a)) for h in basis]
        rows.extend([images[k][i][j] for k in range(3)]
                    for i in range(2) for j in range(2))
    return tuple(matrix(((x, y), (y, z))) for x, y, z in nullspace(rows, 3))


@dataclass(frozen=True)
class Move:
    source: str
    target: str
    operators: tuple
    sheet: int = 0

    def __post_init__(self):
        if (not isinstance(self.source, str) or not self.source
                or not isinstance(self.target, str) or not self.target):
            raise ValueError('Move endpoints need labels')
        if type(self.sheet) is not int:
            raise TypeError('This sheet sector is an integer winding ledger')
        entries = tuple((name, matrix(value)) for name, value in self.operators)
        if not entries or len({name for name, _ in entries}) != len(entries):
            raise ValueError('Use one invertible transport per named carrier')
        for name, value in entries:
            if not isinstance(name, str) or not name:
                raise ValueError('Transport carriers need names')
            inverse(value)
        object.__setattr__(self, 'operators', tuple(sorted(entries)))

    def then(self, later):
        if self.target != later.source:
            raise ValueError('Move endpoints do not compose')
        earlier_ops, later_ops = dict(self.operators), dict(later.operators)
        if earlier_ops.keys() != later_ops.keys():
            raise ValueError('Moves must transport the same carrier alphabet')
        return Move(self.source, later.target, tuple(
            (name, mul(later_ops[name], value)) for name, value in self.operators),
            self.sheet + later.sheet)

    def reverse(self):
        return Move(self.target, self.source, tuple(
            (name, inverse(value)) for name, value in self.operators), -self.sheet)

    def reframe(self, source_changes, target_changes, *, source_sheet=0, target_sheet=0):
        if type(source_sheet) is not int or type(target_sheet) is not int:
            raise TypeError('Sheet reference changes are integers')
        return Move(self.source, self.target, tuple(
            (name, mul(mul(matrix(target_changes[name]), value),
                       inverse(matrix(source_changes[name]))))
            for name, value in self.operators), self.sheet + target_sheet - source_sheet)


def integrate_path(initial, steps):
    """T_next = rho(U)T - L - open_residue, in declared step order.

    Each step is (Move, tensor ledger, tensor open residue); all tensors
    have the initial type. Returns the terminal value, composite Move,
    and total ledger-plus-residue transported to the terminal carrier.
    """
    value = initial
    accumulated = tensor_scale(initial, 0)
    composite = None
    for move, ledger, residue in steps:
        composite = move if composite is None else composite.then(move)
        correction = tensor_add(ledger, residue)
        value = tensor_sub(transport(value, dict(move.operators)), correction)
        accumulated = tensor_add(transport(accumulated, dict(move.operators)), correction)
    if composite is None:
        raise ValueError('Provide at least one transition')
    return value, composite, accumulated


def return_report(loop, signatures, *, sheet_ledger=0):
    """Audit entire tensor types, full carriers, and an independent sheet register."""
    if loop.source != loop.target:
        raise ValueError('Return audit needs a closed base path')
    if type(sheet_ledger) is not int:
        raise TypeError('Sheet ledger must be an integer fixed by the protocol')
    ops = dict(loop.operators)
    visible = [induced_matrix(slots, ops) == identity(prod(s.dimension for s in slots))
               for slots in signatures]
    carrier_return = all(value == identity(len(value)) for value in ops.values())
    return {'tensor_type_returns': visible, 'full_carrier_return': carrier_return,
            'sheet_residue': loop.sheet - sheet_ledger,
            'full_unledgered_return': carrier_return and loop.sheet == 0,
            'full_return_after_sheet_ledger': carrier_return and loop.sheet == sheet_ledger}
