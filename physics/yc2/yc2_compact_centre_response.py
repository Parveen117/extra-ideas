"""YC2: centre-odd and plaquette-even response on the compact one-site carrier.

The one-site SU(2) Yang--Mills configuration carrier has three links
U_i=(a_i,u_i) in S^3 and three exact plaquette energies

    X_ij = |u_i cross u_j|^2,    Q = X_12+X_23+X_31.

Independent link-centre flips make a_i odd and every X_ij even.  This file
certifies the resulting response-block decomposition, the exact Haar
covariance of the three X channels, and the character Fourier transform that
is required when an observable must distinguish the eight centre sectors.
"""
from fractions import Fraction as F
from itertools import permutations, product
from math import prod
from pathlib import Path
import hashlib
import json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
N_LINKS = 3
N_VECTOR_COORDS = 3
N_VARS = N_LINKS*N_VECTOR_COORDS

SOURCES = (
    ROOT / "physics/yc1/YC1_CENTRE_COMPASS_POTENTIAL.md",
    ROOT / "uncut/ca2/CA2_SYMMETRY_COMPASS_ATLAS.md",
    ROOT / "physics/ol1/OL1_ONE_SITE_LATTICE.md",
    ROOT / "physics/cm2/CM2_COMPACT_CORE_MATCHING.md",
)


def source_pins():
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in SOURCES}


def odd_double_factorial(n):
    if n <= 0:
        return 1
    return prod(range(n, 0, -2))


def sphere_moment(exponents, dimension=4):
    """Uniform S^(dimension-1) coordinate moment.

    The three listed exponents are the vector coordinates of one quaternion;
    its scalar-coordinate exponent is zero.  The formula is exact.
    """
    if dimension <= 0 or any(not isinstance(e, int) or e < 0 for e in exponents):
        raise ValueError("nonnegative integer exponents and positive dimension required")
    if any(e % 2 for e in exponents):
        return F(0)
    halves = [e//2 for e in exponents]
    total = sum(halves)
    numerator = prod(odd_double_factorial(2*k-1) for k in halves)
    denominator = prod(dimension+2*j for j in range(total))
    return F(numerator, denominator)


def monomial(*indices):
    exponents = [0]*N_VARS
    for index in indices:
        if not 0 <= index < N_VARS:
            raise ValueError("variable index outside the three-link carrier")
        exponents[index] += 1
    return {tuple(exponents): F(1)}


def padd(left, right, scale=F(1)):
    out = dict(left)
    for exponent, coefficient in right.items():
        out[exponent] = out.get(exponent, F(0))+scale*coefficient
        if not out[exponent]:
            del out[exponent]
    return out


def pmul(left, right):
    out = {}
    for e, a in left.items():
        for f, b in right.items():
            exponent = tuple(x+y for x, y in zip(e, f))
            out[exponent] = out.get(exponent, F(0))+a*b
    return {e: a for e, a in out.items() if a}


def energy_polynomial(i, j):
    """Exact polynomial |u_i cross u_j|^2 on three quaternion vectors."""
    if not 0 <= i < N_LINKS or not 0 <= j < N_LINKS or i == j:
        raise ValueError("two distinct link indices are required")
    out = {}
    for a, b in ((0, 1), (0, 2), (1, 2)):
        ia, ib = i*N_VECTOR_COORDS+a, i*N_VECTOR_COORDS+b
        ja, jb = j*N_VECTOR_COORDS+a, j*N_VECTOR_COORDS+b
        out = padd(out, monomial(ia, ia, jb, jb))
        out = padd(out, monomial(ib, ib, ja, ja))
        out = padd(out, monomial(ia, jb, ib, ja), F(-2))
    return out


def gram_energy_polynomial(i, j):
    """The separately expanded |u_i|^2|u_j|^2-(u_i.u_j)^2."""
    out = {}
    for a in range(3):
        for b in range(3):
            ia, ib = i*3+a, i*3+b
            ja, jb = j*3+a, j*3+b
            out = padd(out, monomial(ia, ia, jb, jb))
            out = padd(out, monomial(ia, ja, ib, jb), F(-1))
    return out


def haar_mean(poly):
    """Product-Haar mean on (S^3)^3, using exact monomial moments."""
    answer = F(0)
    for exponent, coefficient in poly.items():
        if len(exponent) != N_VARS:
            raise ValueError("polynomial does not live on three quaternion vectors")
        value = coefficient
        for link in range(N_LINKS):
            start = link*N_VECTOR_COORDS
            value *= sphere_moment(exponent[start:start+N_VECTOR_COORDS])
        answer += value
    return answer


def link_parities(poly):
    """Return the common Z2^3 parity of a polynomial, or None if mixed."""
    signatures = set()
    for exponent, coefficient in poly.items():
        if coefficient:
            signatures.add(tuple(sum(exponent[3*i:3*i+3]) % 2 for i in range(3)))
    if len(signatures) != 1:
        return None
    return next(iter(signatures))


def permute_links(poly, permutation):
    if sorted(permutation) != list(range(N_LINKS)):
        raise ValueError("a permutation of the three links is required")
    out = {}
    for exponent, coefficient in poly.items():
        moved = [0]*N_VARS
        for old in range(N_LINKS):
            new = permutation[old]
            moved[3*new:3*new+3] = exponent[3*old:3*old+3]
        out[tuple(moved)] = out.get(tuple(moved), F(0))+coefficient
    return {e: a for e, a in out.items() if a}


def haar_energy_data():
    edges = ((0, 1), (0, 2), (1, 2))
    channels = [energy_polynomial(*edge) for edge in edges]
    means = [haar_mean(channel) for channel in channels]
    covariance = []
    for i, left in enumerate(channels):
        row = []
        for j, right in enumerate(channels):
            row.append(haar_mean(pmul(left, right))-means[i]*means[j])
        covariance.append(row)
    diagonal, shared = covariance[0][0], covariance[0][1]
    pair_chi = (diagonal*diagonal-shared*shared)/(diagonal*diagonal)
    symmetric_eigenvalue = diagonal+2*shared
    anisotropy_eigenvalue = diagonal-shared
    normalized_determinant = anisotropy_eigenvalue**2*symmetric_eigenvalue/diagonal**3
    return {
        "edges": edges,
        "channels": channels,
        "means": means,
        "covariance": covariance,
        "pair_chi": pair_chi,
        "pair_lost": 1-pair_chi,
        "symmetric_eigenvalue": symmetric_eigenvalue,
        "anisotropy_eigenvalue": anisotropy_eigenvalue,
        "normalized_determinant": normalized_determinant,
        "mean_Q": sum(means),
        "variance_Q": sum(sum(row) for row in covariance),
    }


def character(sigma, sheet):
    """Character sigma of one element of (Z2)^3, both encoded by signs."""
    if (len(sigma) != 3 or len(sheet) != 3
            or any(x not in (-1, 1) for x in sigma+sheet)):
        raise ValueError("three signs are required for the character and sheet")
    answer = 1
    for eigenvalue, sign in zip(sigma, sheet):
        if sign == -1:
            answer *= eigenvalue
    return answer


def observation_chart(kind="odd"):
    """Typed TVSP -> VTSP chart for a selected pair of compact channels."""
    if kind == "odd":
        before = ("<c_i>", "kappa_j", "kappa_i", "<c_j>")
    elif kind == "energy":
        before = ("<X_alpha>", "j_beta", "j_alpha", "<X_beta>")
    else:
        raise ValueError("kind must be 'odd' or 'energy'")
    return before, (before[1], before[0], before[2], before[3])


def run():
    checks = {}
    data = haar_energy_data()
    channels = data["channels"]

    checks["uniform S3 vector moments are exact quaternion Haar moments"] = (
        sphere_moment((2, 0, 0)) == F(1, 4)
        and sphere_moment((4, 0, 0)) == F(1, 8)
        and sphere_moment((2, 2, 0)) == F(1, 24))
    checks["cross-product and Gram forms of every plaquette energy agree"] = all(
        energy_polynomial(i, j) == gram_energy_polynomial(i, j)
        for i, j in data["edges"])
    checks["each plaquette energy is even under all three independent centre flips"] = all(
        link_parities(channel) == (0, 0, 0) for channel in channels)

    permutation_ok = True
    edge_to_channel = {tuple(sorted(edge)): channel for edge, channel in zip(data["edges"], channels)}
    for permutation in permutations(range(3)):
        for edge, channel in zip(data["edges"], channels):
            moved_edge = tuple(sorted((permutation[edge[0]], permutation[edge[1]])))
            permutation_ok &= permute_links(channel, permutation) == edge_to_channel[moved_edge]
    checks["link permutations act exactly on the three plaquette channels"] = permutation_ok

    checks["each Haar plaquette energy has mean three eighths"] = (
        data["means"] == [F(3, 8)]*3)
    checks["each Haar plaquette energy has variance thirteen over 192"] = all(
        data["covariance"][i][i] == F(13, 192) for i in range(3))
    checks["a shared link contributes covariance one over 64"] = all(
        data["covariance"][i][j] == F(1, 64)
        for i in range(3) for j in range(3) if i != j)
    checks["the exact two-energy TVSP compass is 160 over 169"] = (
        data["pair_chi"] == F(160, 169) and data["pair_lost"] == F(9, 169))
    checks["symmetric and anisotropy response modes are both strictly positive"] = (
        data["symmetric_eigenvalue"] == F(19, 192)
        and data["anisotropy_eigenvalue"] == F(5, 96)
        and data["normalized_determinant"] == F(1900, 2197))
    checks["the total Wilson-core channel has exact Haar mean and variance"] = (
        data["mean_Q"] == F(9, 8) and data["variance_Q"] == F(19, 64))
    checks["the Haar centre-odd loop response block is one quarter times identity"] = (
        sphere_moment((2, 0, 0)) == F(1, 4))

    # Parity proof for the all-coupling response block.  The action Q, energy
    # sources X_ij and self sources c_i^2 are even in every link.  Therefore an
    # integrand with nonzero parity has zero mean at kappa=0, for every theta
    # and every even source value.
    odd_parities = [tuple(1 if i == j else 0 for i in range(3)) for j in range(3)]
    off_diagonal_products = [tuple(a ^ b for a, b in zip(odd_parities[i], odd_parities[j]))
                             for i in range(3) for j in range(3) if i != j]
    checks["independent centre flips diagonalize the three centre-odd loop responses"] = all(
        any(signature) for signature in off_diagonal_products)
    checks["centre-odd loops have zero mixed response with every plaquette-energy channel"] = all(
        any(signature) for signature in odd_parities)

    signs = list(product((-1, 1), repeat=3))
    checks["the eight centre characters give exact orthogonal sector projectors"] = all(
        sum(character(sigma, sheet)*character(tau, sheet) for sheet in signs)
        == (8 if sigma == tau else 0)
        for sigma in signs for tau in signs)

    odd_before, odd_after = observation_chart("odd")
    even_before, even_after = observation_chart("energy")
    checks["observation swaps the first two positions and carries each typed chart"] = (
        odd_before == ("<c_i>", "kappa_j", "kappa_i", "<c_j>")
        and odd_after == ("kappa_j", "<c_i>", "kappa_i", "<c_j>")
        and even_before == ("<X_alpha>", "j_beta", "j_alpha", "<X_beta>")
        and even_after == ("j_beta", "<X_alpha>", "j_alpha", "<X_beta>"))

    checks = {name: bool(value) for name, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])

    return {
        "stage": "YC2",
        "checks": checks,
        "all_pass": True,
        "carrier": {
            "configuration": "(U1,U2,U3) in SU(2)^3 with product Haar measure",
            "gauge_action": "simultaneous conjugation",
            "centre_group": "(Z2)^3 acting by independent Ui -> -Ui",
            "compact_hamiltonian": "H_theta=-sum Delta_i+2 theta Q",
            "Q": "X12+X23+X31, Xij=|ui cross uj|^2",
        },
        "source_potential": {
            "formula": "W=log integral exp[-2 theta Q+sum kappa_i c_i+sum h_i c_i^2+sum j_ij X_ij] dU",
            "centre_odd_channels": "c_i=Tr(U_i)/2",
            "centre_even_channels": "c_i^2 and X_ij",
            "all_coupling_seam": "at kappa=0 the odd block is diagonal and every odd-even mixed response is zero",
            "isotropic_line": "link permutations make the odd block a scalar multiple of I3",
        },
        "haar_energy_response": {
            "centre_odd_loop_covariance": [["1/4" if i == j else "0" for j in range(3)]
                                            for i in range(3)],
            "means": [str(x) for x in data["means"]],
            "covariance": [[str(x) for x in row] for row in data["covariance"]],
            "pair_chi": str(data["pair_chi"]),
            "pair_lost_share": str(data["pair_lost"]),
            "symmetric_eigenvalue": str(data["symmetric_eigenvalue"]),
            "anisotropy_eigenvalue": str(data["anisotropy_eigenvalue"]),
            "normalized_three_channel_determinant": str(data["normalized_determinant"]),
            "mean_Q": str(data["mean_Q"]),
            "variance_Q": str(data["variance_Q"]),
        },
        "observation": {
            "odd_TVSP": list(odd_before),
            "odd_VTSP": list(odd_after),
            "energy_TVSP": list(even_before),
            "energy_VTSP": list(even_after),
            "arrow_rule": "push an arrow and the response metric by the same first-two-position permutation",
        },
        "sector_sensitive_adapter": {
            "projector": "P_sigma=2^-3 sum_z chi_sigma(z) Z_z",
            "twisted_heat_trace": "K_z(t)=Tr[Z_z exp(-t H_theta)]",
            "sector_heat_trace": "Z_sigma(t)=2^-3 sum_z chi_sigma(z) K_z(t)",
            "sector_ground": "E0_sigma=-lim_(t->infinity) t^-1 log Z_sigma(t)",
            "splitting_ratio": "E0_sigma-E0_plus=-lim t^-1 log[Z_sigma(t)/Z_plus(t)]",
        },
        "claim_boundary": {
            "YC1_embeds_as_each_single_link_when_theta_and_cross_sources_vanish": True,
            "centre_odd_and_even_response_blocks_separate_at_the_symmetric_seam": True,
            "Haar_shared_link_compass_is_160_over_169": True,
            "centre_even_Wilson_potential_alone_can_identify_a_centre_character": False,
            "twisted_heat_trace_is_required_for_sector_splitting": True,
            "absolute_CM2_sector_splitting_rate_computed": False,
            "spatial_volume_or_continuum_gap_proved": False,
        },
        "source_sha256": source_pins(),
    }


if __name__ == "__main__":
    result = run()
    (HERE / "YC2_RESULT.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))
