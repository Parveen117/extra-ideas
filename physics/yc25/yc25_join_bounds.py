"""YC25 exact budget calculator; no spectral truncation or certificate suite."""

from fractions import Fraction as F
from math import factorial


RADIUS = F(1, 32)
EXP_CAP = F(9, 7)
CASES = (
    ("cubes", F(2), 24, F(6), F(27, 5), F(1, 840)),
    ("rectangles, theta<=2", F(15, 8), 48, F(64, 5), F(7, 3), F(1, 1792)),
    ("rectangles, theta<=1", F(15, 8), 48, F(228, 25), F(12, 5), F(1, 1680)),
)


def budget(mu, incidence, seed_coefficient, delta, coupling):
    beta = incidence * abs(coupling)
    q = 4 * beta / mu * EXP_CAP * (2 + 8 * (1 + 2 * RADIUS))
    b = 8 * beta / mu * EXP_CAP
    seed = seed_coefficient * abs(coupling)
    return {
        "seed": seed,
        "q": q,
        "b_G": b,
        "ball_reserve": RADIUS - seed - q * RADIUS,
        "physical_gap_floor": delta * (1 - b),
    }


def main():
    x = F(1, 4)
    exp_bound = sum((x**j / factorial(j) for j in range(4)), F(0))
    exp_bound += x**4 / factorial(4) / (1 - x / 5)
    print(f"exp(1/4) <= {exp_bound} < {EXP_CAP}; reserve={EXP_CAP-exp_bound}")
    for name, mu, incidence, seed, delta, endpoint in CASES:
        print(f"\n{name}: extensive floor={mu}")
        for coupling in (F(0), endpoint):
            result = budget(mu, incidence, seed, delta, coupling)
            print(f"  lambda={coupling}: " + ", ".join(f"{k}={v}" for k, v in result.items()))
    name, mu, incidence, seed, delta, _ = CASES[-1]
    old_window = budget(mu, incidence, seed, delta, F(1, 4000))
    print(f"\n{name}, old cap 1/4000: physical_gap_floor={old_window['physical_gap_floor']}")
    print("\nZero budget: seed=q=b_G=0; reference physical floor is retained.")
    print("State/operator/metric specialization follows from YC25 equations (6),(12)-(15).")


if __name__ == "__main__":
    main()
