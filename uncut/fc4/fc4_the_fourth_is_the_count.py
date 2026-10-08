"""FC4: the fourth component of a reading is its count, not a fourth direction. Exact sympy.
Native dagger on a + bK + cS + dR (cut-complex coefficients): conjugate the coefficients, keep K and S, reverse R."""
import json, os, sympy as sp
I = sp.I; I2 = sp.eye(2); Z = sp.zeros(2)
R = sp.Matrix([[0, -1], [1, 0]]); K = sp.Matrix([[1, 0], [0, -1]]); S = R*K
C = [K, S, I*R]
out = {}
zm = lambda M: M.applyfunc(lambda e: sp.simplify(sp.expand(sp.simplify(e).rewrite(sp.exp)))) == Z
def elem(a, b, c, d): return a*I2 + b*K + c*S + d*R
def dagger(a, b, c, d): return (sp.conjugate(a), sp.conjugate(b), sp.conjugate(c), -sp.conjugate(d))
# T1: self-dagger elements form a space of exactly four real components: n, r1, r2, r3 with rho = n + r.C
a, b, c, d = [sp.symbols(f'{x}r', real=True) + I*sp.symbols(f'{x}i', real=True) for x in 'abcd']
cond = [sp.Eq(x, y) for x, y in zip((a, b, c, d), dagger(a, b, c, d))]
free = sp.solve(cond, [sp.Symbol('ai', real=True), sp.Symbol('bi', real=True), sp.Symbol('ci', real=True), sp.Symbol('dr', real=True)], dict=True)[0]
out['T1_four_real_components'] = len(free) == 4 and all(v == 0 for v in free.values())
n, r1, r2, r3 = sp.symbols('n r1 r2 r3', real=True)
rho = n*I2 + r1*C[0] + r2*C[1] + r3*C[2]
# T2: its determinant is n^2 - r.r (IN1-T1; EMK-1 with the R-coefficient iota r3)
out['T2_invariant'] = sp.simplify(rho.det() - (n**2 - r1**2 - r2**2 - r3**2)) == 0
# T3: the count part commutes with every cut; the cuts anticommute among themselves: n is not a fourth direction
out['T3_count_is_central'] = all(zm(I2*Cc - Cc*I2) for Cc in C) and all(zm(C[i]*C[j] + C[j]*C[i]) for i in range(3) for j in range(i))
# and there is no fourth cut (FR1, WD2-T2)
w, x, y, z_ = sp.symbols('w x y z')
X = sp.Matrix([[w, x], [y, z_]])
out['T3_no_fourth_cut'] = sp.solve([e for Cc in C for e in (X*Cc + Cc*X)], [w, x, y, z_], dict=True) == [{w: 0, x: 0, y: 0, z_: 0}]
# T4: the count adds when readings are put together; directions add as directions
m, s1, s2, s3 = sp.symbols('m s1 s2 s3', real=True)
rho2 = m*I2 + s1*C[0] + s2*C[1] + s3*C[2]
tot = rho + rho2
out['T4_count_adds'] = sp.simplify(tot.trace()/2 - (n + m)) == 0
# T5: a reading at rest has only its count; seen by a reader in motion the count grows by cosh and a direction part appears
psi = sp.symbols('psi', real=True)
g = sp.cosh(psi/2)*I2 + sp.sinh(psi/2)*C[0]
moved = (g*(n*I2)*g).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp)))
out['T5_motion'] = zm(moved - (n*sp.cosh(psi)*I2 + n*sp.sinh(psi)*C[0])) and sp.simplify(sp.expand(((g*rho*g).det() - rho.det()).rewrite(sp.exp))) == 0
# T6: the turn that closes a reading is the product of the three cuts: nothing outside the three is used
out['T6_turn_is_product_of_cuts'] = zm(C[0]*C[1]*C[2] - I*I2)
# T7: one reader's own time: at rest det = n^2: the invariant is the reader's own count squared
out['T7_own_time_is_own_count'] = sp.simplify(rho.det().subs({r1: 0, r2: 0, r3: 0}) - n**2) == 0
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'FC4_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
