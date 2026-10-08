"""HX1: the fall is a turning of the radial cut toward the third cut by an imaginary angle;
crossing the horizon adds a real quarter turn and exchanges the two cuts. Exact sympy."""
import json, os, sympy as sp
I = sp.I; I2 = sp.eye(2); Z = sp.zeros(2)
R = sp.Matrix([[0, -1], [1, 0]]); K = sp.Matrix([[1, 0], [0, -1]]); S = R*K
C1, C2, C3 = K, S, I*R
psi, phi = sp.symbols('psi phi', real=True)
sm = lambda M: M.applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp)))
out = {}
G = C1*C3                                                   # generator of turning C1 toward C3
out['T1_generator'] = G*G == -I2 and sp.simplify(R + I*C3) == Z and sm(C1*C2*C3 - I*I2) == Z
Exp = lambda a: sp.cos(a)*I2 + sp.sin(a)*G                  # native Exp of a generator with square -1
turned = lambda a: C1*Exp(a)                                # = cos(a) C1 + sin(a) C3
out['T1_turning'] = sm(turned(phi) - (sp.cos(phi)*C1 + sp.sin(phi)*C3)) == Z
# T2: the frame cut cosh(psi) K + sinh(psi) R is the radial cut turned by the imaginary angle -i psi
frame = sp.cosh(psi)*K + sp.sinh(psi)*R
out['T2_fall_is_imaginary_turn'] = sm(frame - turned(-I*psi)) == Z
# T3: past the horizon the angle gains a real quarter turn: the two cuts exchange
inside = sm(turned(sp.pi/2 - I*psi))
out['T3_exchange'] = sm(inside - (sp.cosh(psi)*C3 + I*sp.sinh(psi)*C1)) == Z
out['T3_same_as_QT2'] = sm(inside - frame.subs(psi, psi + I*sp.pi/2).applyfunc(sp.expand_trig)) == Z
# T4: real angle: the same generator closes; pi gives the opposite cut, 2 pi the cut itself
out['T4_real_turn_closes'] = sm(turned(sp.pi) + C1) == Z and sm(turned(2*sp.pi) - C1) == Z and sm(turned(sp.pi/2) - C3) == Z
# T5: open and closed readings share the generator; their angles differ by the factor iota = C1 C2 C3
out['T5_one_generator'] = sm(turned(-(C1*C2*C3)[0, 0]*psi) - frame) == Z
# T6: square of the turned cut is 1 for every complex angle: it stays a cut all the way through
a = sp.symbols('a')
out['T6_stays_a_cut'] = sm(turned(a)*turned(a) - I2) == Z
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'HX1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
