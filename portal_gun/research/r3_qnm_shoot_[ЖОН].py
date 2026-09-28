"""
R3-QNM [ЖОН] — независимый first-principles расчёт QNM скаляра (fast Newton)
================================================================================
Метрика MT-типа проекта:  A=1,  b(r)=b0^2/r  ('ellis'-флейринг)
  ds^2 = -dt^2 + (1-b0^2/r^2)^{-1} dr^2 + r^2 dOmega^2,  r>=b0
Тортоида:  dr*/dr = (1-b0^2/r^2)^{-1/2}  =>  r = sqrt(r*^2 + b0^2)
Потенциал скаляра (A=1, l):  V(r*) = l(l+1)/(r*^2 + b0^2)

Метод: решение ВУ psi''+(w^2-V)psi=0 на [0,R], симметрия -> чётная/нечётная
моды; условие «излучение на +inf»:  F(w) := psi'(R) - i*w*psi(R) = 0,
корни ищем комплексным Ньютоном F(w)=0.  b0 масштаб: w*b0 = f(l) only.

Литературная сверка: drainhole-Эллис (g_thth=r^2+b0^2) даёт l=1 фундаментальную
~ Re(w*b0)≈0.54, Im≈-0.07. Здесь проверяем, ТА ЖЕ ли цифра у ФЛЕЙРИНГ-ВАРИАНТА
проекта (C=r). Если отличается -> расхождение цифр, пометить (правило команды).
"""
import numpy as np
from scipy.integrate import solve_ivp

b0n = 1.0
Rmax = 45.0


def V(x, l):
    return l * (l + 1) / (x**2 + b0n**2)


def F(w, l, parity):
    def rhs(rst, y):
        psi, psip = y
        return [psip, -(w**2 - V(rst, l)) * psi]
    if parity == 'even':
        y0 = [1.0 + 0.0j, 0.0j]
    else:
        y0 = [0.0j, 1.0 + 0.0j]
    sol = solve_ivp(rhs, [0, Rmax], y0, rtol=1e-11, atol=1e-13,
                    max_step=0.15, dense_output=False)
    psi, psip = sol.y[0][-1], sol.y[1][-1]
    return psip - 1j * w * psi


def dF(w, l, parity, eps=1e-6):
    return (F(w + eps, l, parity) - F(w - eps, l, parity)) / (2 * eps)


def newton(w0, l, parity, itmax=40, tol=1e-12):
    w = complex(w0)
    for it in range(itmax):
        fn = F(w, l, parity)
        df = dF(w, l, parity)
        if abs(fn) < tol * max(1.0, abs(w)):
            return w, abs(fn), it
        if df == 0:
            return None
        dw = -fn / df
        w = w + dw
        if abs(dw) < tol * max(1.0, abs(w)):
            return w, abs(fn), it
    return w, abs(F(w, l, parity)), itmax


if __name__ == '__main__':
    import time
    t0 = time.time()
    starts = [1.2-0.1j, 1.35-0.2j, 1.41-0.1j, 1.5-0.3j, 1.6-0.2j]
    for l in [1, 2]:
        found = {}
        for parity in ['even', 'odd']:
            vals = []
            for s in starts:
                r = newton(s, l, parity)
                if r:
                    w, aF, nit = r[0], r[1], r[2]
                    key = (round(w.real, 4), round(w.imag, 4))
                    if key not in found:
                        vals.append((w, aF, nit))
            # отчёт
            print(f"l={l} {parity}:")
            for (w, aF, nit) in vals:
                print(f"   w={w.real:+.6f}{w.imag:+.6f}i   |F|={aF:.1e}  iters={nit}")
    print("time %.1f s" % (time.time() - t0))