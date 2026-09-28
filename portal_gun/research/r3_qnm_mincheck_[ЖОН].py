"""
R3-QNM [ЖОН] — минимальный численный пересчёт резонанса (один прогон, дешёвый)
===============================================================================
Цель: независимая проверка лит. цифр Ellis-QNM на геометрии проекта.
  V(r*) = L/(r*^2+b0^2), L=l(l+1), b0=1.
ВУ: psi''+(w^2-V)psi=0, [0,R], чётная мода (psi'(0)=0). Условие вылета на +inf:
  F(w)=psi'(R)-i w psi(R)=0, корни комплексным Ньютоном.
Всего: несколько десятков КОРОТКИХ RK4-интегралов. Таймаут 60 с.
"""
import cmath, time
b0n, Rmax, N = 1.0, 35.0, 700

def rk4(w, L, parity):
    h = Rmax / N
    psi, dpsi = (1.0+0.0j, 0.0j) if parity=='even' else (0.0j, 1.0+0.0j)
    w2 = w*w
    for i in range(N):
        x = i*h
        k1 = dpsi
        l1 = -(w2 - L/(x*x+b0n*b0n))*psi
        k2 = dpsi + 0.5*h*l1
        l2 = -(w2 - L/((x+0.5*h)**2+b0n*b0n))*(psi+0.5*h*k1)
        k3 = dpsi + 0.5*h*l2
        l3 = -(w2 - L/((x+0.5*h)**2+b0n*b0n))*(psi+0.5*h*k2)
        k4 = dpsi + h*l3
        l4 = -(w2 - L/((x+h)**2+b0n*b0n))*(psi+h*k3)
        psi += h/6.0*(k1+2*k2+2*k3+k4)
        dpsi += h/6.0*(l1+2*l2+2*l3+l4)
    return dpsi - 1j*w*psi

def F(w, L, parity):
    return rk4(w, L, parity)

def newton(w0, L, parity, itmax=30):
    w = complex(w0)
    for it in range(itmax):
        fn = F(w, L, parity)
        eps = 1e-7*max(1.0, abs(w))
        df = (F(w+eps, L, parity)-F(w-eps, L, parity))/(2*eps)
        if abs(fn) < 1e-9*max(1.0, abs(w)):
            return w, abs(fn), it
        if df == 0:
            return None
        dw = -fn/df
        w += dw
        if abs(dw) < 1e-11*max(1.0, abs(w)):
            return w, abs(F(w, L, parity)), it
    return w, abs(F(w, L, parity)), itmax

if __name__ == '__main__':
    t0 = time.time()
    for l in [1, 2]:
        L = l*(l+1)
        for parity in ['even', 'odd']:
            # стартовые точки вокруг ожидаемой лит. зоны и выше
            for s in [0.5-0.1j, 0.55-0.08j, 0.6-0.1j, 1.3-0.15j, 1.5-0.1j]:
                r = newton(s, L, parity)
                if r:
                    w, aF, n = r
                    print("l=%d %s  start=%.2f%+.2fi -> w=%.6f%+.6fi  |F|=%.1e  (it=%d)  [%.1fs]"
                          % (l, parity, s.real, s.imag, w.real, w.imag, aF, n, time.time()-t0))
    print("TOTAL %.1fs" % (time.time()-t0))