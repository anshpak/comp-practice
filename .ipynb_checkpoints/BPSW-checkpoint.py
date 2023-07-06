from math import sqrt
from math import gcd
from sympy import jacobi_symbol

def bpsw(n):
    if n == 2:
        return True
    if n < 2 or n%2 == 0: #or is_power_of_number(n):
        return False
    dd = 5
    while True:
        g = gcd(n, abs(dd))
        if 1<g and g<n:
            return False
        if jacobi_symbol(dd, n) == -1 and g == 1:
            return True
        dd = -dd+2 if dd<0 else -dd-2
    p=1
    q=(p*p-dd)//4
    d=n+1
    s=0
    while (d & 1) == 0:
        s += 1
        d >>= 1
    u=1
    v=p
    u2m=1
    v2m=p
    qm=q
    qm2=q*2
    qkd=q
    for mask in range(2, d+1, 2):
        u2m = (u2m * v2m) % n
        v2m = (v2m * v2m) % n
        while v2m < qm2:
            v2m += n
        v2m -= qm2
        qm = (qm * qm) % n
        qm2 = qm * 2
        if d & mask:
            t1 = (u2m * v) % n
            t2 = (v2m * u) % n
            t3 = (v2m * v) % n
            t4 = (((u2m * u) % n) * dd) % n
            u = t1 + t2
            if u & 1:
                u += n
            u = (u >> 1) % n
            v = t3 + t4
            if v & 1:
                v += n
            v = (v >> 1) % n
            qkd = (qkd * qm) % n
    if u == 0 or v == 0:
        return True
    qkd2 = qkd*2
    for r in range(1, s):
        v = (v * v) % n - qkd2
        if v < 0:
            v += n
        if v < 0:
            v += n
        if v >= n:
            v -= n
        if v >= n:
            v -= n
        if v == 0:
            return True
        if r < s-1:
            qkd = (qkd * qkd) % n
            qkd2 = qkd * 2
    return False

def is_power_of_number(n):
    root = int(sqrt(n))
    for i in range(2, root + 1):
        temp = n
        count = 0
        while temp % i == 0:
            temp //= i
            count += 1
        if count > 1:
            return True
    return False