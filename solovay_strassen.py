import random
import math
import sympy

def fast_power(a, d, m):
    b = 1
    while d>0:
        if d & 1:
            b = b*a % m
        a = a*a % m
        d = d >> 1
    return b

def solovay_strassen(n, k):
    for i in range(1,k+1):
        a = random.randint(2,n-1)
        if math.gcd(a,n) > 1:
            return False
        j = fast_power(a,(n-1)//2,n)
        jacobi_symb = sympy.jacobi_symbol(a, n) % n
        if jacobi_symb != j:
            print(jacobi_symb)
            return False
    return True