import random
import math
import sympy
from math import sqrt

class PrimeTest:
    @staticmethod
    def trial_division(n):
        d = 2
        while d <= n ** 0.5:
            if n % d == 0:
                return False
            d += 1
        return True

    @staticmethod
    def miller(n, k = 10):
        if n == 2 or n == 3:
            return True
        if n <= 1 or n % 2 == 0:
            return False

        s = 0
        q = n - 1
        while q % 2 == 0:
            s += 1
            q //= 2

        for _ in range(k):
            a = random.randrange(2, n - 1)
            x = pow(a, q, n)
            if x == 1 or x == n - 1:
                continue

            for _ in range(s - 1):
                x = pow(x, 2, n)
                if x == n - 1:
                    break
            else:
                return False
        return True
    
    @staticmethod
    def __fast_power(a, d, m):
        b = 1
        while d > 0:
            if d & 1:
                b = b * a % m
            a = a * a % m
            d = d >> 1
        return b
    
    @staticmethod
    def __gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    @staticmethod
    def solovay_strassen(n, k):
        for i in range(1, k+1):
            a = random.randint(2,n-1)
            if PrimeTest.__gcd(a,n) > 1:
                return False
            j = PrimeTest.__fast_power(a, (n - 1) // 2, n)
            jacobi_symb = sympy.jacobi_symbol(a, n) % n
            if jacobi_symb != j:
                print(jacobi_symb)
                return False
        return True

    @staticmethod
    def __st(n):
        s = 0
        t = n - 1
        while t % 2 == 0:
            s += 1
            t = t // 2
        return s, t

    @staticmethod
    def __witness_q(a, n, s = 0, t = 0):
        if n % 2 == 0:
            return False
        if s == 0:
            s, t = PrimeTest.__st(n)    
        if PrimeTest.__gcd(a,n) > 1:
            return False
        b = pow(a,t,n)
        if b == 1 or b == n - 1:
                       return True
        for _ in range(1, s):
            b = b**2 % n
            if b == n - 1:
                return True
        return False

    @staticmethod
    def miller_rabin(n, k = 10):
        from random import randint
        if n % 2 == 0:
            return False
        s, t = PrimeTest.__st(n)
        for _ in range(k):
            a = randint(2, n-1)
            if not PrimeTest.__witness_q(a, n, s, t):
                return False
        return True