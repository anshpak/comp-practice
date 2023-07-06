import random
import math
import sympy
from math import sqrt
import copy
from jacobi import *

class PrimeTest:
    @staticmethod
    def trial_divisions(n):
        if n == 1:
            return False
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
        if n == 2:
            return True
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
        if n == 1:
            return False
        if n == 2:
            return True
        from random import randint
        if n % 2 == 0:
            return False
        s, t = PrimeTest.__st(n)
        for _ in range(k):
            a = randint(2, n-1)
            if not PrimeTest.__witness_q(a, n, s, t):
                return False
        return True
    
    @staticmethod
    def __prime_factors(n):
        i = 2
        result = []
        while i * i <= n:
            count = 0
            while n % i == 0:
                count += 1
                result.append((i, count))
                n = n // i
            i = i + 1
        if n > 1:
            result.append((n,1))
        return result

    @staticmethod
    def __v(q, t):
        ans = 0
        while(t % q == 0):
            ans +=1
            t//= q
        return ans

    @staticmethod
    def __e(t):
        s = 1
        q_list = []
        for q in range(2, t+2):
            if t%(q-1) == 0 and PrimeTest.trial_divisions(q):
                s *= q ** (1+PrimeTest.__v(q,t))
                q_list.append(q)
        return 2*s, q_list
    
    @staticmethod
    def __smallest_primitive_root(q):
        for r in range(2, q):
            s = set({})
            m = 1
            for i in range(1, q):
                m = (m*r) % q
                s.add(m)
            if len(s) == q-1:
                return r
        return None
    
    @staticmethod
    def __calc_f(q):
        g = PrimeTest.__smallest_primitive_root(q)
        m = {}
        for x in range(1,q-1):
            m[pow(g,x,q)] = x
        f = {}
        for x in range(1,q-1):
            f[x] = m[ (1-pow(g,x,q))%q ]
        return f
    
    @staticmethod
    def __calc_J_ab(p, k, q, a, b):
        j_ret = JacobiSum(p,k,q)
        f = PrimeTest.__calc_f(q)
        for x in range(1,q-1):
            pk = p**k
            if (a*x+b*f[x]) % pk < j_ret.m:
                j_ret.coef[(a*x+b*f[x]) % pk] += 1
            else:
                r = (a*x+b*f[x]) % pk - p**(k-1)
                while r>=0:
                    j_ret.coef[r] -= 1
                    r-= p**(k-1)
        return j_ret


    @staticmethod
    def __calc_J(p, k, q):
        return PrimeTest.__calc_J_ab(p, k, q, 1, 1)

    @staticmethod
    def __calc_J3(p, k, q):
        j2q = PrimeTest.__calc_J(p, k, q)
        j21 = PrimeTest.__calc_J_ab(p, k, q, 2, 1)
        j_ret = j2q * j21
        return j_ret

    @staticmethod
    def __calc_J2(p, k, q):
        j31 = PrimeTest.__calc_J_ab(2, 3, q, 3, 1)
        j_conv = JacobiSum(p, k, q)
        for i in range(j31.m):
            j_conv.coef[i*(p**k)//8] = j31.coef[i]
        j_ret = j_conv * j_conv
        return j_ret
    
    @staticmethod
    def __APRtest_step4a(p, k, q, N):    
        J = PrimeTest.__calc_J(p, k, q)
        s1 = JacobiSum(p,k,q).one()
        for x in range(p**k):
            if x % p == 0:
                continue
            t = J.sigma_inv(x)
            t = t.modpow(x, N)
            s1 = s1 * t
            s1.mod(N)

        r = N % (p**k)
        s2 = s1.modpow(N//(p**k), N)
        J_alpha = JacobiSum(p,k,q).one()
        for x in range(p**k):
            if x % p == 0:
                continue
            t = J.sigma_inv(x)
            t = t.modpow((r*x)//(p**k), N)
            J_alpha = J_alpha * t
            J_alpha.mod(N)

        S = (s2 * J_alpha).mod(N)
        exist, h = S.is_root_of_unity(N)

        if not exist:
            return False, None
        else:
            if h%p!=0:
                l_p = 1
            else:
                l_p = 0
            return True, l_p
    
    @staticmethod
    def __APRtest_step4b(p, k, q, N):
        J = PrimeTest.__calc_J3(p, k, q)
        s1 = JacobiSum(p,k,q).one()
        for x in range(p**k):
            if x % 8 not in [1,3]:
                continue
            t = J.sigma_inv(x)
            t = t.modpow(x, N)
            s1 = s1 * t
            s1.mod(N)

        r = N % (p**k)
        s2 = s1.modpow(N//(p**k), N)

        J_alpha = JacobiSum(p,k,q).one()
        for x in range(p**k):
            if x % 8 not in [1,3]:
                continue
            t = J.sigma_inv(x)
            t = t.modpow((r*x)//(p**k), N)
            J_alpha = J_alpha * t
            J_alpha.mod(N)

        if N%8 in [1,3]:
            S = (s2 * J_alpha ).mod(N)
        else:
            J2_delta = PrimeTest.__calc_J2(p,k,q)
            S = (s2 * J_alpha * J2_delta).mod(N)

        exist, h = S.is_root_of_unity(N)

        if not exist:
            return False, None
        else:
            if h%p!=0 and (pow(q,(N-1)//2,N) + 1)%N==0:
                l_p = 1
            else:
                l_p = 0
            return True, l_p
    
    @staticmethod
    def __APRtest_step4c(p, k, q, N):
        J2q = PrimeTest.__calc_J(p, k, q)
        s1 = (J2q * J2q * q).mod(N)
        s2 = s1.modpow(N//4, N)

        if N%4 == 1:
            S = s2
        elif N%4 == 3:
            S = (s2 * J2q * J2q).mod(N)

        exist, h = S.is_root_of_unity(N)

        if not exist:
            return False, None
        else:
            if h%p!=0 and (pow(q,(N-1)//2,N) + 1)%N==0:
                l_p = 1
            else:
                l_p = 0
            return True, l_p

    @staticmethod
    def __APRtest_step4d(p, k, q, N):
        S2q = pow(-q, (N-1)//2, N)
        if (S2q-1)%N != 0 and (S2q+1)%N != 0:
            return False, None
        else:
            if (S2q + 1)%N == 0 and (N-1)%4==0:
                l_p=1
            else:
                l_p=0
            return True, l_p
    
    @staticmethod
    def __APRtest_step4(p, k, q, N):
        if p>=3:
            result, l_p = PrimeTest.__APRtest_step4a(p, k, q, N)
        elif p==2 and k>=3:
            result, l_p = PrimeTest.__APRtest_step4b(p, k, q, N)
        elif p==2 and k==2:
            result, l_p = PrimeTest.__APRtest_step4c(p, k, q, N)
        elif p==2 and k==1:
            result, l_p = PrimeTest.__APRtest_step4d(p, k, q, N)
        return result, l_p
    
    def APRtest(N):
        t_list = [2, 12, 60, 180, 840, 1260, 1680, 2520, 5040, 15120, 55440, 110880, 720720, 1441440, 4324320, 24504480, 73513440]
        if N==1:
            return False
        if N==2 or N==3:
            return True
        for t in t_list:
            et, q_list = PrimeTest.__e(t)
            if N < et*et:
                break
        else:
            return False

        # Step 1
        g = PrimeTest.__gcd(t*et, N)
        if g != 1:
            return False

        # Step 2
        l = {}
        fac_t = PrimeTest.__prime_factors(t)
        for p, k in fac_t:
            if p>=3 and pow(N, p-1, p*p)!=1:
                l[p] = 1
            else:
                l[p] = 0

        # Step 3 & Step 4
        for q in q_list:
            if q == 2:
                continue
            fac = PrimeTest.__prime_factors(q-1)
            for p,k in fac:

                # Step 4
                result, l_p = PrimeTest.__APRtest_step4(p, k, q, N)
                if not result:
                    return False
                elif l_p==1:
                    l[p] = 1

        # Step 5
        for p, value in l.items():
            if value==0:
                count = 0
                i = 1
                found = False
                while count < 30:
                    q = p*i+1
                    if N%q != 0 and PrimeTest.trial_divisions(q) and (q not in q_list):
                        count += 1
                        k = PrimeTest.__v(p, q-1)
                        # Step 4
                        result, l_p = PrimeTest.__APRtest_step4(p, k, q, N)
                        if not result:
                            return False
                        elif l_p == 1:
                            found = True
                            break
                    i += 1
                if not found:
                    return False

        # Step 6
        r = 1
        for t in range(t-1):
            r = (r*N) % et
            if r!=1 and r!= N and N % r == 0:
                return False
        return True