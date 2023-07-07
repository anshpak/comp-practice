from tests import PrimeTest

class Numbers:

    @staticmethod
    def factorize(n):
        if n == 0:
            return (0, 0)
        s = 0
        while n % 2 == 0:
            n //= 2
            s += 1
        return (n, s)

    @staticmethod
    def prot(p):
        l = p - 1
        k, n = Numbers.factorize(l)
        if k % 2 != 0 and k < 2**n:
            return True
        else:
            return False

    @staticmethod
    def is_mills_number(n):
        A = 1.3063778838630806904686144926
        for i in range(1,8):
            if round(A**(3**i)) == n:
                return True
        return False

    @staticmethod
    def is_mersenne_number(n):
        n += 1
        p = 0
        while n != 1:
            if n & 1:
                return False
            n >>= 1
            p += 1
        if not PrimeTest.trial_divisions(p):
                return False
        return True

    @staticmethod
    def is_kallen_number(p):
        e_list = [1, 141, 4713, 5795, 6611, 18496, 32292, 32469, 59656, 90825, 262419, 361275, 481899, 1354828, 6328548, 6679881 ]
        for n in e_list:
            if n * 2**n + 1 == p:
                return True
        return False