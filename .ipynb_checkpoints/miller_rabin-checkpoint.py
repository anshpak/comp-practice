def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def st(n):
    s = 0
    t = n - 1
    while t % 2 == 0:
        s += 1
        t = t // 2
    return s, t

def witness_q(a, n, s=0, t=0):
    if n % 2 == 0:
        return False
    if s == 0:
        s, t = st(n)    
    if gcd(a,n) > 1:
        return False
    b = pow(a,t,n)
    if b == 1 or b == n-1:
        return True
    for _ in range(1,s):
        b = b**2 % n
        if b == n-1:
            return True
    return False

def miller_rabin(n, k = 10):
    from random import randint
    if n % 2 == 0:
        return False
    s,t = st(n)
    for _ in range(k):
        a = randint(2,n-1)
        if not witness_q(a,n,s,t):
            return False
    return True