import random

def miller(n, k=10):
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