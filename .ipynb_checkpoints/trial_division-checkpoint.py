def trial_division(n):
    if n == 1 or n == 0:
        return False
    d = 2
    while d<= n**0.5:
        if n % d == 0:
            return False
        d += 1
    return True