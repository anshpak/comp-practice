def trial_division(n):
    d = 2
    while d<= n**0.5:
        if n % d == 0:
            return False
        d += 1
    return True