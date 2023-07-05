import copy

class JacobiSum:
    def __init__(self, p, k, q):
        self.p = p
        self.k = k
        self.q = q
        self.m = (p-1)*p**(k-1)
        self.pk = p**k
        self.coef = [0]*self.m
        

    def one(self):
        self.coef[0] = 1
        for i in range(1,self.m):
            self.coef[i] = 0
        return self
    

    def mul(self, jac):
        m = self.m
        pk = self.pk
        j_ret=JacobiSum(self.p, self.k, self.q)
        for i in range(m):
            for j in range(m):
                if (i+j)% pk < m:
                    j_ret.coef[(i+j)% pk] += self.coef[i] * jac.coef[j]
                else:
                    r = (i+j) % pk - self.p ** (self.k-1)                    
                    while r>=0:
                        j_ret.coef[r] -= self.coef[i] * jac.coef[j]
                        r-= self.p ** (self.k-1)
        return j_ret


    def __mul__(self, right):
        if type(right) is int:
            j_ret=JacobiSum(self.p, self.k, self.q)
            for i in range(self.m):
                j_ret.coef[i] = self.coef[i] * right
            return j_ret
        else:
            return self.mul(right)
        
    
    def modpow(self, x, n):
        j_ret=JacobiSum(self.p, self.k, self.q)
        j_ret.coef[0]=1
        j_a = copy.deepcopy(self)
        while x>0:
            if x%2==1:
                j_ret = (j_ret * j_a).mod(n)
            j_a = j_a*j_a
            j_a.mod(n)
            x //= 2
        return j_ret
    
    
    def mod(self, n):
        for i in range(self.m):
            self.coef[i] %= n
        return self
    

    def sigma(self, x):
        m = self.m
        pk = self.pk
        j_ret=JacobiSum(self.p, self.k, self.q)
        for i in range(m):
            if (i*x) % pk < m:
                j_ret.coef[(i*x) % pk] += self.coef[i]
            else:
                r = (i*x) % pk - self.p ** (self.k-1)                    
                while r>=0:
                    j_ret.coef[r] -= self.coef[i]
                    r-= self.p ** (self.k-1)
        return j_ret
    
                
    def sigma_inv(self, x):
        m = self.m
        pk = self.pk
        j_ret=JacobiSum(self.p, self.k, self.q)
        for i in range(pk):
            if i<m:
                if (i*x)%pk < m:
                    j_ret.coef[i] += self.coef[(i*x)%pk]
            else:
                r = i - self.p ** (self.k-1)
                while r>=0:
                    if (i*x)%pk < m:
                        j_ret.coef[r] -= self.coef[(i*x)%pk]
                    r-= self.p ** (self.k-1)
        return j_ret
    

    def is_root_of_unity(self, N):
        m = self.m
        p = self.p
        k = self.k
        one = 0
        for i in range(m):
            if self.coef[i]==1:
                one += 1
                h = i
            elif self.coef[i] == 0:
                continue
            elif (self.coef[i] - (-1)) %N != 0:
                return False, None
        if one == 1:
            return True, h
        for i in range(m):
            if self.coef[i]!=0:
                break
        r = i % (p**(k-1))
        for i in range(m):
            if i % (p**(k-1)) == r:
                if (self.coef[i] - (-1))%N != 0:
                    return False, None
            else:
                if self.coef[i] !=0:
                    return False, None
        return True, (p-1)*p**(k-1)+ r