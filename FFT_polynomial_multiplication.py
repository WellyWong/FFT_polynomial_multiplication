from cmath import exp
from math import pi

#A simple class to simulate n-th root of unity. It is implemented merely for FFT and FPM algorithm
class NthRootOfUnity:
    def __init__(self, n, k=1):
        self.k = k
        self.n = n

    def __pow__(self, other):
        if type(other) is int:
            n = NthRootOfUnity(self.n, self.k * other)
            return n

    def __eq__(self, other):
        if other == 1:
            return abs(self.n) == abs(self.k)

    def __mul__(self, other):
        return exp(2 * 1j * pi * self.k/self.n) * other

    def __repr__(self):
        return str(self.n) + '-th root of unity to the ' + str(self.k)

    @property
    def th(self):
        return abs(self.n//self.k)


# The Fast Fourier Transform Algorithm
# Input: A, An array of integers of size n representing a polynomial
#        omega, a root of unity
# Output: [A(omega), A(omega^2), ..., A(omega^(n-1))]
# Complexity: O(n log n)
def FFT(A, omega):
    if omega == 1:
        return [sum(A)]
    B = [[],[]]
    i = 0
    for a in A:
        B[i % 2].append(a)
        i += 1
    o2 = omega**2
    C1 = FFT(B[0], o2)
    C2 = FFT(B[1], o2)
    C3 = [None] * omega.th
    for i in range(omega.th//2):
        C3[i] = C1[i] + omega**i * C2[i]
        C3[i+omega.th//2] = C1[i] - omega**i * C2[i]
    return C3

# The Fast Polynomial Multiplication Algorithm
# Input: A,B, two arrays of integers representing polynomials
# Output: Coefficient representation of AB
# Complexity: O(n log n)
def FPM(A, B):
    n = 1 << (len(A)+len(B)-2).bit_length()
    #print('n = ', n)
    o = NthRootOfUnity(n)
    #print('Nth Root of Unity = ', o)
    AT = FFT(A, o)      #FFT input: coefficient representation of polynomial, output: value representation of polynomial
    #print('Value representation of A:\n', AT)
    #print('len(AT) = ', len(AT))
    BT = FFT(B, o)
    #print('Value representation of B:\n', BT)
    #print('len(BT) = ', len(BT))
    C = [AT[i] * BT[i] for i in range(n)]
    print('Value representation of C:\n', C)
    #print('len(C) = ', len(C))
    #nm = (len(A) + len(B) - 1)
    D = [int((a/n).real) for a in FFT(C, o**-1)]
    print('Coefficient representation of C:\n', D)
    while True:
        if D[-1] != 0:
            return D
        else:
            del D[-1]

a = [1, 3, 1, 2, 2, 0, 5]           # (1 + 3x + x^2 + 2x^3 + 2x^4 + 5x^6)
b = [10, 0, 0, 0, 2, 7, 2, 0, 1]    # (10 + 2x^4 + 7x^5 + 2x^6 + x^8)
#5x^14 + 12x^12 + 37x^11 + 15x^10 + 21x^9 + 21x^8 + 17x^7 + 75x^6 + 13x^5 + 22x^4 + 20x^3 + 10x^2 + 30x^1 + 10x^0
print(FPM(a, b))

import cmath
"""
c1 = -3 + 1j
c2 = 11 + 7j
print(c1 * c2)

import cmath
import math
def omega(p, q):
   ''' The omega term in DFT and IDFT formulas'''
   return cmath.exp((2.0 * cmath.pi * 1j * q) / p)

print(1j**2)
print(omega(8, -4))
print(complex(3, -3))
z = complex(math.cos(math.pi / 6), math.sin(math.pi / 6))
print(z)
z = -0.5 + (math.sqrt(3)/2)*1j
print(z*z*z)

#n = 1<<(len(A)+len(B)-2).bit_length()

n = 1 << (5+4-2).bit_length()
print(n)
"""

#https://www.youtube.com/watch?v=Ty0JcR6Dvis    FFT example, unraveling the recursion
def fft2(P):
    # P: [p0, p1, ...pn-1]  coefficient representation of a polynomial
    n = len(P)      #n is a power of 2
    if n == 1:
        return P
    w = cmath.exp((2.0 * cmath.pi * 1j) / n)    #this is omega
    P_even = P[0::2]
    P_odd = P[1::2]
    y_even = fft2(P_even)
    y_odd = fft2(P_odd)
    y = [0] * n
    for j in range(n//2):
        y[j] = y_even[j] + w**j * y_odd[j]
        y[j + n//2] = y_even[j] - w**j * y_odd[j]
    return y


P = [5, 3, 2, 1]    #P(x) = 5 + 3x + 2x^2 + x^3
print(fft2(P))

def ifft2(P):       #this is not working for now
    # P: [p0, p1, ...pn-1]  coefficient representation of a polynomial
    n = len(P)      #n is a power of 2
    if n == 1:
        return P
    w = (cmath.exp((-2.0 * cmath.pi * 1j) / n)) / n   #this is omega
    P_even = P[0::2]
    P_odd = P[1::2]
    y_even = fft2(P_even)
    y_odd = fft2(P_odd)
    y = [0] * n
    for j in range(n//2):
        y[j] = y_even[j] + w**j * y_odd[j]
        y[j + n//2] = y_even[j] - w**j * y_odd[j]
    return y

def ifft(X):
    #The complex conjugate of a complex number is obtained by changing the sign of its imaginary part.
    result = fft2([x.conjugate() for x in X])
    return [x.conjugate()/len(X) for x in result]

def multiply_polynomials(p, q):
    x = len(p) + len(q)
    n = 1

    while n < x:
        n = n * 2

    for i in range(n - len(p)):
        p.append(0)
    for i in range(n - len(q)):
        q.append(0)

    pfft = fft2(p)
    qfft = fft2(q)

    c = []
    for i in range(n):
        c.append(pfft[i] * qfft[i])

    d = ifft(c)
    result = []
    for i in range(len(d)):
        result.append(round(d[i].real))

    string = ""
    n = len(result)
    for i in range(n):
        if result[n - 1 - i] != 0:
            if result[n - 1 - i] == 1:
                string += ("x^" + str(n - 1 - i) + " + ")
            elif (n - 1 - i) == 0:
                string += str(result[n - 1 - i])
            else:
                string += (str(result[n - 1 - i]) + "x^" + str(n - 1 - i) + " + ")
    #string = string[:-3]

    return string

r = multiply_polynomials(a, b)
print(r)

