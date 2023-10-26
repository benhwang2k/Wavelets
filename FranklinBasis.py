import copy
import math
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from functools import partial

"""
This file will contain the necessary tools to approximate functions using the Franklin system.
In extension, I look to develop the tools to construct an empirical version distribution of some data
in the Franklin space. The special thing about this empirical distribution is that it will be continuous
instead of the classic dirac delta version which is not a pdf and the haar (histogram) version which is 
discontinuous at the boundary of the boxes.
"""


def plot_f(f):
    x = [0] * 1001
    y = [0] * 1001
    for i in range(1001):
        x[i] = i * (1 / 1000)
        y[i] = f(x[i])
    plt.plot(x, y)


class FranklinFunc:
    """
    Since each function of the Franklin basis is piecewise linear with linear pieces between each point (i/(2**J))
    i in [0, (2**J)], we can describe each function completely with only the function values at each of those points.
    """

    def __init__(self, J):
        """
        :param J: J describes the resolution (1 / (2**J))
        """
        self.J = J
        self.dx = (1 / (2 ** J))
        self.N = 2 ** J + 1  # this is the number of points (including the extra endpoint at x = 1)
        self.f = [1.] * self.N

    def integrate(self):
        tot = 0
        for i in range(self.N - 1):
            tot += self.dx * (self.f[i] + self.f[i + 1]) / 2
        return tot

    def inner_product(self, v: 'FranklinFunc'):
        tot = 0
        for i in range(self.N - 1):
            m1 = (self.f[i + 1] - self.f[i]) / self.dx
            m2 = (v.f[i + 1] - v.f[i]) / self.dx
            a = self.f[i]
            c = v.f[i]
            tot += a * c * self.dx + 0.5 * (a * m2 + c * m1) * (self.dx ** 2) + (1 / 3.) * (m1 * m2) * (self.dx ** 3)
        return tot

    def normalize(self):
        factor = self.inner_product(self)
        factor = factor ** 0.5
        for i in range(len(self.f)):
            self.f[i] = self.f[i] / factor

    def plot(self):
        x = [0] * self.N
        for i in range(self.N):
            x[i] = i * self.dx
        plt.plot(x, self.f)

    def mult_function(self, y):
        for i in range(self.N):
            self.f[i] = self.f[i] * y(i * self.dx)

    def mult_const(self, c):
        for i in range(len(self.f)):
            self.f[i] = self.f[i] * c

    def mult(self, g: 'FranklinFunc'):
        for i in range(len(self.f)):
            self.f[i] = self.f[i] * g.f[i]

    def add(self, g: 'FranklinFunc'):
        for i in range(len(self.f)):
            self.f[i] = self.f[i] + g.f[i]

    def subtract(self, g: 'FranklinFunc'):
        for i in range(len(self.f)):
            self.f[i] = self.f[i] - g.f[i]

    def get_fcn(self):
        def f(x):
            if x == 1:
                return self.f[self.N-1]
            else:
                i = math.floor((x / 1) * (2 ** self.J))
                pcnt = (x - i * self.dx) / (self.dx)
                return self.f[i] + ((self.f[i + 1] - self.f[i]) / self.dx) * pcnt * self.dx
        return f




def orthogonalize(v: FranklinFunc, basis):
    f = FranklinFunc(v.J)
    f.mult(v)
    for base in basis:
        temp = FranklinFunc(v.J)
        temp.mult(base)
        proj = temp.inner_product(v)
        temp.mult_const(proj)
        f.subtract(temp)
    f.normalize()
    return f


class FranklinBasis:
    """
    This class will give the set of functions that is used to construct the basis of functions for the Franklin system.
    The chosen resolution parameter J will define the granularity of the approximations possible from this basis.
    """

    def __init__(self, J):
        """
        :param J: Integer J defines the resolution of the approximation. 1/(2**J) is the width of resolution grains.
        """
        self.J = J
        self.dx = (1 / (2 ** J))
        self.N = 2 ** J

        """ initialize the basis with the constant function """
        self.raw_basis = [FranklinFunc(J)]
        self.basis = [FranklinFunc(J)]
        """ add the first basis function f(x) = x after orthogonalization"""
        fnext = FranklinFunc(J)
        fnext.mult_function(lambda x: x)
        fnext = orthogonalize(fnext, self.basis)
        self.basis.append(fnext)
        fnext = FranklinFunc(J)
        fnext.mult_function(lambda x: x)
        self.raw_basis.append(fnext)
        """ add the rest of the functions up to the desired resolution"""
        j = 1
        while j <= J:
            num = 1
            while num < (2 ** j):
                an = num / (2 ** j)
                fnext = FranklinFunc(J)
                fnext.mult_function(lambda x: (x - an) if x >= an else 0)
                self.raw_basis.append(fnext)
                fnext = FranklinFunc(J)
                fnext.mult_function(lambda x: (x - an) if x >= an else 0)
                fnext = orthogonalize(fnext, self.basis)
                self.basis.append(fnext)

                num += 2
            j += 1

    def approximate_function(self, fnc):
        a = [1.]*self.N         # the coefficient weights for the approximation
        f_tot = FranklinFunc(self.J)    # the approximation
        f_tot.mult_function(lambda x: 0)    # initialize to 0

        """ Obtain the weights and add the result to the total"""
        for i in range(self.N):
            fi = FranklinFunc(self.J)
            fi.mult(self.basis[i])
            fi_func = fi.get_fcn()
            (a[i], _) = integrate.quad(lambda x: fi_func(x) * fnc(x), 0, 1)
            fi.mult_const(a[i])
            f_tot.add(fi)
        print(f"the coefficients are {a}")
        return f_tot

approximator = FranklinBasis(4)
f_true = lambda x: math.sin(2*math.pi*x)
approx = approximator.approximate_function(f_true)
approx.plot()
plt.show()




# # Test functions
# f0 = lambda x: 1
#
#
# def f1(x):
#     return (3. ** 0.5) * (2 * x - 1)
#
#
# def f2(x):
#     if x <= 0.5:
#         return (3. ** 0.5) * (1 - 4 * x)
#     else:
#         return (3. ** 0.5) * (4 * x - 3)
#
#
# def f3(x):
#     if x <= 0.25:
#         return (3. ** 0.5) * (1 / (19)) * (10 - 76 * x)
#     elif 0.25 < x <= 0.5:
#         return (3. ** 0.5) * (1 / (19)) * (-22 + 52 * x)
#     else:
#         return (3. ** 0.5) * (1 / (19)) * (10 - 12 * x)
#
#
# def f3test(x):
#     v = (x - 1 / 4) if x >= 1 / 4 else 0
#     v -= (9 / 32)
#     v -= ((9 * 3 ** 0.5) / (64)) * f1(x)
#     v -= (1 / (16 * 3 ** 0.5)) * f2(x)
#     return v
#
#
# def f3normalized(x):
#     (intg, _) = integrate.quad(lambda x: f3test(x) ** 2, 0, 1)
#     return f3test(x) / (intg ** 0.5)


