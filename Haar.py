import math
import random
import matplotlib.pyplot as plt


def plot_f(f):
    # plot
    N = 1000
    x = [0] * N
    y = [0] * N
    for i in range(N):
        x[i] = (1 / N) * i
        y[i] = f(x[i])
    plt.plot(x, y)



# plot the Gaussian Distribution and its samples
N = 10000
x = [0.]*N
y = [0.]*N
z = [0]*N
for i in range(N):
    x[i] = i/N
    y[i] = (1/((2*math.pi*0.1**2)**0.5))*math.exp(-((x[i] - 0.5)**2)/(2*0.1**2))
    sample = random.gauss(0.5, 0.1)
    box = round(N*round(sample, 2))
    z[box] += (1 / N) * 100

# plt.plot(x, y)
# plt.plot(x, z, 'r-')
# plt.show()

# now for the Haar basis
def phi(j, k):
    return lambda x : (2**(j/2))*(1 if ((2**j)*x - k) >= 0 and ((2**j)*x-k) <= 1 else 0)

def psi(j, k):
    return lambda x : (2**(j/2))*(1 if ((2**j)*x - k ) >= 0 and ((2**j)*x-k) <= 0.5 else (-1 if 1 >= ((2 ** j) * x - k) >= 0.5 else 0))



# Q = [0,1], l=1
J = 8
j0 = 5
resolution = 2**(-J)
k = [0]*(2**j0)
Phi = []
for i in range(2**j0):
    k[i] = i
    Phi.append(phi(j0, k[i]))
Psi = []
for j in range(j0, J):
    Psi.append([])
    for i in range(2**j):
        k_i = i
        Psi[j-j0].append(psi(j,k_i))


# now we will approximate the sin function
s = lambda x : 100*math.pi*math.sin(math.pi*x)
# sample = lambda a, b: 0.5*math.cos((math.pi)*a) - 0.5*math.cos((math.pi)*b)
# sample = lambda a, b: 0.5*(s(a) + s(b))/((2**(j0/2)))
sample = lambda a, b: 0.5*(s(a) + s(b))
N = 10000
resolution = 500
samples = []
for i in range(resolution):
    dx = 1/resolution
    tot = 0
    for k in range(resolution):
        tot += sample(k*dx, (k+1)*dx)
    for num in range(round((sample(i*dx, (i+1)*dx)/tot)*N)):
        samples.append(i*dx + (0.5)*dx)

def indicator(x, a, b):
    return 1 if a <= x <= b else 0

alpha = []
for i in range(2**j0):
    # alpha.append((1/N) * sum([Phi[i](x) for x in samples]))
    alpha.append((1/N) * sum([Phi[i](x) for x in samples]))

print(alpha)
def f_a(x):
    tot = 0
    for i in range(2**j0):
        tot += alpha[i]*Phi[i](x)
    return tot
plot_f(f_a)
plot_f(lambda x : 0.5
                  *math.pi*math.sin(math.pi*x))
plt.show()



def f_aprx(x):
    res = 10
    fx = 0
    a = 0
    tot = 0
    for i in range(2**j0):
        box = (i * (1 / (2 ** j0)), (i + 1) * (1 / (2 ** j0)))
        tot += sample(box[0], box[1])
    for i in range(2**j0):
        box = (i * (1 / (2 ** j0)), (i + 1) * (1 / (2 ** j0)))
        a += Phi[i](x) * sample(box[0], box[1])/tot
        print(f"added = {Phi[i](x) * sample(box[0], box[1])/tot}")
        print(f"phi coeff = {(2**(j0/2))}")
    # for i in range(2**j0):
    #     box = (i*(1/(2**j0)), (i+1)*(1/(2**j0)))
    #     # a += Phi[i](x) * alpha[i]
    #     a += Phi[i](x) * sample(box[0], box[1])
    #     tot += sample(box[0], box[1])
    # a = a / tot
        # print(f"alpha = {alpha[i]} sin(x) = {0.5*math.sin(math.pi*(i*(1/(2**j0) + (0.5/(2**0)))))}")
    fx += a
    #
    # for j in  range(j0, j0+1):
    #     for i in range(2**j):
    #         b = 0
    #         left = i * (1 / (2 ** j))
    #         for dx in range(res):
    #             width = (1/res)*(1 / (2 ** j))
    #             b += Psi[j - j0][i](left + (dx*width)) * (sample(left, left + width))
    #         print(f"b= {b}")
    #         fx += b*Psi[j-j0][i](x)
    # fx += b
    return fx

# plot_f(f_aprx)
# plot_f(lambda x : 0.5*math.pi*math.sin(math.pi*x))
# def dist(x):
#     dx = 1/20
#     i = math.floor(x/dx)
#     return sample(i*dx, (i+1)*dx)
# plot_f(dist)
# plt.show()

# v0 = lambda x : 1
# v1 = lambda x : x
# v2 = lambda x : 0 if x < 0.5 else (x - 0.5)
#
# dx = 0.001
# f0 = lambda x : 1
# f1 = lambda x : (v1(x) - f0(x)*(sum([dx*f0(0.001*t)*v1(0.001*t) for t in range(1000)])))/(-sum([dx*v1(0.001*t)*v1(0.001*t) for t in range(1000)])**0.5)
#
#
# print(f"{(-sum([dx*v1(0.001*t)*v1(0.001*t) for t in range(1000)]))}")
# plot_f(v1)
#
# plt.figure()
#
# plot_f(f1)
# plot_f(f_test)
# plot_f(f_test2)
# plt.show()



# # plot
# N = 1000
# x = [0]*N
# y = [0]*N
# for r in range(5,10):
#     for i in range(N):
#         x[i] = (1/N)*i
#         y[i] = Psi[r](x[i])
#     plt.plot(x,y)
# plt.show()




