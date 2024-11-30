import numpy as np
from matplotlib import pyplot as plt

def step(x, yl, yr):
    y = []
    for i in x:
        if i <= 0:
            y.append(yl)
        else:
            y.append(yr)
    return np.array(y)

def mnk (u_prev, alpha):
    return np.dot(u_prev, alpha)

def first_order(num):
    if num == 1:
        return (0.5, 0)
    elif num == 2:
        return (2/3, 0)
    elif num == 3:
        return (0, 1/3)
    elif num == 4:
        return (0, 0.2)

def second_order(num):
    if num == "min_osc":
        return (77/156, -5/156)
    elif num == 1:
        return (0, 1/15)
    elif num == 2:
        return (0.6, -4/75)

def third_order():
    return (192/531, -1/177)

def alpha(a00, a11):
    a21 = -1 + 2*a00 + 5*a11
    a10 = 2 - 3*a00 - 6*a11
    return(np.array([a21, a10, a00, a11]))

def analitical(c, rho, t, x, ul, ur, pl, pr):
    u_an = np.zeros((201), np.double)
    p_an = np.zeros((201), np.double)
    for i in range(0, len(x)-1):
        if (x[i] < -c*t):
            u_an[i] = ul
            p_an[i] = pl
        elif (x[i] > c*t):
            u_an[i] = ur
            p_an[i] = pr
        else:
            u_an[i] = (ul + ur)/2 - (pr-pl)/(2*rho*c)
            p_an[i] = (pl + pr)/2 - (ur - ul) * rho*c / 2
    return(u_an, p_an)

def ordinary_scheme(u, inv):
    #a00, a11 = first_order(1)
    a00, a11 = second_order("min_osc")
    for i in range(1, len(u)-1):
        if (inv == "y"):
            u[i+1,0] = u[i,0];
            u[i+1,1] = u[i,1];
            for j in range(2, len(x)-1):
                u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha(a00, a11))
            u[i+1, len(x)-1] =  u[i, len(x)-1]
            
        elif (inv == "z"):
            u[i+1,0] = u[i,0];
            for j in range(1, len(x)-2):
                u[i+1,j] = mnk(np.array([u[i-1,j+2], u[i,j+1], u[i,j], u[i-1,j-1]]), alpha(a00, a11))
            u[i+1, len(x)-1] = u[i, len(x)-1];
            u[i+1, len(x)-2] = u[i, len(x)-2];

def hybrid_scheme(u, inv):
    a00, a11 = second_order(1)
    a00_, a11_ = third_order()
    alpha1 = alpha(a00, a11)
    alpha2 = alpha(a00_, a11_)
    for i in range(1, len(u)-1):
        if (inv == "y"):
            u[i+1,0] = u[i, 0];
            u[i+1,1] = u[i, 1];
            for j in range(2, len(x)-1):
                u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha1)
                if not (u[i+1,j]>=min(u[i,j-1], u[i,j+1]) and u[i+1,j]<=max(u[i,j-1], u[i,j+1])):
                    u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha2)
            u[i+1, len(x)-1] = u[i, len(x)-1];

        elif (inv == "z"):
            u[i+1,0] = u[i,0];
            for j in range(1, len(x)-2):
                u[i+1,j] = mnk(np.array([u[i-1,j+2], u[i,j+1], u[i,j], u[i-1,j-1]]), alpha1)
                if not (u[i+1,j]>=min(u[i,j-1], u[i,j+1]) and u[i+1,j]<=max(u[i,j-1], u[i,j+1])):
                    u[i+1,j] = mnk(np.array([u[i-1,j+2], u[i,j+1], u[i,j], u[i-1,j-1]]), alpha2)
            u[i+1, len(x)-1] = u[i, len(x)-1]
            u[i+1, len(x)-2] = u[i, len(x)-2]

courant = 1/2
h = 0.01
tau = courant*h

x = np.linspace(-2, 2, 201)
t = np.linspace(0, 100*tau, 100)
u = np.zeros((100, 201), np.double)
p = np.zeros((100, 201), np.double)

ul = 1.0
ur = 0.0
pl = 5.0
pr = 2.0

rho = 0.25
c = 2

u[0] = step(x, ul, ur)
u[1] = step(x-t[1], ul, ur)

p[0] = step(x, pl, pr)
p[1] = step(x-t[1], pl, pr)

y = np.zeros((100, 201), np.double)
z = np.zeros((100, 201), np.double)


y = u + p / (rho*c)
z = u - p / (rho*c)

#ordinary_scheme(y, "y")
#ordinary_scheme(z, "z")

hybrid_scheme(y, "y")
hybrid_scheme(z, "z")

u = (y + z)/2
p = (y - z)*rho*c/2 

u_an, p_an = analitical(c, rho, t[len(t) - 1], x, ul, ur, pl, pr)

fig, axs = plt.subplots(2, 1, sharex=True)
axs[0].plot(x, u[len(u)-1])
axs[0].plot(x, u_an)
axs[0].set_title("u")
axs[1].plot(x, p[len(p)-1])
axs[1].plot(x, p_an)
axs[1].set_title("p")

plt.show()
