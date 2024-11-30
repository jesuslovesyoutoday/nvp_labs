import numpy as np
from matplotlib import pyplot as plt

def u0(x):
    y = np.zeros((201), np.double)
    for i in range(len(x)):
        if (x[i] >= 0.4 and x[i] <= 0.6):
            y[i] = np.sqrt(1 - 100*np.power((x[i] - 0.5), 2))
    return y

def mnk (u_prev, alpha):
    return np.dot(u_prev, alpha)

def alpha(a00, a11):
    a21 = -1 + 2*a00 + 5*a11
    a10 = 2 - 3*a00 - 6*a11
    return(np.array([a21, a10, a00, a11]))

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

def ordinary_scheme(u):
    for i in range(1, len(u)-1):
        u[i+1,0] = 0;
        u[i+1,1] = 0;
        for j in range(2, len(x)-1):
            a00, a11 = first_order(1)
            u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha(a00, a11))

def hybrid(u):
    for i in range(1, len(u)-1):
        u[i+1,0] = 0;
        u[i+1,1] = 0;
        for j in range(2, len(x)-1):
            a00, a11 = second_order(1)
            a00_, a11_ = second_order(2)
            alpha1 = alpha(a00, a11)
            alpha2 = alpha(a00_, a11_)
            u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha1)
            if not (u[i+1,j]>=min(u[i,j-1], u[i,j+1]) and u[i+1,j]<=max(u[i,j-1], u[i,j+1])):
                u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha2)

def three_hybrid(u):
    for i in range(1, len(u)-1):
        u[i+1,0] = 0;
        u[i+1,1] = 0;
        for j in range(2, len(x)-1):
            a00, a11 = second_order(1)
            a00_, a11_ = second_order(2)
            a00__, a11__ = third_order()
            alpha1 = alpha(a00, a11)
            alpha2 = alpha(a00_, a11_)
            alpha3 = alpha(a00__, a11__)
            u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha1)
            if not (u[i+1,j]>=min(u[i,j-1], u[i,j+1]) and u[i+1,j]<=max(u[i,j-1], u[i,j+1])):
                u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha2)
            if not (u[i+1,j]>=min(u[i,j-1], u[i,j+1]) and u[i+1,j]<=max(u[i,j-1], u[i,j+1])):
                u[i+1,j] = mnk(np.array([u[i-1,j-2], u[i,j-1], u[i,j], u[i-1,j+1]]), alpha3)

courant = 1/2
h = 0.01
tau = courant*h

x = np.linspace(0, 2, 201)
t = np.linspace(0, 100*courant*h, 100)
u = np.zeros((100, 201), np.double)
u[0] = u0(x)
u[1] = u0(x-t[1])

#ordinary_scheme(u)
three_hybrid(u)
#hybrid(u)
   
u_an = u0(x - t[len(t)-1])

print(np.max(np.abs((u[len(u)-1]-u_an))))

plt.plot(x, u_an, label='analitical')
plt.plot(x, u[len(u)-1], label='numerical')
plt.title('hybrid 2+2+3')
plt.legend()
plt.show()

