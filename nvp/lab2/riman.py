import numpy as np
from math import sqrt
import csv
from matplotlib import pyplot as plt

def VFromW(W):
    V = []
    for j in range(len(W)):
            V.append([W[j][0], W[j][1]/W[j][0], 
                      W[j][2]/W[j][0]])
    return V

def plot(X, V, gamma):
    rho = []
    u = []
    p = []
    e = []
    for x in V:
        rho.append(x[0])
        u.append(x[1])
        e.append(x[2])
        p.append((gamma-1)*x[0]*x[2])
    ax1 = plt.subplot(2,2,1)
    ax1.scatter(X, rho)
    ax1.set_title('rho')
    ax1.spines['left'].set_position('center')
    plt.grid()

    ax2 = plt.subplot(2,2,2)
    ax2.scatter(X, u)
    ax2.set_title('u')
    ax2.spines['left'].set_position('center')
    plt.grid()

    ax3 = plt.subplot(2,2,3)
    ax3.scatter(X, e)
    ax3.set_title('e')
    ax3.spines['left'].set_position('center')
    plt.grid()

    ax4 = plt.subplot(2,2,4)
    ax4.scatter(X, p)
    ax4.set_title('p')
    ax4.spines['left'].set_position('center')
    plt.grid()
    plt.show()

    

def Omega(w, gamma):
    u = w[1]/w[0]
    c = sqrt(w[2]*(gamma-1)*gamma/w[0])
    return np.array([[-u*c, c, gamma-1], 
                    [-c**2, 0, gamma-1], [u*c, -c, gamma-1]])

def Lambda(w, gamma):
    u = w[1]/w[0]
    c = sqrt(w[2]*(gamma-1)*gamma/w[0])
    lambdas = np.array([u+c, u, u-c])
    return(np.diag(lambdas), np.max(np.abs(lambdas)))

def KIR(gamma, u0, rho0, p0, T, nodes, h, moment, x):
    n = int((nodes+0)/2)
    a = np.array([rho0[0], rho0[0]*u0[0], p0[0]/(gamma-1)])
    b = np.array([rho0[1], rho0[1]*u0[1], p0[1]/(gamma-1)])
    wa0 = np.full((n, 3), a)
    wb0 = np.full((n, 3), b)
    w0 = np.concatenate((wa0, wb0), axis=0)

    W = [w0]
    V = []
    
    t = 0
    tau = 1e-2
    
    while (t<=T):
        t += tau
        w = []
        for i in range(1, nodes-1):
            om = Omega(w0[i], gamma)
            om_1 = np.linalg.inv(om)
            l, lm = Lambda(w0[i], gamma)
            if (tau*lm/h) > 1:
                while(tau*lm/h > 1):
                    tau = tau/10
            a = l.dot(om)
            A = om_1.dot(a)
            a = np.abs(l).dot(om)
            A_ = om_1.dot(a)
            w_i = (w0[i] - tau*A.dot((w0[i+1] - w0[i-1])/(2*h))
                 + tau*A_.dot((w0[i+1] - 2*w0[i] + w0[i-1]))/(2*h))
            w.append(w_i)
        w.append(w_i)
        w = [w[0]] + w
        w = np.array(w)
        w0 = w
        W.append(w)
        if(abs(t - moment) < tau):
            plot(x, VFromW(w), gamma)
            print(t)

    for i in range(len(W)):
        v = VFromW(W[i])
        V.append(v)
    
    return np.array(V)

gamma = 5/3
u0 = [0, 0]
rho0 = [13, 1.3]
p0 = [10*101325, 101325]

T = 0.02
L = 10
nodes = 100
moment = 0.015
h = L/(nodes-0)
x = np.linspace(-L, L, nodes)
xa = np.linspace(-L, 0, int((nodes+0)/2))
xb = np.linspace(0, L, int((nodes+0)/2))

V = KIR(gamma, u0, rho0, p0, T, nodes, h, moment, x)

