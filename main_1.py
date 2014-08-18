from sympy import Function, Symbol, solve, diff, plot, integrate, Integral, pprint, solve_linear_system, Matrix
import matplotlib.pyplot as plt
import math

x = Symbol('x')
u = Function('u')


a = 0 
b = 0.5
u0 = 291
u_sr = 278
lamb = 0.885
c = 882
ro = 2100
time = 3600

sigma = 0.5
gamma = 8.15
gamma1 = gamma*(b-a)/lamb

tau = 0.0001
N = 50
h = 1.0/N
print tau, 0.5 * h ** 2

def outsizing_u(u):
    return (u - u_sr) * 1.0 / abs(u_sr)

def sizing_u(u):
    return u * abs(u_sr) + u_sr

def outsizing_time(t):
    return t * (lamb/(c*ro) / (b-a) ** 2)

def sizing_time(t):
    return t / (lamb/(c*ro) / (b-a) ** 2)

def metod_progonky(a,b,c,f):
    row = len(b)
    P = [0]*(row)
    Q = [0]*(row)
    xi = [0]*(row)
    yi = [0]*(row)
    P[0] = -c[0]*1.0/b[0]
    Q[0] = f[0]*1.0/b[0]
    for i in range(1, row-1):
        P[i] = c[i]*1.0/(-b[i]-a[i-1]*P[i-1])
        Q[i] = (a[i-1]*Q[i-1]-f[i])*1.0/(-b[i]-a[i-1]*P[i-1])
    xi[0] = (a[row-2]*Q[row-2]-f[row-1])*1.0/(-b[row-1]-a[row-2]*P[row-2])
    for j in range(1,row):
        xi[j] = P[row-j-1]*xi[j-1]+Q[row-j-1]
    for j in range(row):
        yi[j] = xi[row-j-1]
    return yi

    
 # a=b[:]

def main_method(number_of_points):
    n = number_of_points
    h = 1.0/n

    v0 = outsizing_u(u0)
    
    new_equation_copy = []   
    new_equation = [0, 0, 0, 0]
    new_equation_copy.append([v0]*(n+1))

    


    for k in range(0, int(outsizing_time(time)/tau)): 
        new_equation[0] = []
        new_equation[1] = []
        new_equation[2] = []
        new_equation[3] = [] 

        new_equation[1].append(-sigma*1.0/h)
        new_equation[2].append(sigma*1.0/h)
        new_equation[3].append(  -(1-sigma)*(-(new_equation_copy[k][1]-new_equation_copy[k][0])*1.0/h)+gamma1*new_equation_copy[k][0] )

        for i in range(1, n):
            new_equation[0].append(sigma*1.0/h**2)
            new_equation[1].append(-1.0/tau-2*sigma*1.0/h**2)
            new_equation[2].append(sigma*1.0/h**2)
            new_equation[3].append(  -(1-sigma)*(new_equation_copy[k][i+1]-2*new_equation_copy[k][i]+new_equation_copy[k][i-1])/h**2-new_equation_copy[k][i]/tau  )

        new_equation[0].append(sigma*1.0/h)  
        new_equation[1].append(sigma*1.0/h)
        new_equation[3].append( (1-sigma)*( -(new_equation_copy[k][n]-new_equation_copy[k][n-1])*1.0/h)+gamma1*new_equation_copy[k][n]    )

        xi = metod_progonky(new_equation[0], new_equation[1], new_equation[2], new_equation[3])
        new_equation_copy.append(xi)

    return new_equation_copy
    # return ([a + i * h  for i in range(n + 1)], new_equation_copy[j] for j in range(len(new_equation_copy)))
n = 100
real_h = (b-a)*1.0/n
# print metod_progonky([1, 1], [2,1,1], [1,1], [4,4,3]);
result = main_method(n)
for i in range(0, len(result)):
        plt.plot([a + j * real_h for j in range(n+1)], [sizing_u(j) for j in result[i]], "g--")
plt.show()
