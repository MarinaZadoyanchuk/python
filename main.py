from sympy import Function, Symbol, solve, diff, plot, integrate, Integral, pprint, solve_linear_system, Matrix
import matplotlib.pyplot as plt
import math

x = Symbol('x')
u = Function('u')

def metod_progonky(matrix):
    td = len(matrix[0])
    row = len(matrix)
    P = [0]*(row)
    Q = [0]*(row)
    xi = [0]*(row)
    P[0] = -matrix[0][1]/matrix[0][0]
    Q[0] = matrix[0][td-1]/matrix[0][0]
    for i in range(1, row):
        P[i] = matrix[i][i+1]/(-matrix[i][i]-matrix[i][i-1]*P[i-1])
        Q[i] = (matrix[i][i-1]*Q[i-1]-matrix[i][td-1])/(-matrix[i][i]-matrix[i][i-1]*P[i-1])
    xi[0] = (matrix[row-1][td-3]*Q[row-2]-matrix[row-1][td-1])/(-matrix[row-1][td-2]-matrix[row-1][td-3]*P[row-2])
    for j in range(1,row):
        xi[j] = P[row-j-1]*xi[j-1]+Q[row-j-1]
    return xi
    


def get_data():
    a = [5, 3, 2, 4]
    n = [7, 5, 2]
    k = [4, 3]
    b = [1, 4, 2]
    p = [1, 3]
    c = [2, 3, 3]
    q = [2, 2]
    d = [1, 3, 1]
    _a = 0
    _b = 1

    solution = a[0]*x**n[0]+a[1]*x**n[1]+a[2]*x**n[2]+a[3];

    alpha = a[0]*_a**n[0]+a[1]*_a**n[1]+a[2]*_a**n[2]+a[3];
    beta = (a[0]*n[0]*_a**(n[0]-1)+a[1]*n[1]*_a**(n[1]-1)+a[2]*n[2]*_a**(n[2]-1))
    gamma = (a[0]*_b**n[0]+a[1]*_b**n[1]+a[2]*_b**n[2]+a[3])
    delta = -(a[0]*n[0]*_b**(n[0]-1)+a[1]*n[1]*_b**(n[1]-1)+a[2]*n[2]*_b**(n[2]-1))


    kx = b[0]*x**k[0]+b[1]*x**k[1]+b[2]
    px = c[0]*x**p[0]+c[1]*x**p[1]+c[2]
    qx = d[0]*x**q[0]+d[1]*x**q[1]+d[2]

    # l_colloc = -(kx * u(x).diff(x)).diff(x) + px * u(x).diff(x) + qx * u(x)
    # f_colloc = (-kx * solution.diff(x, 2) -
    #      kx.diff(x) * solution.diff(x) +
    #      px * solution.diff(x) +
    #      qx * solution)

    l_ritz = -(kx * u(x).diff(x)).diff(x) + qx * u(x)
    f_ritz = (-kx * solution.diff(x, 2) -
         kx.diff(x) * solution.diff(x) +
         qx * solution)

    return {"k": kx, "q": qx,"p": px, "f":-f_ritz,"alpha": alpha, "beta": beta, "gamma": gamma, "delta": delta, "a": _a, "b": _b, 'solution': solution}

    

def integro_interpol(number_of_points):
    n = number_of_points
    data = get_data()
    k, q, p, f, a, b, alpha, beta, gamma, delta = (data[i] for i in ('k', 'q', 'p', 'f','a','b','alpha', 'beta', 'gamma', 'delta'))
    # k1 = kx
    # k2 = px - kx.diff(x)
    # mu = (1/k1)*math.e**integrate(k2/k1,x)
    # k = k1*mu
    # q = -qx*mu
    # f = fx*mu

    h = (b-a)*1.0/n
    beta = (beta*k.subs(x, a)*1.0) / alpha
    delta = (delta * k.subs(x, b) * 1.0) / gamma


    ai = [k.subs(x, a + h * i + h/2) for i in range(number_of_points)]


    y = [Symbol('y' + repr(i)) for i in range(number_of_points + 1)]
    system = []

    new_equation = [0] * (n + 1)
    new_equation[0] = (beta + 0.5 * h * q.subs(x, a) + ai[0] / h)
    new_equation[1] = -ai[0] / h
    new_equation.append(-0.5*h*f.subs(x,a));
    system.append(new_equation)
    for i in range(1, number_of_points):
        new_equation = [0] * (n + 1)
        new_equation[i - 1] = ai[i-1] / h ** 2
        new_equation[i] = -ai[i] / h ** 2 - ai[i-1]/h**2 - q.subs(x, i * h + a)
        new_equation[i + 1] = ai[i] / h ** 2
        new_equation.append(f.subs(x, i*h))
        system.append(new_equation)
    new_equation = [0] * (n + 1)
    new_equation[n-1] = -ai[n-1] / h
    new_equation[n] = (ai[n-1] / h + delta + 0.5 * q.subs(x, b) * h)
    new_equation.append(-0.5*h*f.subs(x,b))
    # print system
    system.append(new_equation)
    xi = metod_progonky(system)

    return ([a + i * h  for i in range(n + 1)], [xi[n-j] for j in range(n+1)])

number = 200
real_solution = get_data()['solution']
x_real = [n*1.0/number for n in range(number)]
y_real = [real_solution.subs(x,n) for n in x_real]
plt.plot(x_real, y_real, 'g--')

f = integro_interpol(number)
plt.plot(f[0], f[1], 'r-')
w = 32
print "xi".rjust(w), "yi".rjust(w),  "u(xi)".rjust(w), "yi-u(xi)".rjust(w)
for i in range(201):
    print `f[0][i]`.rjust(w), `f[1][i]`.rjust(w), `real_solution.subs(x, f[0][i])`.rjust(w), `(f[1][i]-real_solution.subs(x,f[0][i]))`.rjust(w)

plt.legend ( ("Tochnyj", "IIM") )
plt.show()
