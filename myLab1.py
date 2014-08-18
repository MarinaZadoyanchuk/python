from sympy import Function, Symbol, solve, diff, plot, integrate, Integral, pprint, solve_linear_system, Matrix
import matplotlib.pyplot as plt

x = Symbol('x')
u = Function('u')


def get_data(colloc):
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

    # alpha = [solution.subs(x, _a),
    #          solution.diff(x).subs(x, _a),
    #          solution.subs(x, _b),
    #          -solution.diff(x).subs(x, _b)
             #]
    alpha = a[0]*_a**n[0]+a[1]*_a**n[1]+a[2]*_a**n[2]+a[3];
    beta = a[0]*n[0]*_a**(n[0]-1)+a[1]*n[1]*_a**(n[1]-1)+a[2]*n[2]*_a**(n[2]-1)
    gamma = a[0]*_b**n[0]+a[1]*_b**n[1]+a[2]*_b**n[2]+a[3]
    delta = -(a[0]*n[0]*_b**(n[0]-1)+a[1]*n[1]*_b**(n[1]-1)+a[2]*n[2]*_b**(n[2]-1))

    # function_factory = lambda b, k:  b[0] * x ** k[0] + b[1] * x ** k[1] + b[2]

    kx = b[0]*x**k[0]+b[1]*x**k[1]+b[2]
    px = c[0]*x**p[0]+c[1]*x**p[1]+c[2]
    qx = d[0]*x**q[0]+d[1]*x**q[1]+d[2]

    l_colloc = -(kx * u(x).diff(x)).diff(x) + px * u(x).diff(x) + qx * u(x)
    f_colloc = (-kx * solution.diff(x, 2) -
         kx.diff(x) * solution.diff(x) +
         px * solution.diff(x) +
         qx * solution)

    l_ritz = -(kx * u(x).diff(x)).diff(x) + qx * u(x)
    f_ritz = (-kx * solution.diff(x, 2) -
         kx.diff(x) * solution.diff(x) +
         qx * solution)

    if colloc:
        return {"equation": l_colloc - f_colloc, "alpha": alpha, "beta": beta, "gamma": gamma, "delta": delta, "a": _a, "b": _b, 'solution': solution}
    else:
        return {"k": kx, "q": qx, "f":f_ritz,"alpha": alpha, "beta": beta, "gamma": gamma, "delta": delta, "a": _a, "b": _b, 'solution': solution}

def get_basis(alpha,beta, gamma, delta, a, b, n):
    A = Symbol('A')
    B = Symbol('B')
    phi = [(x-a)**2*(x -B), (x-b)**2*(x-A)]
    statements = [alpha * u(x).diff(x) - beta * u(x), gamma * u(x).diff(x) +delta * u(x)]
    eq = statements[0].subs(u(x), phi[1]).doit().subs(x, a)
    A_sol = solve(eq, A)[0]
    eq = statements[1].subs(u(x), phi[0]).doit().subs(x, b)
    B_sol = solve(eq, B)[0]
    basis = [phi[0].subs(B, B_sol), phi[1].subs(A, A_sol)]
    cur_phi = (x - a) ** 2 * (x - b) ** 2
    for i in range(n - 2):
        basis.append(cur_phi)
        cur_phi *= x - a
    return basis

def collocation_method(number_of_points):
    data = get_data(True)
    a, b = data['a'], data['b']
    alpha, beta, gamma, delta = data['alpha'], data['beta'], data['gamma'], data['delta']
    c = [Symbol('c' + repr(i)) for i in range(number_of_points)]
    basis = get_basis(alpha, beta, gamma, delta, a, b, number_of_points)
    u = sum([c[i] * basis[i] for i in range(number_of_points)])
    equation = data['equation'].subs(Function('u')(x), u).doit()
    points = [float(i + 1) / (number_of_points + 1) * (b - a) + a
              for i in range(number_of_points)]
    system = [equation.subs(x, i) for i in points]
    solutions = solve(system)
    u = u.subs(solutions.items())
    return u

def ritz_method(number_of_points):
    data = get_data(False)
    k, q, f,a,b, alpha, beta, gamma, delta = (data[i] for i in ('k', 'q', 'f','a','b','alpha', 'beta', 'gamma', 'delta'))
    phi = get_basis(alpha, beta, gamma, delta, a, b, number_of_points)

    # phi = [x**i for i in range(number_of_points)]
    beta = -beta* 1.0 / alpha * k.subs(x, a)
    delta = delta*1.0/ gamma * k.subs(x, b)
    g = lambda u, v: (integrate(k * u.diff(x) * v.diff(x) + q * u * v, (x, a,b)) 
        + beta * u.subs(x, a) * v.subs(x, a) + delta * u.subs(x, b) * v.subs(x, b) )
    l = lambda u: integrate(f * u, (x, a, b))
    c = [Symbol('c' + repr(i)) for i in range(number_of_points)]
    system = []
    for j in range(number_of_points):
        system.append([g(phi[i], phi[j]) for i in range(number_of_points)])
        system[-1].append(l(phi[j]))
    solution = solve_linear_system(Matrix(system), *c)
    print(solution.items());
    # for j in range(number_of_points):
    #     system.append(sum([c[i] * g(phi[i], phi[j]) for i in range(number_of_points)]) - l(phi[j]))
    # solution = solve(system)
    u = sum([c[i] * phi[i] for i in range(number_of_points)])
    u = u.subs(solution.items())
    return u

g = get_data(True)['solution']
f = ritz_method(20)
# print(f);
c = collocation_method(5);
print(c);
# f = collocation_method(5)
# pylab.plot(f,(x, 0, 1))
# plt.plot(g)
# pylab.plot(c, (x, 0, 1))
x1 = range (0,200)
data1 = [ 1-n/200.0 for n in x1]
dataf = [f.subs(x,n) for n in data1]
datag = [g.subs(x,n) for n in data1]
datac = [c.subs(x,n) for n in data1]
plt.plot(data1, dataf, 'r-')
plt.plot(data1, datag, 'g--')
plt.plot(data1, datac, 'b')
plt.legend ( ("Ritza", "Kolokacii", "Tochnyj") )
plt.show()
