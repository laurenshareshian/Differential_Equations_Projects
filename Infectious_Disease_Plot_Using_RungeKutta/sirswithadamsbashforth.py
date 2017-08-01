# SIRS model solved using Runge Kutta method
# the SIRS model models the spread of a disease through a population

import matplotlib.pyplot as plt

# differential equations
# disease parameters
beta = 0.6
lamb = 0.34
mu = 0.1
def f1(s, i): return -beta * s * i + mu * (1 - s - i)
def f2(s, i): return beta * s * i - lamb * i

# Runge Kutta method routine
def rk4(t0,s0,i0,tend,n):
    h = (tend - t0) / float(n)
    s = []
    i = []
    t = []
    s.append(s0)
    i.append(i0)
    t.append(t0)

    for m in range(3):
        ks1 = f1(s[m], i[m])
        ks2 = f1(s[m] + .5 * ks1 * h, i[m] + 0.5 * ks1 * h)
        ks3 = f1(s[m] + .5 * ks2 * h, i[m] + 0.5 * ks2 * h)
        ks4 = f1(s[m] + ks3 * h, i[m] + ks3 * h)

        ki1 = f2(s[m], i[m])
        ki2 = f2(s[m] + .5 * ki1 * h, i[m] + 0.5 * ki1 * h)
        ki3 = f2(s[m] + .5 * ki2 * h, i[m] + 0.5 * ki2 * h)
        ki4 = f2(s[m] + ki3 * h, i[m] + ki3 * h)

        t0 = t[m] + h
        s0 = s[m] + h / 6 * (ks1 + 2 * ks2 + 2 * ks3 + ks4)
        i0 = i[m] + h / 6 * (ki1 + 2 * ki2 + 2 * ki3 + ki4)

        t.append(t0)
        s.append(s0)
        i.append(i0)

    return t, s, i

# adams-bashforth
def adams_bash(t, s, i, tend, n):
    h = (tend - t[0]) / float(n)
    for m in range(4,n):
        t.append(t[m-1]+h)
        s.append(s[m-1]+h/24*(55*f1(s[m-1],i[m-1])-59*f1(s[m-2],i[m-2])+37*f1(s[m-3],i[m-3])-9*f1(s[m-4],i[m-4])))
        i.append(i[m-1]+h/24*(55*f2(s[m-1],i[m-1])-59*f2(s[m-2],i[m-2])+37*f2(s[m-3],i[m-3])-9*f2(s[m-4],i[m-4])))

    return t, s, i

def main():

    # initial conditions
    s0 = 0.9
    i0 = 0.1
    t0 = 0.0
    tend = 100
    n = 100


# use runge-kutta to get next 3 values
    t, s, i = rk4(t0,s0,i0,tend,n)

# use adams-bashforth to get the rest
    t, s, i = adams_bash(t,s,i, tend, n)

    # get recoverd population
    r = []
    for j in range(len(i)):
        r.append(1 - s[j] - i[j])

    for m in range(len(i)):
        print(t[m], s[m], i[m])

    # plot
    plt.figure(1)
    plt.subplot(211)
    plt.plot(s, i, 'o')
    plt.xlabel('susceptible')
    plt.ylabel('infected')
    plt.title('Susceptible versus Infected Population Percentages using RK4')
    plt.grid(True)
    plt.subplot(212)
    line_1, = plt.plot(t, s, label='susceptible')
    line_2, = plt.plot(t, i, label='infected')
    line_3, = plt.plot(t, r, label='recovered')
    plt.xlabel('days')
    plt.ylabel('percent')
    plt.legend([line_1, line_2, line_3], ['susceptible', 'infected', 'recovered'])
    plt.grid(True)
    plt.show()

if __name__ =='__main__':
    main()