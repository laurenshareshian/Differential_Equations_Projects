# SIRS model solved using Euler's method
# the SIRS model models the spread of a disease through a population

import matplotlib.pyplot as plt

# differential equations
# disease parameters
beta = 0.6
lamb = 0.34
mu = 0.1
def f1(s, i): return -beta * s * i + mu * (1 - s - i)
def f2(s, i): return beta * s * i - lamb * i

# Euler's method routine
def euler(t0,s0,i0,tend,n):
    h = (tend - t0) / float(n)
    s = []
    i = []
    t = []
    s.append(s0)
    i.append(i0)
    t.append(t0)

    for m in range(n):

        t0 = t0 + h
        stemp = s0
        s0 = s0 + h*f1(s0, i0)
        i0 = i0 + h*f2(stemp, i0)

        t.append(t0)
        s.append(s0)
        i.append(i0)

    return t, s, i


def main():

    # initial conditions
    s0 = 0.9
    i0 = 0.1
    t0 = 0.0
    tend = 100
    n = 100

    t, s, i = euler(t0,s0,i0,tend,n)

    # get recoverd population
    r = []
    for j in range(len(i)):
        r.append(1 - s[j] - i[j])

    # plot
    plt.figure(1)
    plt.subplot(211)
    plt.plot(s, i, 'o')
    plt.xlabel('susceptible')
    plt.ylabel('infected')
    plt.title('Susceptible versus Infected Population Percentages using Eulers method')
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