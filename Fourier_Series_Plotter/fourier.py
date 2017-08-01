from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np
#this program graphs the fourier approximation and spectrum plot of a periodic function

def fourier_coefficients(f, p, n):
    a_list = []
    b_list = []
    energy_list = []
    a0 = 1 / p * integrate.quad(f, -p/2, p/2)[0]
    energy_list.append(2 * a0**2)

    for k in range(1,n+1):
        a = lambda x: np.cos(k*x)
        b = lambda x: np.sin(k*x)
        a_prod = lambda x: f(x)*a(x)
        b_prod = lambda x: f(x)*b(x)
        ak = 2/p*integrate.quad(a_prod, -p/2, p/2)[0]
        bk = 2/p * integrate.quad(b_prod, -p/2, p/2)[0]
        a_list.append(ak)
        b_list.append(bk)
        energy_list.append(ak**2+bk**2)
    return a0, a_list, b_list, energy_list

def fourier_function(x, n, a0, a_list, b_list):
   func = np.array([a_list[k-1]*np.cos(k*x)+b_list[k-1]*np.sin(k*x) for k in range(1,n+1)])
   return a0+func.sum()

# What p-periodic function do you want to find the Fourier Coefficients for?
def f(x):
    p=2*np.pi
    if -p/2 <= x < 0:
        return 0
    if 0 <= x <= p/2:
        return 1

def main():
    #How many harmonics do you want?
    n=3

    #What is the period of your function?
    p=2*np.pi

    a0, a_list, b_list, energy_list = fourier_coefficients(f, p, n)

    time = np.linspace(0, 10, 101)
    function_values = np.array([f(x) for x in time])
    fourier_approx = np.array([fourier_function(x,n, a0, a_list, b_list) for x in time])

    plt.plot(time,fourier_approx)
    plt.plot(time, function_values)
    plt.title('Fourier Approximation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

    plt.plot([0,0], [0, energy_list[0]], color = 'black')
    for k in range(1,n+1):
        plt.plot([k,k], [0, energy_list[k]], color = 'black')
    plt.plot([0,0], [0, energy_list[0]], color = 'black')
    plt.title('Energy Spectrum Plot')
    plt.xlabel('k')
    plt.ylabel('A_k^2')
    plt.show()



if __name__ == "__main__":
    main()

