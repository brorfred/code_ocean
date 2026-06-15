import numpy as np
import matplotlib.pyplot as plt
import math
NT = 200
T = np.linspace(0,50,NT)

Z = 28
w = 15

def fN(x, w=15):
    return (1 - ( (x-Z)/w)**2 )*0.59*np.exp(0.0633*x)

def blackford(x, Q10=2):
    return Q10**((x-10)/10) - Q10**((x-Z)/3)

def blackford_derivative(T, Q10=2):
    return (np.log(Q10) / 10.0) * Q10**((T - 10) / 10) - (np.log(Q10) / 3) * Q10**((T - Z) / 3)


#plt.plot(T,fN(T))
#plt.xlim(0,50)
#plt.ylim(0,6);

def blackford_noise(T, Q10=2):

    blackfordm = np.zeros(NT)
    blackfordstd = np.zeros(NT)
    fNm = np.zeros(NT)
    fNstd = np.zeros(NT)
    Tstd = 5

    for y in range(0,NT):
        Trand = np.random.normal(loc=y*50/NT, scale=Tstd, size=200000)
        GRB =  blackford(Trand)
        #blackfordm[y] = np.mean( GRB[GRB > 0] )
        blackfordm[y] = np.mean( GRB )
        #blackfordstd[y] = np.std( blackford(Trand) )
        GRN =  fN(Trand)
        fNm[y] = np.mean( GRN )
        #fNstd[y] = np.std( fN(Trand) )
    return blackfordm
def plot(ax=None):
    if ax is None:
        plt.figure(figsize=(8,5))
        ax = plt.gca()
    #plt.plot(T,fN(T))
    ax.plot(T,blackford(T), lw = 3, label = 'Blackford')
    ax.plot(T,10*blackford_derivative(T), lw = 3, label = 'Blackford Derivative*10')
    ax.plot(T,blackford_noise(T), lw = 3, c = 'red', label = 'Blackford with noise')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Growth rate')
    ax.legend()
    ax.set_xlim(0,40)
    ax.set_ylim(0,4);
    plt.savefig("figs/Blackford.pdf")
    if ax is None:
        plt.show()
