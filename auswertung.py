import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

# Konstanten

l = 1.54286

Ndiamant = np.array([3, 8, 11, 16, 19, 24, 27, 32, 35, 40, 43, 48, 51])

Ncaesiumchlorid = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26])

#########################
#                       #
#        Probe 1        #
#                       #
#########################

######################
# Einlesen der Daten #
######################

r , rerr = np.genfromtxt('werte1.txt', unpack=True)

rges = unp.uarray(r, rerr)

theta = rges/2

print('theta')
print(theta)

d = l/(2 * unp.sin(2*np.pi*theta/360))

d1 = l/(2 * unp.sin(2*np.pi*theta[0]/360))

print('d in Angström')
print(d)

verhaeltnis = d1 / d

print('d1/d')
print(verhaeltnis)

a = d * (Ndiamant)**(0.5)

print('a in Angström')
print(a)

theta = unp.nominal_values(theta)
a1 = unp.nominal_values(a)



def f(theta,m , b):
    return m * (np.cos((2*np.pi)/360*theta))**2  + b

params, covariance = curve_fit(f, theta, a1)


errors = np.sqrt(np.diag(covariance))

print('m =', params[0], '\pm', errors[0], 'in Angström')
print('b =', params[1], '\pm', errors[1], 'in Angström')

x = np.linspace(0, 1)

plt.figure(figsize=(10,8))

plt.errorbar(np.cos((2*np.pi)/360*theta)*np.cos((2*np.pi)/360*theta), unp.nominal_values(a), yerr=unp.std_devs(a), fmt='bx', label='Messwerte')
#plt.plot(np.cos((2*np.pi)/360*theta)*np.cos((2*np.pi)/360*theta), a, 'bx', label='Messwerte')
plt.plot(x, params[0]*x+params[1], 'r-', label='Ausgleichsgerade')
plt.legend(loc=2)
plt.savefig('probe1.pdf')



#########################
#                       #
#        Probe 2        #
#                       #
#########################

######################
# Einlesen der Daten #
######################

r, rerr = np.genfromtxt('werte2.txt', unpack=True)

rges = unp.uarray(r, rerr)

theta = rges/2

print('theta')
print(theta)

d = l/(2 * unp.sin(2*np.pi*theta/360))

d1 = l/(2 * unp.sin(2*np.pi*theta[0]/360))

print('d in Angström')
print(d)

verhaeltnis = d1 / d

print('d1/d')
print(verhaeltnis)

a = d * (Ncaesiumchlorid)**(0.5)

print('a in Angström')
print(a)

theta = unp.nominal_values(theta)
a1 = unp.nominal_values(a)



def f(theta,m , b):
    return m * (np.cos((2*np.pi)/360*theta))**2  + b

params, covariance = curve_fit(f, theta, a1)


errors = np.sqrt(np.diag(covariance))

print('m =', params[0], '\pm', errors[0], 'in Angström')
print('b =', params[1], '\pm', errors[1], 'in Angström')

x = np.linspace(0, 1)

plt.figure(figsize=(10,8))

plt.errorbar(np.cos((2*np.pi)/360*theta)*np.cos((2*np.pi)/360*theta), unp.nominal_values(a), yerr=unp.std_devs(a), fmt='bx', label='Messwerte')
#plt.plot(np.cos((2*np.pi)/360*theta)*np.cos((2*np.pi)/360*theta), a, 'bx', label='Messwerte')
plt.plot(x, params[0]*x+params[1], 'r-', label='Ausgleichsgerade')
plt.legend(loc=2)
plt.savefig('probe2.pdf')
