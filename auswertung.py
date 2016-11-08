import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit
import scipy.odr.odrpack as odrpack

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True


plt.rcParams.update({'font.size': 18})



# Konstanten

l = 1.54286

Ndiamant = np.array([3, 4, 8, 11, 16, 19, 24, 27, 32, 35,36])

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
#r = 180-r
#r=r[::-1]

rges = unp.uarray(r, rerr)

theta = rges/2

print('theta')
print(theta)

d = l/(2 * unp.sin(2*np.pi*theta/360))

print('d in Angström')
print(d)

verhaeltnis = d[0] / d

print('d1/d')
print(verhaeltnis)

################plot
a = d * (Ndiamant)**(0.5)

print('a in Angström')
print(a)

theta1 = unp.nominal_values(theta)
theta2 = unp.std_devs(theta)
a1 = unp.nominal_values(a)
a2 = unp.std_devs(a)

#print(theta1, theta2)


def f(theta1, m, b):
    return m * (np.cos((2*np.pi)/360*theta1))**2  + b

params, covariance = curve_fit(f, theta1, a1)

errors = np.sqrt(np.diag(covariance))

print('m =', params[0], '\pm', errors[0], 'in Angström')
print('b =', params[1], '\pm', errors[1], 'in Angström')

################################################################################
##
##    Test orth. regression
##
################################################################################

def f(B, theta1):
    return B[0] * (np.cos((2*np.pi)/360*theta1))**2  + B[1]

x = theta1
y = a1
sx = theta2
sy = a2

linear = odrpack.Model(f)

myData = odrpack.RealData(x, y, sx=sx, sy=sy)
myOdr = odrpack.ODR(myData, linear , beta0=[1,1])
myOdr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
out = myOdr.run()
#out.pprint()

coeff = out.beta[::-1]
err   = out.sd_beta[::-1]

print(coeff[0], err[0])
print(coeff[1], err[1])

################################################################################



#print(np.cos(np.deg2rad(theta2))**2)

x = np.linspace(0, 1, 10000)

plt.figure(figsize=(10,8))

plt.errorbar(np.cos((2*np.pi)/360*theta1)*np.cos((2*np.pi)/360*theta1), unp.nominal_values(a),xerr=2*np.cos(np.deg2rad(theta1))*np.sin(np.deg2rad(theta1))*np.deg2rad(theta2), yerr=unp.std_devs(a), fmt='bx', label='Messwerte')
#plt.plot(np.cos((2*np.pi)/360*theta1)**2, coeff[1]*x+coeff[0], 'bx', label='Messwerte')
plt.plot(x, coeff[1]*x+coeff[0], 'r-', label='Ausgleichsgerade')

plt.xlabel(r"$\mathrm{cos^2(\theta)}$", fontsize=20)
plt.ylabel(r"$\mathrm{a\,\, /\,\, 10^{-10}\,m}$", fontsize=20)
plt.legend(loc=2)
plt.savefig('probe1.pdf', bbox_inches='tight', pad_inches=0.5)



#########################
#                       #
#        Probe 2        #
#                       #
#########################

######################
# Einlesen der Daten #
######################

r, rerr = np.genfromtxt('werte22.txt', unpack=True)

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

theta1 = unp.nominal_values(theta)
theta2 = unp.std_devs(theta)
a1 = unp.nominal_values(a)
a2 = unp.std_devs(a)
print(theta1)
print((np.cos((2*np.pi)/360*theta1))**2)
#print(theta1, theta2)


def f(theta1, m, b):
    return m * (np.cos((2*np.pi)/360*theta1))**2  + b


params, covariance = curve_fit(f, theta1, a1)


errors = np.sqrt(np.diag(covariance))

print('m =', params[0], '\pm', errors[0], 'in Angström')
print('b =', params[1], '\pm', errors[1], 'in Angström')


################################################################################
##
##    Test orth. regression
##
################################################################################

def f(B, theta1):
    return B[0] * (np.cos((2*np.pi)/360*theta1))**2  + B[1]

x = theta1
y = a1
sx = theta2
sy = a2

linear = odrpack.Model(f)

myData = odrpack.RealData(x, y, sx=sx, sy=sy)
myOdr = odrpack.ODR(myData, linear , beta0=[1,1])
myOdr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
out = myOdr.run()
#out.pprint()

coeff = out.beta[::-1]
err   = out.sd_beta[::-1]

print(coeff[0], err[0])
print(coeff[1], err[1])

################################################################################



x = np.linspace(0, 1)

plt.figure(figsize=(10,8))

plt.errorbar(np.cos((2*np.pi)/360*theta1)**2, unp.nominal_values(a), xerr=2*np.cos(np.deg2rad(theta1))*np.sin(np.deg2rad(theta1))*np.deg2rad(theta2), yerr=unp.std_devs(a), fmt='bx', label='Messwerte')
#plt.plot(np.cos((2*np.pi)/360*theta)*np.cos((2*np.pi)/360*theta), a, 'bx', label='Messwerte')
plt.plot(x, coeff[1]*x+coeff[0], 'r-', label='Ausgleichsgerade')
plt.xlabel(r"$\mathrm{cos^2(\theta)}$", fontsize=20)
plt.ylabel(r"$\mathrm{a\,\, /\,\, 10^{-10}\,m}$", fontsize=20)
plt.legend(loc=2)
plt.savefig('probe2.pdf', bbox_inches='tight', pad_inches=0.5)
