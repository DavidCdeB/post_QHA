import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sympy as sym
from matplotlib import cm
from sympy.solvers import solve
#from sympy import Symbol, diff, sin

# Function to fit:
def func(X, a0, a1, a2, a3, a4, a5):
     x, y = X
     return a0 + a1*y + a2*x + a3*y**2 + a4*x**2  + a5*x*y 

# Load data:
y_data, z_data, x_data  = np.loadtxt('./1__scatter_xyz_sort_y_wise.dat').T
y_data_2, z_data_2, x_data_2  = np.loadtxt('./2__scatter_xyz_sort_y_wise.dat').T

# Calling non linear curve_fit
popt, pcov = curve_fit(func, (x_data, y_data), z_data) 
popt_2, pcov_2 = curve_fit(func, (x_data_2, y_data_2), z_data_2) 

print 'popt = ', popt
print 'pcov = ', pcov

perr = np.sqrt(np.diag(pcov))

print 'pcov = ', pcov
print 'perr = np.sqrt(np.diag(pcov))', perr


#print 'a0 = ', a0
a0 =popt[0] 
a1 =popt[1] 
a2 =popt[2] 
a3 =popt[3] 
a4 =popt[4]    
a5 =popt[5]  

print 'a0 = ', a0

a0_s2 =popt_2[0] 
a1_s2 =popt_2[1] 
a2_s2 =popt_2[2] 
a3_s2 =popt_2[3] 
a4_s2 =popt_2[4]    
a5_s2 =popt_2[5]  

print 'a0_s2 = ', a0_s2

print """ 

The equations are the following:

G_I  (T, P) = a0 + a1*y + a2*x + a3*y**2 + a4*x**2  + a5*x*y
G_II (T, P) = a0_s2 + a1_s2*y + a2_s2*x + a3_s2*y**2 + a4_s2*x**2  + a5_s2*x*y

"""
print('G_I  (T, P) = ({a0}) + ({a1})*P + ({a2})*T  ({a3})*P**2  ({a4})*T**2  + ({a5})*T*P'.format(a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, a5 = a5))

print """
"""
print('G_II  (T, P) = ({a0_s2}) + ({a1_s2})*P + ({a2_s2})*T  ({a3_s2})*P**2  ({a4_s2})*T**2  + ({a5_s2})*T*P'.format(a0_s2 = a0_s2, a1_s2 = a1_s2, a2_s2 = a2_s2, a3_s2 = a3_s2, a4_s2 = a4_s2, a5_s2 = a5_s2))

print """
"""

print """
G_I  (T, P) = G_II (T, P)

"""

x_mesh = np.linspace(10.0000000000000, 2000.0000000000000, 20)
x_mesh_2 = np.linspace(10.0000000000000, 2000.0000000000000, 20)
print x_mesh[0]
print x_mesh[-1]
y_mesh = np.linspace(-4.4121040129800, 10.8555489379000, 20)
y_mesh_2 = np.linspace(8.0622039627300, 17.6458151433000, 20)
print y_mesh[0]
print y_mesh[-1]

xx, yy = np.meshgrid(x_mesh, y_mesh)
xx_2, yy_2 = np.meshgrid(x_mesh_2, y_mesh_2)

z_fit = a0 + a1*yy + a2*xx + a3*yy**2 + a4*xx**2  + a5*xx*yy		
z_fit_2 = a0_s2 + a1_s2*yy_2 + a2_s2*xx_2 + a3_s2*yy_2**2 + a4_s2*xx_2**2  + a5_s2*xx_2*yy_2		


## Solving the intersection:
print """ 

The equations are the following:

z_I  (x, y) = a0 + a1*y + a2*x + a3*y**2 + a4*x**2  + a5*x*y
z_II (x, y) = a0_s2 + a1_s2*y + a2_s2*x + a3_s2*y**2 + a4_s2*x**2  + a5_s2*x*y

"""
print('z_I  (x, y) = ({a0}) + ({a1})*y + ({a2})*x  ({a3})*y**2  ({a4})*x**2  + ({a5})*x*y'.format(a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, a5 = a5))

print """
"""
print('z_II  (x, y) = ({a0_s2}) + ({a1_s2})*y + ({a2_s2})*x  ({a3_s2})*y**2  ({a4_s2})*x**2  + ({a5_s2})*x*y'.format(a0_s2 = a0_s2, a1_s2 = a1_s2, a2_s2 = a2_s2, a3_s2 = a3_s2, a4_s2 = a4_s2, a5_s2 = a5_s2))

print """
"""

print """
The intersection of both surfaces is satisfied when:

z_I(x, y) = z_II(x, y)

In other words, I am looking for the expression of the function y=y(x)

"""
# Setting "x" and "y" to be symbolic:
x, y = sym.symbols('x y', real=True)

def z_I(x,y):
        return   a0 + a1*y + a2*x + a3*y**2 + a4*x**2  + a5*x*y

def z_II(x,y):
        return   a0_s2 + a1_s2*y + a2_s2*x + a3_s2*y**2 + a4_s2*x**2  + a5_s2*x*y

sol = sym.solve(z_I(x,y) - z_II(x,y), y)
print 'sol =', sol


# For obtaining the plot of the two branches y=y(x), we need np.sqrt
def y_sol_1(x):
    return 0.000319359080035813*x - 1.22230952828787e-15*np.sqrt(-1.07919313606384e+24*x**2 + 2.00910207755286e+28*x - 1.12101975048632e+30) + 10.6162640815323


def y_sol_2(x):
   return 0.000319359080035813*x + 1.22230952828787e-15*np.sqrt(-1.07919313606384e+24*x**2 + 2.00910207755286e+28*x - 1.12101975048632e+30) + 10.6162640815323



print ' y_sol_1(10.0) = ', y_sol_1(10.0+0j)
print ' y_sol_2(10.0) = ', y_sol_2(10.0+0j)


cross = sym.solve(sol[0]-sol[1])
print ' cross = ', cross


##### Plotting:
# Use this to turn on matplotlib 1.5 defaults:
#style.use('classic')

# set "fig" and "ax" varaibles
fig = plt.figure()
ax = fig.gca(projection='3d')

# If you set "classic", use also this, 
# in order to capture z-labels right:
#z_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
#ax.zaxis.set_major_formatter(z_formatter)

# Plot the original function
ax.plot_surface(xx, yy, z_fit, color='y', alpha=0.5)
ax.plot_surface(xx_2, yy_2, z_fit_2, color='g', alpha=0.5)

# Plot the initial scattered points
ax.scatter(x_data, y_data, z_data, color='r', marker='o') # 'ro') #color='r', marker='o')
ax.scatter(x_data_2, y_data_2, z_data_2, '^s') #color='r', marker='o')


ax.set_xlabel('T (K)')
ax.set_ylabel('P (GPa)')
ax.set_zlabel('\nGibbs free energy / F.unit (a.u.)', linespacing=3)

# New figure for the y=y(x) function:
fig = plt.figure()
x = np.linspace(10.0, 2000.0, 10000)
plt.plot(x, y_sol_1(x))
plt.plot(x, y_sol_2(x))
plt.xlabel('T (K)')
plt.ylabel('P (GPa)')
plt.title('Exact expression of P=P(T)\nas a result of making $G^{I}(T,P)=G^{II}(T,P)$')
tics_shown =  [10, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250]
plt.xticks(tics_shown)
plt.grid()

# New figure for the y=y(x) function in circle:
fig = plt.figure()
x_circle = np.linspace(10.0, 2000.0*100, 10000*100)
plt.plot(x_circle, y_sol_1(x_circle))
plt.plot(x_circle, y_sol_2(x_circle))
plt.xlabel('T (K)')
plt.ylabel('P (GPa)')
plt.title('Exact expression of P=P(T)\nas a result of making $G^{I}(T,P)=G^{II}(T,P)$')
#tics_shown =  [10, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250]
#plt.xticks(tics_shown)
plt.grid()



# New figure for the G1(T,P)-G2(T,P) vs P at a fixed T:
#fig = plt.figure()
#plt.plot(press,diff)
#plt.xlabel('P (GPa)')
#plt.ylabel('$G^{I}(T,P) - G^{II}(T,P)$ / F.unit  (a.u.)')
#plt.title('At T={fixed_T}K'.format(fixed_T=fixed_T))
#plt.grid()
#

#===============
#  1st subplot: G1(T,P)
#===============
# set up the axes for the first plot
fig = plt.figure()
ax = fig.add_subplot(2, 3, 1, projection='3d', title = '$G^{I}(T,P)$')
ax.plot_surface(xx, yy, z_fit, cmap = cm.viridis, antialiased=False, alpha = 0.5)


#===============
# 2nd subplot: G2(T,P)
#===============
# set up the axes for the second plot
ax = fig.add_subplot(2, 3, 2, projection='3d', title = '$G^{II}(T,P)$')

#ax = fig.gca(projection='3d')
ax.plot_surface(xx, yy, z_fit_2, cmap = cm.viridis, antialiased=False, alpha = 0.5)


#===============
# 3rd subplot: G1(T,P)-G2(T,P)
#===============
# set up the axes for the second plot
#ax = fig.add_subplot(2, 3, 3, projection='3d', title='G(T,P)-G(T,P)')

# set "fig" and "ax" varaibles
ax = fig.add_subplot(2, 3, 3, projection='3d', title='$G^{I}(T,P) - G^{II}(T,P)$')

ax.plot_surface(xx, yy, z_fit - z_fit_2, rstride=1, cstride=1, cmap = cm.viridis, antialiased=False, alpha = 0.5)

#ax.set_xlabel('T (K)')
#ax.set_ylabel('P (GPa)')
#ax.set_zlabel('\n$G^{I}(T,P) - G^{II}(T,P)$  / F.unit (a.u.)', linespacing=3)

#===============
# 4th subplot: contour
#===============
# set up the axes for the second plot
#ax = fig.add_subplot(2, 3, 4, projection='3d', title="contour: levels = [0]")
ax = fig.add_subplot(2, 3, 4, projection='3d', title="contour")


ax.contour(xx, yy, z_fit - z_fit_2, levels = [0])
ax.set_zlim(zmin = -2)

#===============
# 5th subplot: projection of the contour on the xy plane
#===============
# set up the axes for the second plot
#ax = fig.add_subplot(2, 3, 5, projection='3d', title="              projection of the contour on the xy plane: zdir='z', offset=-2, levels = [0]")

ax = fig.add_subplot(2, 3, 5, projection='3d', title="projection of the contour on the xy plane")


ax.contour(xx, yy, z_fit - z_fit_2, zdir='z', offset=-2, levels = [0])
ax.set_zlim(zmin = -2)


#plt.legend()
plt.show()

