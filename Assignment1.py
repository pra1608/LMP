# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:16:53 2020

@author: Ravi Raj
"""

import numpy as np
import matplotlib.pyplot as plt
import math



# Defining material properties
# Chosen Material is copper
k= 385 # (W/mK)
rho= 8960 #(kg/m**3)
C=376.812 #(J/Kg-K)
alpha= k/(rho*C)





# Choosing Laser fluence and pulse duration
Las_Ful = 10**4 # ( J/m^2)
T_P= 10*10**(-9) # 10 ns
f= 10**6 # 1 Mega-hertz


d= 100*10**(-6) # beam diameter (m)
area= (np.pi)*((d**2)/4) # beam area (m^2)
E= Las_Ful*area # Pulse Energy
P= E* f # Laser Power for unit cross sectional area
print(f"laser power is equal to {P} W")
I_avg= P/area
I_max= I_avg/(f*T_P)







# Fresnel Equation for copper and air interface at 1064 nm wavelength
# Ref (https://refractiveindex.info/?shelf=3d&book=metals&page=copper)
lmda= 1064*10**(-9)  # nm
n_c= 0.37863
k_c= 7.0660
n_a= 1.00027398
k_a= 0
# Reflectivity Calculation
R= (((n_c- n_a)**2+(k_c- k_a)**2)/((n_c+ n_a)**2+ (k_c+ k_a)**2))
# Absorption Coefficient
ab= (4*(np.pi)* k_c)/ lmda
print(f"Absorption Coefficient is equal to {ab}")
print(f"Reflectivity is equal to {R}")





I= (1-R)*I_max
p_d = np.sqrt(alpha*T_P)
print(f"Penetration depth is equal to {p_d}")

q_abs = I/p_d # internal volumetric heat generation due to absorption of laser







# Defining number of nodes and steps for space and time
# space node calculation
m=30  # number of elements/ number of nodes -1 in space
a=0 # top node co-ordinate where laser is irradiated
b= 5*10**(-6)  # bottom node co-ordinate (5 um)
delx= (b-a)/(m)  # space element size
x= np.linspace(a,b,m+1)




# time step calculation
# time steps
t0=0  # initial time
tf=T_P# final time (s)
delt= (((1/3)*((delx)**2))/(alpha))  # time step
n=int((tf-t0)/delt)  # number of time steps






# Heat source calculation
Q = np.zeros(m+1)
for i in range(0,m+1):
        Q[i] =   q_abs * np.exp((-1)*(ab*(x[i])))
       

 


#Using equation 6 Heat fluid flow lecture note. 
# Defining A & B constants
A= (alpha*delt)/(delx*delx)
b= ((alpha*delt)/(k))




# defining zero matrices to contain temperature data 
D= np.zeros([n+1,m+1])
K= np.zeros([m+1,m+1])
L= np.zeros([n+1,m+1])
M= np.zeros([m+1,m+1])
N= np.zeros([n+1,m+1])







# Defining Starting Condition (supposing that the material is at a particular temperature)
TS=10  # in degree celcius
T=np.zeros(m+1)
for i in range (m+1):
    T[i]=TS    












# Explicit Case    
# Filling A matrix with grid fourier numbers in explicit case
K[0,0]=1
K[m,m]=1
for i in range (1,m):
    j=i-1
    if j<m-1:
        K[i,j]= A
        K[i,j+1]= 1- 2*A
        K[i,j+2]= A
 
# filling a matrix containg temperature by FDM method explicit condition Using Matrix
Tn=T
Tn1=Tn
for j in range (n+1):
    
    L[j,:]= Tn1
    Tn=Tn1
    Tn1= np.dot(K,Tn) + b*Q
    Tn1[m]=Tn1[m-1]














# Implicit Case    
# Filling A matrix with grid fourier numbers in implicit case
M[0,0]=1
M[m,m]=1
for i in range (1,m):
    j=i-1
    if j<m-1:
        M[i,j]= (-1)* A
        M[i,j+1]= 1 + 2*A
        M[i,j+2]= (-1)* A
 
# filling a matrix containg temperature by FDM method implicit condition Using Matrix
Tn=T
Tn1=Tn
for j in range (n+1):
    
    N[j,:]= Tn1
    Tn=Tn1
    Tn1= np.dot((np.linalg.inv(M)),Tn) + b* Q
    Tn1[m]=Tn1[m-1]











# Closed form solution    
def ierfc (u):
    return (np.exp(-(u**2))/(np.sqrt(np.pi)))- u*(1- math.erf(u))

# filling a matrix containg temperature data by closed form solution

for j in range (1,n+1):
    for i in range (m+1):
        D[j,i]= TS + (2*I/k)* ((np.sqrt(alpha*j*delt))*ierfc((i*delx)/(2*(np.sqrt(alpha*j*delt)))))









#plotting figures
        
        
        
# plot closed form case
fig2 = plt.figure(figsize=(6,5))
left, bottom, width, height= 0.05,0.05,0.8,0.8 
ax2 = fig2.add_axes([left, bottom, width, height]) 

    # Plotting figure with labels, title and legends
plt.plot(x, D[int(n/6),:], "red",x, D[int(2*n/6),:], "orange",x, D[int(3*n/6),:], "yellow")
plt.plot(x, D[int(4*n/6),:], "green",x, D[int(5*n/6),:], "blue",x, D[n,:], "black")
plt.legend([(t0+delt*(n/6)),(t0+delt*(2*n/6)), (t0+delt*(3*n/6)), (t0+delt*(4*n/6)), (t0+delt*(5*n/6)), (t0+delt*(6*n/6)) ])
ax2.set_title('Plot of "T" vs "x" with increasing time for copper rod \n with exponetial heat source  \n using closed form solution')
ax2.set_xlabel('x-axis')
ax2.set_ylabel('T-axis')
plt.show()










# plot explicit case
fig4 = plt.figure(figsize=(6,5))
left, bottom, width, height= 0.05,0.05,0.8,0.8 
ax4 = fig4.add_axes([left, bottom, width, height]) 

    # Plotting figure with labels, title and legends
plt.plot(x, L[int(n/6),:], "red",x, L[int(2*n/6),:], "orange",x, L[int(3*n/6),:], "yellow")
plt.plot(x, L[int(4*n/6),:], "green",x, L[int(5*n/6),:], "blue",x, L[n,:], "black")
plt.legend([(t0+delt*(n/6)),(t0+delt*(2*n/6)), (t0+delt*(3*n/6)), (t0+delt*(4*n/6)), (t0+delt*(5*n/6)), (t0+delt*(6*n/6)) ])
ax4.set_title('Plot of "T" vs "x" with increasing time for copper rod \n with exponetial heat source within the pulse time \n using FDM Explicit Matrix method')
ax4.set_xlabel('x-axis')
ax4.set_ylabel('T-axis')
plt.show()









# plot implicit case
fig6 = plt.figure(figsize=(6,5))
left, bottom, width, height= 0.05,0.05,0.8,0.8 
ax6 = fig6.add_axes([left, bottom, width, height]) 

    # Plotting figure with labels, title and legends
plt.plot(x, N[int(n/6),:], "red",x, N[int(2*n/6),:], "orange",x, N[int(3*n/6),:], "yellow")
plt.plot(x, N[int(4*n/6),:], "green",x, N[int(5*n/6),:], "blue",x, N[n,:], "black")
plt.legend([(t0+delt*(n/6)),(t0+delt*(2*n/6)), (t0+delt*(3*n/6)), (t0+delt*(4*n/6)), (t0+delt*(5*n/6)), (t0+delt*(6*n/6)) ])
ax6.set_title('Plot of "T" vs "x" with increasing time for copper rod \n with exponetial heat source within the pulse time \n using FDM Implicit Matrix method')
ax6.set_xlabel('x-axis')
ax6.set_ylabel('T-axis')
plt.show()