from scipy.integrate import odeint
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from scipy import integrate


OmegaM = 0.315
Omega_RA = 5.042e-5
Omega_DA = 0.68494958
H0 = 67.5
h = 0.67
G = 543.86    # since rho_cr = (3H0**2)/(8*pi*G) = 1
Rho_Cr = 1                    
T0 = 2.725

#Functions of a

# This is hubble rate in terms of conformal time

def H(a):
    return H0* math.sqrt((OmegaM/a) + (Omega_RA/a**2) + (Omega_DA*a**2))


def rho_DM(a):

    return (OmegaM * Rho_Cr /(a**3))

def rho_RA(a):

    return ((Rho_Cr * 4.15 * (10**-5))/((h**2)*(a**4)))


#This function takes in a k mode, then evolves delta through a and return the final value of delta (at a = -1), 

def get_delta_squared(k):
    
   
    
    a_list = np.logspace(-8, -1,num = 1000, endpoint=True)
    

    PHI_0 = 1/(k**3)

    DELTA_0 = 1.5 * PHI_0 

    THETA0_0 = 0.5 * PHI_0 

    THETA1_0 = (-k*PHI_0)/(6*H(10e-8))

    V_0 = (k*PHI_0)/(2*H(10e-8))


    
    def new_model3(state, a):
    
        phi,delta,theta0,theta1,v = state




        dphida = ((4*math.pi*G*a)/(3*(H(a)**2)))*(rho_DM(a)*delta + 4 * rho_RA(a)* theta0)-(phi)/(a)-(3*phi*k**2)/(a*(H(a))**2)

        ddeltada = -3 * dphida - ((1j)*k*v)/(a*H(a))


        dtheta0da = - dphida - (k*theta1)/(a*H(a))


        dtheta1da = (((k*theta0)/3) - ((k*phi)/3))/(a*(H(a)))

        dvda = (-v*H(a) + k*phi)/(a*H(a))



        return [dphida,ddeltada,dtheta0da,dtheta1da, dvda]




    init_state = [PHI_0 ,DELTA_0, THETA0_0, THETA1_0, V_0]


    solver = odeint(new_model3,init_state ,a_list)
    
    
    
    
    
    return ((solver[999][1]-solver[0][1]))			# function return whatever we're looking for 
	





# delsquare = []						#list to store all the delta squares


# k_list = np.logspace(-3,0,num = 1000, endpoint=True)


# for s in k_list:
#     delsquare.append(get_delta_squared(s))		#evaluate delta**2 for each k

    
# print(len(delsquare))
# print(len(k_list))

    
# plt.figure(figsize=(10,10))

# plt.xlabel('log(k)', size = 20)

# plt.ylabel('P(k)',size = 20)

# plt.plot(k_list,delsquare)

# # plt.plot(a_list,b)

# # plt.axvline(x= 10**-3.8, color = 'black', label ='$a_{eq}$')

# plt.xscale('log')
# plt.yscale('log')


# plt.legend()
# plt.show()
    
    
    
    
    
    
    
    
