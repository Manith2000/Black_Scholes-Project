# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 18:15:08 2020

@author: Manith Adikari
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#Black_Scholes function to output prices for a Call Option 
def Black_Scholes(Company,Stock,Strike_Price,Years_to_Expiry,w):
#---------------User Input Data----------------
    K = Strike_Price        #strike price in $
    E = Years_to_Expiry     #time for option to expire in years
    Xb = Stock              #Ideally, positive infinity
    
    #--------------- Other Data ------------
    #Time array where each time interval is a month
    dt = 0.02
    t = np.arange(0,E+dt,dt) 
    Nt =len(t)
    #Spatial Range
    dx = 0.2
    xb = 2*Xb #Should ideally by positive infinity (details in document)
    xa = -xb #The opposite edge will be equal to xb
    x = np.arange(xa,xb+dx,dx)
    Nx = len(x)
    #Had to limit the xb to a certain value to reduce computing time
    #Ensured that r = 0.5 the max value to achieve accurate values by Crank-Nicolson
    #----------------Crank-Nicholson-------------
    
    #Non-dimensionalise x (non-dimensional x = x/L)
    xnd = x/(xb-xa)
    
    #make output final matrix with values for time and space
    u = np.zeros((Nx,Nt))
    
    #Make Inverse matrix. Identical for all time steps.
    
    A = np.zeros((Nx-2,Nx-2)) 
    #subtract two since we have boundary conditions
    
    r = dt/(dx**2)
    
    for i in range(0,Nx-2):
        A[i,i]= 2+2*r
        if i != 0:
           A[i,i-1]= -r
        elif i != (Nx -2):
            A[i,i+1] = -r
          
    #now invert this matrix
    M = np.linalg.inv(A)
    
    #Call Option Non-dimensionalised boundary and initial conditions 
    U0 = 200 #Arbitrary U value
    u[:,0] = (K/U0)*(np.exp(xnd)-1)
    u[0,:] = 0
    u[Nx-1,:] = 1
    
    for j in range(1,Nt): #for all except initial time
        #Make matrix that will be multiplied by inverse
        B = np.ndarray((Nx-2,1)) #We subtract 1 from Nt due to initial conditions
        B[0,0] = u[2,j-1] + u[0,j-1] + u[0,j] #first B element
        B[Nx-3,0] = u[Nx-1,j-1] + u[Nx-3,j-1] + u[Nx-1,j] #last B element
        B[1:-1,0] = u[3:-1,j-1] + u[1:-3,j-1] #remaining B elements
        u[1:-1,j] = (np.matmul(M,B).ravel())  #Multiply the inverse matrix & B
        #ravel simply turns the 1D (2,) matrix to 2D (2,1) matrix
    
    #------------------------Plots---------------------
    Week = 1 #no. of time intervals equivalent to one week
    #Plot of Option Values for the inputed week
    plt.figure(1, figsize=(8,6))
    plt.plot(x,u[:,w*Week]*U0, label= str(Company))
    plt.xlim(0,Xb*2) #set graph limits (min. and max values)
    plt.ylim(0,u[Nx-10,w*Week]*U0) #Negative values are removed; not realistic
    plt.ylabel("Option Value")
    plt.xlabel("Asset price ($)")
    plt.legend()
    plt.title("Option value vs Asset price at Week " + str(w))
    
    #Plot of Option Values for the following month from inputted week
    plt.figure(2, figsize=(8,6))
    plt.plot(x[:-10],u[:-10,w*Week]*U0, label=str(Company) + ' at Week ' + str(w)) 
    plt.plot(x[:-10],u[:-10,(w+1)*Week]*U0, label=(Company) + ' at Week ' + str(w+1))
    plt.plot(x[:-10],u[:-10,(w+2)*Week]*U0, label=(Company) + ' at Week ' + str(w+2))
    plt.plot(x[:-10],u[:-10,(w+3)*Week]*U0, label=(Company) + ' at Week ' + str(w+3)) 
    #array is sliced for 0:-10 to remove the limit at boundary condition (not accurate at boundaries)
    plt.ylabel("Option Value")
    plt.xlabel("Asset price ($)")
    plt.xlim(0,None)
    plt.ylim(0,u[-10,w*Week]*U0) 
    plt.legend()
    plt.title("Option Value vs Asset Price for following Month")

    #Surface Plot of Option Values by time by asset (i.e. stock) price
    plt.figure(3, figsize=(8,6))
    from mpl_toolkits.mplot3d import Axes3D
    X, T = np.meshgrid(t,x)
    ax = plt.axes(projection = '3d') 
    ax.plot_surface(X,T,u, antialiased=True)
    plt.xlabel("Time")
    plt.ylabel("Asset price")
    ax.set_zlabel("Option Value")
    plt.title("Option Value vs Asset Price Until Option Maturity")

#Apply function on following stocks. 


#------------Before Coronavirus: Stock Price Data from NASDAQ as of 1:55PM 20/02/2020

#Apple Share Price $320.30, Strike Price $100, Expiry Time 2 yrs, 5th Week
Black_Scholes('Apple',320.30,200,2,5)

#Microsoft Share Price $188.70, Strike Price $100, Expiry Time 2 yrs, 5th Week
Black_Scholes('Microsoft',188.70,100,2,5)

#------------After Coronavirus confirmed as Pandemic by WHO: Stock Price Data from NASDAQ as of 5:50PM 03/12/2020
Black_Scholes('Apple after COVID19',258.22,200,2,5)
Black_Scholes('Microsoft after COVID19',146.03,100,2,5)

















