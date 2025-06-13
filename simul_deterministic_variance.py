#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
point here is to check when sb super large what happens
"""

import numpy as np
import math
from scipy.stats import norm
from timeit import timeit
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import animation, rc
from matplotlib.animation import FuncAnimation
from numpy.random import choice
import seaborn as sns

from IPython.display import HTML, Image
rc('animation', html='html5')
gif=False

alllambdas=[]
allsigmas=[]

lambda_step=0.25
    
Is=[]
lambdas=[]
Us=[]

mean=[]
var=[]

prev=1e10

for UU in range(4,5,2):
    Lambda_ode=0.25
    while Lambda_ode < 4.02 and Lambda_ode<=prev:
        print("Lambda_ode="+str(Lambda_ode))
        sbeta=1
        Vs=1/sbeta
        gmax=4
        n=120
        
        U=UU
        fl=0.4
        Lambda=math.sqrt(Lambda_ode)  #SD , math.sqrt(variance)
        beta0=4
        b=3
        d=1
        alpha=1
        
        I0=I0=b/(alpha+d+U * fl) - (d/beta0)
        sigma0=math.sqrt(1) #SD , math.sqrt(variance)  STANDING VARIANCE
        g0=0
        
        n_timesteps=10000
        dt=0.002
        
        
        
        
        def grid_to_geno(x):
            return(-gmax + ((2*gmax) / (n))*x )
        
        def geno_to_grid(x):
            if ((n/2)) + (n/(2*gmax)) * x  != int(((n/2)) + (n/(2*gmax)) * x ):
                print('WARNING, GRID TOMBE PAS JUSTE')
            return(  int(((n/2)) + (n/(2*gmax)) * x ))
        
        def dgeno_to_dist(x):
            return(  round(( n / (2*gmax) ) * abs(x)) )
        
        
        '''grid_to_geno(0)
        grid_to_geno(n)
        geno_to_grid(-gmax)
        geno_to_grid(gmax)
        geno_to_grid(0)
        '''
        def make_normalisation(l):
            #scale of the grid to compute the normalisation constant with SD l
            cst=1 #check that n * cst is a round number
            
            x, y = np.meshgrid(np.linspace(-int(gmax*cst), int(gmax*cst), int(n*cst+1)),
               					np.linspace(int(-gmax*cst), int(gmax*cst), int(n*cst+1)))  
            dstx=np.sqrt(x**2)
                 
            # lower normal part of gaussian
            probas=(1/( l * np.sqrt(2 * np.pi)) * np.exp(-(0.5) * (dstx/ l )**2))
            norma=sum(probas[1])
            return(norma)
        
        
        def make_mat0(gbar):
            normalisation=make_normalisation(Lambda)
            
            
            ### ADD MUT
            
            
            # create mat and S
            M=[0]*(n_timesteps+1)
            M[0]=np.zeros([n+1,n+1])
            M[0][geno_to_grid(gbar),geno_to_grid(gbar)]=I0
            
            S=[0]*(n_timesteps+1)
            S[0]= (alpha + d + U*fl)/beta0
            
            
            
            
            
            # create beta mat
            x, y = np.meshgrid(np.linspace(-int(gmax), int(gmax), int(n+1)),
               					np.linspace(int(-gmax), int(gmax), int(n+1)))  
            x, y = np.meshgrid(np.linspace(-gmax, gmax, int(n+1)),
               					np.linspace(-gmax, gmax, int(n+1)))  
            dst2 = x**2+y**2
            beta_mat=beta0-(dst2 / (2 * Vs))
            
            beta_mat[beta_mat<0] = 0  # CHANGE HERE TO TEST WHAT HAPPENS IF BETA CAN GET NEGATIVE
                            
            #np.tensordot(D, M, 2)
            #mat dx
            
            
            x, y = np.meshgrid(np.linspace(-int(gmax), int(gmax), int(n+1)),
               					np.linspace(int(-gmax), int(gmax), int(n+1)))        
            dx=np.absolute(x[:,:,None,None]-x)
            dy=np.absolute(y[:,:,None,None]-y)
            distx=dx
            disty=dy
            ### Mat D : (i,j,k,l) p(mutation) from (i,j) to (k,l)
            
            Dx=(1/(Lambda * np.sqrt(2 * np.pi)) * np.exp(-(0.5) * (distx/Lambda)**2))
            Dy=(1/(Lambda * np.sqrt(2 * np.pi)) * np.exp(-(0.5) * (disty/Lambda)**2))
            
            D = Dx * Dy / (normalisation**2)
            
            # Standing variance
            
            Dinitx=(1/(sigma0 * np.sqrt(2 * np.pi)) * np.exp(-(0.5) * (distx/sigma0)**2))
            Dinity=(1/(sigma0 * np.sqrt(2 * np.pi)) * np.exp(-(0.5) * (disty/sigma0)**2))    
            Dinit= Dinitx * Dinity / (make_normalisation(sigma0)**2)
            M[0]=np.tensordot(Dinit,M[0],2)
            return(M,S,D,beta_mat,x,y)
            
            
        ### CHECK D ET DINIT
            
            
            
            
        M,S,D,beta_mat,x,y=make_mat0(g0)    
            
            
            
        for t in range(0,n_timesteps):
            M[t+1]=M[t]    +    dt * ( U*(1-fl) * (np.tensordot(D,M[t],2)-M[t])  +   M[t] * (beta_mat * S[t] - alpha - U * fl -d) )
            S[t+1]=S[t]  +  S[t] * dt * (  - d - np.sum(M[t] * beta_mat)) + dt * b
            if round(S[t],4)==round(S[t-5],4)==round(S[t-10],4)==round(S[t-20],4):
                break
        tmax=t+1
        
        I=[np.sum(x) for x in M]
        Is+=[I[tmax]]
        lambdas+=[Lambda_ode]
        Us+=[U]

        toti = np.sum(M[t])
        meansq = np.sum(x*x*M[t]) / toti
        meann = np.sum(x*M[t]) / toti
        mean+=[meann]
        var += [meansq - meann**2]
        
        if I[tmax] < 1e-3:
            print('wha')
            prev=Lambda_ode
        Lambda_ode+=lambda_step
        
        print(4*M[t][:,100].sum()/M[t].sum())
        print('var')
        print(np.log(meansq - meann**2))
      
    
  

cols=[]
epsilon=0.001
for i in Is:
    if i>epsilon:
        cols+=["k"]
    else:
        cols+=["none"]
        
X,Y = np.meshgrid(lambdas,Us)

print(lambdas)


print(np.log(var))
    
  
    
  #-0.37, -0.10, 0.03, 0.12
    