import math
import random
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import copy
from matplotlib.ticker import MaxNLocator

nodemo=True

np.random.seed(1710)
random.seed(1710)

# Creation of a pool of n phenotypes


def pool_of_pheno(I, var0, L, mean0):
    pheno = []
    for i in range(int(I)):
        phen_i = []
        for l in range(L):
            # Corrected line: wrap the value in a list
            phen_i += [np.random.normal(mean0, math.sqrt(var0))]
        pheno.append(phen_i)
    return pheno


def mutation(pheno, U, mut, lam, dt, L, f):
    for i in range(len(pheno) - 1, -1, -1):
        nb_mut = np.random.poisson((1 - f) * U * dt)
        # Initialize newpheno as a copy of the original phenotype
        newpheno = copy.deepcopy(pheno[i])
        while nb_mut > 0 and mut == 1:
            for l in range(L):
                newpheno[l] = np.random.normal(newpheno[l], math.sqrt(lam))
            nb_mut = nb_mut - 1
        pheno[i] = newpheno
    return pheno

# Calculation of r and w for each phenotype


def each_beta(dI, pheno, beta0, Vs, L):
    beta_i = []
    for i in range(int(dI)):
        xo = sum([x**2 for x in pheno[i]])
        #beta_x = beta0*math.exp((-(abs(xo))/(2*Vs)))      #GAUSSIAN FGM
        beta_x = beta0 -xo/(2*Vs)                          #QUADRATIC FGM
        beta_x = max(0.000000001,beta_x)                             #QUADRATIC FGM
        beta_i.append(beta_x)
    beta_pop = np.mean(beta_i)
    return(beta_i, beta_pop)

#Demography and evolution



def growth_evol(dI, pheno, dt, beta_pop, beta_weight, dS,U):
    nb = np.random.poisson((1/n) * dS * dI * beta_pop * dt)
    nd = np.random.poisson(dI*(d+alpha+f*U)*dt)
    if dI+nb-nd < 1 :
        return([],0,0,0,0)
    select_drift = np.random.multinomial(nb, beta_weight)
    pheno_gen2 = []
    for i in range(len(pheno)):
        stack = [pheno[i]]*select_drift[i]
        if stack != []:
            pheno_gen2 += stack
    if nd<len(pheno):
        for i in range(nd):
            pheno.remove(random.choice(pheno))
    else:
        pheno=[]
    pheno += pheno_gen2
    dI = len(pheno)
    return(pheno, dI,nb,nb,nd)


tmax = 10000
dt = 0.1
mut = 1  # Boolean : Mutation or not
U = 0  # Mutation rate
lam = 0.5
# var0 = 10e-10 #Diversity of original population

sbeta=1
Vs = 1/sbeta #keep it at 1
It = []
St = []
t = []
ind = []
beta_t = []

L = 2#n dimansions for UC
# epidemio0

b = 2
d = 1
x = 1


alpha =4
f = 0.4

beta0 = 4

oldU=1e-2 #STARTING U BEFORE TREATMENT



beta0-(sbeta*L*2.8**2)/(2)


if f>0:
    A=(1-f) * lam * b * L*L
    B=8 * f * (beta0 * b - d * (d+alpha))/sbeta
    uc= sbeta*(B+A-math.sqrt((A* (A + 2 * B ))))/(8*d*f*f)
else:
    uc = (4* (beta0 * b - d*(d+alpha))**2) / (b*d*L*L*sbeta*lam)

# Epidemiological parameters
n = [1000]
nb_sim = 10

stress = np.arange(0.6, 3, 0.1)
for i in range(len(stress)):
    stress[i] = round(stress[i], 2)

dmut=0.1
mut_rate = np.arange(-2, 1.5, dmut)
for i in range(len(mut_rate)):
    mut_rate[i] = 10 ** mut_rate[i]
    


for n in n:
    df = pd.DataFrame(index=stress, columns=mut_rate, dtype=float)
    df = df.applymap(lambda x: 0)
    for X in range(len(stress)):
        print(str(X)+"/"+str(len(stress)))
        for Y in range(len(mut_rate)):
            Y=round(Y,4)
            for sim in range(nb_sim):
                #stressX = math.sqrt(stress[X])
                stressX = stress[X]
                U = mut_rate[Y]
                mu = lam*(1-f)*oldU
                gamma = math.sqrt(4*Vs*(alpha+d+U*f)*beta0 + (mu*L**2/4))
                S0 = ((alpha+d+oldU*f)/beta0)+(((L*math.sqrt(mu))*(2*gamma + L*math.sqrt(mu)))/(8*Vs*beta0**2))
                if nodemo:
                    S0=b/d   #nodemo
                I0 = (b/(alpha+d+oldU*f))-(d/beta0)-(d/alpha+d+oldU*f)*(((L*math.sqrt(mu))*(2*gamma + L*math.sqrt(mu)))/(8*Vs*beta0**2))
                dI = I0*n

                dS = S0*n
                var0 = math.sqrt(mu*Vs/S0)
                time = 0
                pheno = pool_of_pheno(dI, var0, L, mean0=stressX)
                oripheno=copy.deepcopy(pheno)
                for ti in range(tmax):
                    if 1 > dI :
                        break
                    time = ti*dt
                    pheno = mutation(pheno, U, mut, lam, dt, L, f)
                    beta = each_beta(dI, pheno, beta0, Vs, L)
                    beta_pop = beta[1]
                    beta_i = beta[0]
                    beta_weight = np.array(beta_i)/sum(beta_i)
                    grw_evo = growth_evol(
                        dI, pheno, dt, beta_pop, beta_weight, dS,U)
                    pheno = grw_evo[0]

                    
                    dI = grw_evo[1]

                    dS = dS + dt * (-dS * ((beta_pop*dI)/n + d) + b*n)    #old
                    if b*n - dS*d>0:
                        bdS = np.random.poisson(dt * (b*n - dS*d))
                    else:
                        bdS = np.random.poisson(-dt * (b*n - dS*d))
                    
                    
                    dS = dS + bdS - grw_evo[2]    #S + birth death - new infected
                    if nodemo:
                        dS=b*n/d  #nodemo


                    
                    It.append(dI)
                    St.append(dS)
                    t.append(time)
                    
                    if ti>20:
                        if It[-1]>It[-5]>It[-10]>It[-20]:
                            break


                if It[-1] < 1:
                    # print("dead")
                    df.loc[stress[X], mut_rate[Y]] += 0.0
                else:
                    # print("survive")
                    #print('len_pheno ='+str(len(pheno)))
                    df.loc[stress[X], mut_rate[Y]] += 1.0


    name = "/Users/martin/Documents/lethal/figures/rescue/nodemorescue_"+str(n)+"f"+str(f)
    #df.to_csv(name+".csv")
    df=df/nb_sim
    sns.heatmap(df.T, vmax=1)
    #plt.yticks(ticks = [1,11,21,31], labels = ["$10^{-2}$", "$10^{-1}$","1","$10^{1}$"]) MAKE THIS AUTO SOMETIMES
    plt.yticks(ticks = [0.5,0.5+(1/dmut),0.5+(2/dmut),0.5+(3/dmut)], labels = ["$10^{-2}$", "$10^{-1}$","1","$10^{1}$"])
    plt.xticks(ticks = [1,11,21], labels = ["1","2","3"],rotation=0)

    plt.plot([0,25],[0.5+((2+math.log10(uc))/dmut),0.5+((2+math.log10(uc))/dmut)], "w--", linewidth=2)
    #ax = plt.figure().gca()
    #ax.set_xticks([1,2,3])
    plt.title("n="+str(n)+"  , f="+str(f))
    plt.xlabel(r"Initial stress $\overline{x}(0)$", fontsize=15)
    plt.ylabel("Mutation rate $U$", fontsize=15)
    plt.gca().invert_yaxis()
    plt.savefig(name+".pdf")
    plt.show()

print("DONE")









if 0:
    n=[100,200,500]
    for nn in n:
        df=pd.read_csv("/Users/martin/Documents/lethal/figures/rescue/rescue_"+str(nn)+"f"+str(f)+".csv", header=0, index_col=0)
        df=df/nb_sim
        sns.heatmap(df.T, vmax=1)
        plt.yticks(ticks = [1,11,21,31], labels = ["$10^{-2}$", "$10^{-1}$","1","$10^{1}$"])

        ax = plt.figure().gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.savefig("/Users/martin/Documents/lethal/figures/rescue/rescue_"+str(nn)+"f"+str(f)+".png")
        plt.show()
        
        

if f>0:
    A=(1-f) * lam * b * L*L
    B=8 * f * (beta0 * b - d * (d+alpha))/sbeta
    uc= sbeta*(B+A-math.sqrt((A* (A + 2 * B ))))/(8*d*f*f)
else:
    uc = (4* (beta0 * b - d*(d+alpha))**2) / (b*d*L*L*sbeta*lam)



f=0

rmax = beta0 *S0 - (d+alpha+U*f)
N0 = I0
betastart = beta0-(sbeta*L* ((2.8**2) + var0)) / (2)
r0 = betastart * S0 - (d+alpha+U*f)
rd=-r0
yd=rd/rmax




#sbeta*(B+sbeta*A-sbeta*math.sqrt((A* (A + 2 * B / sbeta))))/(8*d*f*f)
