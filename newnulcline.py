import numpy as np
import matplotlib
from matplotlib import pyplot as plt
def f1(x1,x2):
    g10=500.0 #protein production in unbound state, repressive state active
    g11=0.0  #non-zero value will represent protein production in bound state
    V=100.0  #System volume
    h=10.0
    f=1.0
    z1=(h*x2*x2)/(f+h*x2*x2)  #no co-operative/binary binding
    fofx1=-x1+(g11*z1+g10*(1-z1))/V
    return(fofx1)
#---------------------------------------------------------------------------
def f2(x1,x2): #function of second gene
    g10=500.0
    g11=0.0
    V=100.0
    h=10.0
    f=1.0
    z1=(h*x1*x1)/(f+h*x1*x1)
    fofx2=-x2+(g11*z1+g10*(1-z1))/V
    return(fofx2)
#-----------------------------------------------------------------------------
def rootfinder1(x1,x2):
    N=10.0
    max_x2 = max(x2)
    min_x2 = min(x2)
    di = (max_x2-min_x2)/N
    b=0.0
    bracket=np.zeros(0)
    j=0
    x1_rootfound=np.zeros(0)
    for k in x1:
        a=0.0           #initial value of x2
        func1=f1(k,a)          #x1:protein 1 conc,x2: protein 2 conc
        b = 0.0
        while b <= max_x2:
              func2=f1(k,b)
              if func1*func2 < 0:
                 bracket=np.insert(bracket,[0],[b]) # root value of x2(x2*)
                 x1_rootfound=np.insert( x1_rootfound,[0],[k])
                 j=j+1
                 break
              else:
                 b = b + di
#-----------------------------------------------------------------------------
    tol=0.000001
    di=0.0001
    root=[]
    p=0
    a=0.0
    k=x1_rootfound[p]
    func1 = f1(k,a)
    for j in bracket:
        b=j+di
        c = (a+b)/2.0 #bisection
        func = f1(k,b)
        while j <=max(bracket):
              func = f1(k,c)
              if abs(func) <tol:
                 root.append(c)
                 k=x1_rootfound[p]
                 p=p+1
                 break
              else:
                 if func1*func <0:
                    b=c
                    c = (a+b)/2.0
                 else:
                    a=c
                    c = (a+b)/2.0
              if p ==len(bracket):
                 break
    return(x1_rootfound,root)   
#------------------------------------------------------------------------------
def rootfinder2(x1,x2):   #rootfind for second function
    N=10.0
    max_x1 = max(x1)
    min_x1 = min(x1)
    di = (max_x1-min_x1)/N
    b=0.0
    bracket=np.zeros(0)
    j=0
    x2_rootfound=np.zeros(0)
    for k in x2:
        a=0.0           #initial value of x2
        func1=f2(a,k)          #x1:protein 1 conc,x2: protein 2 conce
        b = 0.0
        while b <= max_x1:
              func2=f2(b,k)
              if func1*func2 < 0:
                 bracket=np.insert(bracket,[0],[b]) # root value of x2(x2*)
                 x2_rootfound=np.insert( x2_rootfound,[0],[k])
                 j=j+1
                 break
              else:
                 b = b + di
#------------------------------------------------------------------------------
    tol=0.0000001
    di=0.0001
    root=[]
    p=0
    a=0.0
    k=x2_rootfound[p]
    func1 = f2(a,k)
    for j in bracket:
        b=j+di
        c = (a+b)/2.0 #bisection
        func = f2(b,k)
        while j <=max(bracket):
              func = f2(c,k)
              if abs(func) <tol:
                 root.append(c)
                 k=x2_rootfound[p]
                 p=p+1
                 break
              else:
                 if func1*func <0:
                    b=c
                    c = (a+b)/2.0
                 else:
                    a=c
                    c = (a+b)/2.0
              if p ==len(bracket):
                 break
    return(x2_rootfound,root)
#------------------------------------------------------------------------------
x1=np.arange(0.0,10.0,0.005)
x2=np.arange(0.0,10.0,0.005)
x1_rootfound,x2root=rootfinder1(x1,x2)
x1_up=np.arange(0.,10.0,0.005)
x2_up=np.arange(0.,10.0,0.005) #not using the *root* set x2,instead fresh set
x2_rootfound,x1root=rootfinder2(x1_up,x2_up)

plt.plot(x1_rootfound,x2root,'-',label=r'$dp_1/dt=0$')
plt.plot(x1root,x2_rootfound,'-',label='$dp_2/dt=0$')
plt.legend()
plt.xlabel(r'$p_{1}$')
plt.ylabel(r'$p_{2}$')
plt.show()
plt.savefig('p1p2_ratio20_g10_6000.pdf')
t(elem,y)
print(len(x1root),len(x2root))
