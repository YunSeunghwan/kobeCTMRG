import math as mt
import numpy as np
from matplotlib import pyplot as plt
import pickle
import itertools as it
import sys
from time import perf_counter as perf
from datetime import datetime
#from scipy.sparse.linalg.eigen.arpack import eigsh as scipylargest_eigh
print('math=mt , numpy=np , matplotlib.pyplot=plt , pickle , itertools=it , sys,','\n',
      'perf_counter=perf , clear_output=clear , scipy.sparse.linalg.eigen.arpack.eigsh = scipylargest_eigh','\n',
      'is imported')

class error(Exception):

    def __init__(self,name = 'unknown'):
        self.name = name
        print(self.name,'error occured!!','on',sys._getframe(1).f_code.co_name,'  ,from   ',sys._getframe(2).f_code.co_name,'\n')

def newton_raphson_method(f,y,x_new=0.9999,prec=4): # Calcutlate f^-1(x) for given y...y=f(x)
	x_given = x_new
	while True: # Iteration of newton raphson method
	    x=x_new
	    x_new= x+(y-f(x))/dif(f,x)
	    d=abs(x_new-x)
	#        print('d',d)
	    if d > 10**10:
	        print('x is',x,'x_new is',x_new,'x_given is',x_given)
	        print('')
	        error('divergence occured.')
	    if d <= 10**-prec: # Return x_new for certain degree of precision
	        return x_new
    
def newton_raphson_method_dif_given(f,fd,y,x_new=0.9999,prec=4): # Calcutlate f^-1(x) for given y...y=f(x)
    x_given = x_new
    while True: # Iteration of newton raphson method
        x=x_new
        x_new= x+(y-f(x))/fd(x)
        d= abs(x_new - x)
#        print('d',d)
        if d > 10**10:
            print('x is',x,'x_new is',x_new,'x_given is',x_given)
            print('')
            error('divergence occured.')
        if d <= 10**-prec: # Return x_new for certain degree of precision
            return x_new

def dif_tanh(x): # Function mt.exp gives result accurato to 11 places
    return 4/(mt.exp(x)+mt.exp(-x))**2

def tanh(x): # Function mt.exp gives result accurato to 11 places
    return (mt.exp(x)-mt.exp(-x))/(mt.exp(x)+mt.exp(-x))

def cosh(x): # Function mt.exp gives result accurato to 11 places
    return (mt.exp(x)+mt.exp(-x))/2

def sinh(x): # Function mt.exp gives result accurato to 11 places
    return (mt.exp(x)-mt.exp(-x))/2

def dif(f,x,precision=8): # Return difference of f on x with certain degree ofprecision
    e = 10**-precision
    return (f(x+e)-f(x))/e

def f_with_dif(f,x,order = 1,precision=6): # Return difference of f on x with certain degree ofprecision
    e = 10**-precision
    fx= f(x)
    fxe = f(x + e)
    if order == 1:
        return np.array([fx, (fxe - fx) / e])
    if order == 2:
        fxee = f(x+2*e)
        return np.array([fx, (fxe - fx) / e, (fxee - 2*fxe + fx) / e**2])

def power_iteration(A , num , method = 'default' ): # A as unit matrix
    if method == 'division' and num>=3:
        def Adot_sbs(num,W,V):
            V=V.reshape(2,-1)#1
            U=t_init(W,V)#2
            for k in range(1,num-2):#3
                U=t_next(num,k,W,U)
            U=U_to_U_dash(U)#4
            W=W.reshape(4,4)
            V_ret=t_last(W,U)#5
            return V_ret
            
        def eigenvalue_sbs(num,W,v):
            Av = Adot_sbs(num,W,v)
            return v.dot(Av)

        def V_0(n):#0
            return np.zeros(2**n)/np.sqrt(2**n)

        def V_dash(v):#1
            return v.reshape(2,-1)

        def t_init(W_given,V_given):#2
            res=W_given.dot(V_given)
            tra=res.reshape(4,2,2,-1).transpose(1,2,0,3)
            trans=np.array([tra[0][0],tra[1][1]]).reshape(2,-1)
            return trans

        def t_next(num,k,W_given,U_given):#3
            res=W_given.dot(U_given).reshape(2,2,2,2**k,2,2,-1).transpose(2,5,0,4,3,1,6)
            tra=np.array([[res[0][0][0][0],res[0][0][1][1]],[res[1][1][0][0],res[1][1][1][1]]])
            tran=tra.transpose(0,2,1,3,4).reshape(2,-1)
            return tran

        def U_to_U_dash(U_given):#4
            return U_given.reshape(2,-1,2).transpose(0,2,1).reshape(4,-1)

        def t_last(W_given,U_given):#5
            res=W_given.dot(U_given)
            tra=res.reshape(2,2,-1,2).transpose(0,3,2,1)
            return np.array([tra[0][0],tra[1][1]]).transpose(1,0,2).reshape(-1)
      
        W = A.transpose(2,3,1,0).reshape(8,2)
        v = np.ones(2**num) / np.sqrt(2**num)
        ev = eigenvalue_sbs(num,W,v)

        while True:
            Av = Adot_sbs(num,W,v)
            v_new = Av / np.linalg.norm(Av)
            ev_new = eigenvalue_sbs(num,W,v_new)
            if np.abs(ev - ev_new) < 10**-2:
                break

            v = v_new
            ev = ev_new

        return ev_new, v_new

    if method == 'default' or num<3:
        n,d = A.shape

        v = np.ones(d) / np.sqrt(d)
        ev = eigenvalue(A, v)

        while True:
            Av = A.dot(v)
            v_new = Av / np.linalg.norm(Av)

            ev_new = eigenvalue(A, v_new)
            if np.abs(ev - ev_new) < 10**-pp:
                break

            v = v_new
            ev = ev_new

        return ev_new, v_new

def bisection_bethod(f,x_0,x_1):
    pass

def least_squares_method(x_list,y_list):#y=ax+b
    if len(x_list)!=len(y_list):
        error('number doesnt match!')
    x=np.sum(np.array(x_list))
    xx=np.sum(np.array(x_list)**2)
    y=np.sum(np.array(y_list))
    xy=np.sum(np.array(x_list)*np.array(y_list))
    n=len(x_list)
    a=(n*xy-x*y)/(n*xx-x*x)
    b=(xx*y-x*xy)/(n*xx-x*x)
    return a,b

def peak_find(f,x_min,x_max,point=4,f_min=None,f_max=None,precision=6,printx=0):
    x_itv=x_max-x_min
    if printx==1:
        print('x_min is ',x_min)
    if x_itv<10**-precision:
        x_output = (x_min+x_max)/2
        y_output = (f_min+f_max)/2
        return x_output,y_output

    x_list=np.ones(point)*x_min+np.arange(point)*x_itv/(point-1)
    f_list=np.zeros(point)
    if f_min!=None:
        f_list[0] = f_min
    else:
        f_list[0] = f(x_list[0])
    if f_max!=None:
        f_list[point-1] = f_max
    else:
        f_list[point-1] = f(x_list[point-1])

    for i in range(1,point-1):
        f_list[i] = f(x_list[i])
    f_peak = max(f_list)
    x_peak = x_list[list(f_list).index(f_peak)]
    x_peak_index = list(x_list).index(x_peak)
    if x_peak_index ==0:
        l_idx = 0
        r_idx = 1
    elif x_peak_index ==point-1:
        r_idx = point-1
        l_idx = point-2
    else:
        l_idx= x_peak_index-1
        r_idx= x_peak_index+1
        

    return peak_find(f,x_min=x_list[l_idx],x_max=x_list[r_idx],f_min=f_list[l_idx],f_max=f_list[r_idx],precision=precision,printx=printx)


if __name__=='__main__':
    x=np.array([2.07944154, 2.30258509 ,2.48490665, 2.63905733 ,2.77258872])
    y=np.array([1.2663554,  1.63720125, 1.62377828 ,1.80379848 ,1.81902006])
    print(least_squares_method(x,y))