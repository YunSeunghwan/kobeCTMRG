import numpy as np
import math as mt
from matplotlib import pyplot as plt
import matplotlib as mpl
from time import perf_counter as pf

def randomwalk(Randomwalk_electron,potential,direction=4):
	pass
	return Randomwalk_electron
def init_Resistance(resistance_size,V):
	[x,y,z]=resistance_size
	x_init,y_init,z_init=0,*map(int,np.ceil([y/2,z/2]))
	x_goal,y_goal,z_goal=x,*map(int,np.ceil([y/2,z/2]))
	Randomwalk_electron=[x_init,y_init,z_init,x_goal,y_goal,z_goal]
	
	return Randomwalk_electron,potential
def arrived_electron(Randomwalk_electron):
	if Randomwalk_electron[0:3]==Randomwalk_electron[3:6]:return 1
	return 0
def init_electorn(randomwalk_electron):
	pass
	return Randomwalk_electron,potential
def randomwalk_resistance(x,y,z=1,V=1,repeat=100):#electron moves from x=0 to x=L(given x)
	resistance_size=[x,y,z]
	Randomwalk_electron,potential=init_resistance(resistance_size,V)
	list_of_time=[]
	for i in range(repeat):
		time=1
		init_eletron(randomwalk_electron)
		while arrived_electron(Randomwalk_electron):
			Randomwalk_electron=randomwalk(Randomwalk_electron,potential)
			time+=1
		list_of_time.append(time)

def teion4syuume():
	def nyu2(m):
		I=2.5
		beta=1#57.4/180*mt.pi#mt.pi/2
		nyu_L=0.975#MHz
		nyu_Q=20#MHz
		A=24*m*(m-1)-4*I*(I+1)+9
		B=0.25*(6*m*(m-1)-2*I*(I+1)+3)
		res=nyu_L+0.5*nyu_Q*(0.5-m)*(3*mt.cos(beta)**2-1)+nyu_Q**2/nyu_L*(1.5*mt.sin(beta)**2*((A+B)*mt.cos(beta)**2-B))
		return res
	def nyu1(m):
		I=2.5
		beta=1#57.4/180*mt.pi#mt.pi/2
		nyu_L=0.975#MHz
		nyu_Q=20#MHz
		res=nyu_L+0.5*nyu_Q*(0.5-m)*(3*mt.cos(beta)**2-1)
		return res
	def plotdifnu():
		mlist=[2.5,1.5,0.5,-0.5,-1.5]
		x_list1=list(map(nyu1,mlist))
		#x_list2=list(map(nyu2,mlist))
		y_list=[1,2,3,4,5]
		plt.stem(x_list1,y_list,'r-')
		#plt.stem(x_list2,y_list,'b-')
		#plt.xlim([0,max(x_list)*1.1])
		plt.ylim([0,max(y_list)*1.1])
		plt.legend(np.arange(10))
		plt.xlabel('$\\nu$[MHz]')
		plt.show()
	def plotnunu():
		x_list=[20,40]
		y_list=[1,1]
		plt.stem(x_list,y_list,'r-')
		plt.xlim([0,max(x_list)*1.1])
		plt.ylim([0,max(y_list)*1.1])
		plt.xlabel('$\\nu$[MHz]')
		plt.title('I=5/2,NQRspectrum')
		plt.tick_params(axis="y", labelleft=False)
		#plt.savefig('NQRS.eps')
		plt.show()
	plotnunu()
	#plotdifnu()
	#print(*map(nyu,[2.5,1.5,0.5,-0.5,-1.5,-2.5]))

def get_C(T):
	setting=Method_getC
	epsilon=1
	k_B=1
	if setting=='reciprocal':
		w=mt.exp(-epsilon/(k_B*T))  #	=1/W(T)
		C=epsilon**2/(1+w)**2*w/(k_B*T**2)

	elif setting=='normal':
		try:
			W=mt.exp(epsilon/(k_B*T))
			C=epsilon**2/(1+W)**2*W/(k_B*T**2)
		except:
			return 1
	elif setting=='difZ':
		k=1
		beta=1/(k*T)
		def Z(beta):
			J=1
			return 8*np.cosh(3*beta*J)**3+8*np.cosh(beta*J)**3
		def difflnZ(beta):
			e=10**-6
			return (np.log(Z(beta+2*e))-2*np.log(Z(beta+e))+np.log(Z(beta)))/e**2
		C=k*beta**2*difflnZ(beta)
	elif setting=='getZ':
		k=1
		beta=1/(k*T)
		Z=get_Z
		def difflnZ(beta):
			e=10**-4
			return (np.log(Z(beta+e))-2*np.log(Z(beta))+np.log(Z(beta-e)))/e**2
		C=k*beta**2*difflnZ(beta)
	return C

def plot(func,Ts,Te,offset=0.01,save=0,show=1):
	x_list=np.arange(Ts,Te,offset)
	y_list=list(map(func,x_list))
	plt.plot(x_list,y_list,'.-')
	plt.xlabel('T')
	plt.ylabel('$C(T)$')
	plt.title('Critical Heat')
	if save!=0:
		plt.savefig('CH{}.eps'.format(save))
	if show==1:
		plt.show()

def get_Z(beta):
	#boundarycondition_default='open'
	method=Method_getZ
	n=size_of_Isingmodel
	Z=0
	def get_W(beta):
		pass
	def add_W(beta):
		pass
	def init_P(beta):
		n=size_of_Isingmodel
		V_now = get_W(beta)
		for step in range(1,n-1):
			V_now = add_W(V_now,beta)
		return V_now
	def add_P(P_now,beta):
		#	s - s - s   
		#	|       |   
		#	l P_now r     
		#	|       |   
		#	e - e - e  =>P_now(s|l|r|e) = P_now(s|l|r|e) @ P(s|e)
		#    	+	                    = P_now(s|l|r|e)
		#	s - s - s	  				= add_P(P_now)
		#	|   P   |       
		#	e - e - e   
		return P_now @ init_P(beta)
	def get_periodicZ(P_now,beta):
		n=size_of_Isingmodel
		#	Tr_lr(V_now)->Tr_se(V_now)
		#	Tr_...:one more P on that direction, ->trace
		P_now=add_P(P_now,beta).reshape()
		
	def first_V(V_now,beta):
		#	s - s 			     s - s
		#	| W | 	 		     |   |
		#	e - e                l   s - s - s
		#	  +            =     |           |
		#	s'- s'- s - s 	     l   V_now   r
		#	|           |        |			 |
		#	l   V_now   r 		 e - e - e - e
		#	|			|
		#	e - e - e - e    
		#	
		#	: V_now(s|e) = W(s|e) & V_now(s'|se)
		#				 = first_V(V_now)
		pass
	def add_V(V_now,step,beta):
		#	s - s - s' 	s'- s 	s - s - s - s'
		#	| V_now | + | W | = |   V_now   |    (.... = init_P)
		#	e - e -	e'  e'- e   e - e - e - e'
		#	
		#	: V_now(ss'|es') = V_now(ss'|ee') & W(s's|e'e)
		#		 			 = add_W(V_now)
		pass
	if method=='eachstate':
		def sp(x):
			return 2*x-1
		num_array=np.arange(n**2)
		spin_array=np.zeros(n**2)
		for state in range(2**(n**2)):
			
			H=0
			for spin in range(n**2):#spin assignment
				spin_array[spin]=int(state/(2**spin))%2
			for horizontal in range(n):#calculate Hamiltonian for horizontal bond
				for eachbond in range(n-1):
					H-=sp(spin_array[n*horizontal+eachbond])*sp(spin_array[n*horizontal+eachbond+1])
				if boundary_condition=='periodic':
					H-=sp(spin_array[n*horizontal+(n-1)])*sp(spin_array[n*horizontal])
			for vertical in range(n):#calculate Hamiltonian for vertical bond
				for eachbond in range(n-1):
					H-=sp(spin_array[n*eachbond+vertical])*sp(spin_array[n*eachbond+vertical+n])
				if boundary_condition=='periodic':
					H-=sp(spin_array[vertical])*sp(spin_array[vertical+n*(n-1)])
			Z+=np.exp(-beta*H)

	elif method=='eachlayer':
		#	s - s - s   
		#	|       |   
		#	l P_now r     
		#	|       |   
		#	e - e - e  =>P_now(s|l|r|e) = P_now(s|l|r|e) @ P(s|e)
		#    	+	                    = add_P(P_now)
		#	s - s - s	  
		#	|   P   |       
		#	e - e - e   
		P_now = init_P(beta)
		for layer in range(n-1):
			P_now = add_P(P_now,beta)
		if method == 'open' : Z=np.sum(P_now)
		if emthod == 'periodic' : Z=get_periodicZ(P_now,beta)

	elif method=='eachspin':
		#	s - s - s' 	s'- s 	s - s - s - s'
		#	| V_now | + | W | = |   V_now   |    (.... = init_P)
		#	e - e -	e'  e'- e   e - e - e - e'
		#	
		#	: V_now(ss'|es') = V_now(ss'|ee') & W(s's|e'e)
		#		 			 = add_W(V_now)
		#	
		#	
		#	
		#	s - s 			     s - s
		#	| W | 	 		     |   |
		#	e - e                l   s - s - s
		#	  +            =     |           |
		#	s'- s'- s - s 	     l   V_now   r
		#	|           |        |			 |
		#	l   V_now   r 		 e - e - e - e
		#	|			|
		#	e - e - e - e    
		#	
		#	: V_now(s|e) = W(s|e) & V_now(s'|se)
		#				 = first_V(V_now)
		#	
		#	
		#	
		#				
		#	 		  s'- s
		#	s - s' +  | W | 
		#	|   |     e - e   
		#	l   s'- s'- s''     
		#	|           |
		#   l   V_now   r
		#   | 			|
		#	e - e - e - e		
		#	
		#	
		#	: V_now(s|l|r|e) = V_now(ss's''|l|r|e) & W(s's|e)
		#					 = add_V(V_now)
		#	
		#	
		#	
		for floor in range(n-1):
			if floor==1:
				V_now=init_P(n,beta)
			else:
				V_now = first_V(V_now,beta)
				for step in range(1,n-1):
					V_now = add_V(V_now,step,beta)
	return Z

def riron4syuume():
	global Method_getZ,boundary_condition,size_of_Isingmodel,Method_getC
	Method_getZ='eachstate'
	boundary_condition='open'
	size_of_Isingmodel=2
	Method_getC='getZ'
	for mc in ['eachstate']:#['eachstate','eachlayer','eachspin']:
		for size in [2,3,4]:
			if mc=='eachstate' and size>5 : continue
			for bc in ['open']:#,'periodic']:
				size_of_Isingmodel=size
				start_time=pf()
				plot(lambda x: get_C(x)/size**2,1,3,0.1,show=0)
				print('method:',mc,'bc:',bc,'size:',size,pf()-start_time,'sec for calculation')
				plt.legend([2,3,4],fontsize=18)
				# plt.savefig('./figure/CT_{}_{}_{}.eps'.format(mc,bc,size))
				# plt.savefig('./figure/CT_{}_{}_{}.png'.format(mc,bc,size))
	plt.show()
	plt.savefig('./figure/CT_eachstate_open_2to4.pdf')
	plt.savefig('./figure/CT_eachstate_open_2to4.eps')
#riron4syuume()
def bussei5syuume():
	x_list=np.arange(-314,314)/100
	#y_list=list(map(lambda x: abs(np.sin(x)),x_list))
	y_list=list(map(lambda x: 1.5+np.sqrt(2.25-2*np.sin(x)**2),x_list))
	yy_list=list(map(lambda x: 1.5-np.sqrt(2.25-2*np.sin(x)**2),x_list))
	
	plt.plot(x_list,y_list)
	plt.plot(x_list,yy_list)
	plt.show()

bussei5syuume()