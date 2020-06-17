from NM import *
LSM=least_squares_method

now=datetime.now()
today=str(now.year)[2:]+str(now.month).zfill(2)+str(now.day).zfill(2)+'_'+str(now.hour).zfill(2)+str(now.minute).zfill(2)

class error(Exception):
	def __init__(self,name = 'unknown'):
		self.name = name
		print(self.name,'error occured!!','on',sys._getframe(1).f_code.co_name,',from',sys._getframe(2).f_code.co_name,'\n')

class CTM():
	def __init__(self,J=1,beta=1,C_scale=1):
		self.J=J
		self.beta=beta
		self.temp=1/beta
		self.C_scale=C_scale
		self.Pnow=0
		self.Cnow=0
		self.boundary_condition="free"
		self.blockspin_m=0
		self.lnC_max=0
		self.lnP_max=0
		self.projC=np.kron(np.ones(4).reshape(4,1),np.eye(4))
		self.projP=np.kron(np.ones(2).reshape(2,1),np.eye(4))
		self.W=0
#			   c
# variables:  d b   ...but same as 4types
#			   a
	def getW(self,beta=0):
		if beta==0:
			beta=self.beta
		abcd=np.array([
			[4,0,0,0],
			[0,-4,0,0],
			[0,0,-4,0],
			[0,0,0,4]],dtype='float64')
		return np.exp(beta*self.J*abcd)

	def setC_init(self):#make C of size W
		self.W=self.getW()
		Cnow=self.W
		if self.boundary_condition=="fixed":
			Cnow[1]=[0,0,0,0]
			Cnow[2]=[0,0,0,0]
			Cnow[3]=[0,0,0,0]
			Cnow.transpose(1,0)
		C=(Cnow.reshape(-1)@self.projC).reshape(2,2)
		#W(ab|cd)->C(a|b) or C(c|d)
		self.C_now=C#C(a|b)
	def setC_added(self):
		C=self.C_now
		size=np.shape(C)[0]
		Pnow=self.P_now
		CP=C@Pnow.reshape(size,-1)
		CPP=CP.reshape(size,2,size).transpose(1,2,0).reshape(-1,size)@Pnow.reshape(size,-1)
		CPPW=CPP.reshape(2,size,2,size).transpose(3,1,0,2).reshape(-1,4)@self.W
		C_added=CPPW.reshape(size,size,2,2).transpose(0,2,1,3).reshape(2*size,2*size)
		self.C_now=C_added#C(a|b)
	def setP_init(self):#make P size of W
		Pnow=self.W
		self.W28=self.W.reshape(2,8)
		if self.boundary_condition=="fixed":
			Pnow[1]=[0,0,0,0]
			Pnow[3]=[0,0,0,0]
		P=(Pnow.reshape(2,-1)@self.projP).reshape(2,2,2)#W(ab|cd)->P(a|c|d)=P(b|d|a)
		self.P_now=P#P(b|d|a)
	def setP_added(self):
		P=self.P_now
		size=np.shape(P)[0]
		PW=P.transpose(0,2,1).reshape(-1,2)@self.W28
		P_added=PW.reshape(size,size,2,2,2).transpose(0,2,3,1,4).reshape(2*size,2,2*size)
		self.P_now=P_added#P(b|d|a)

	def blockspinC(self):
		C = self.C_now
		m = self.blockspin_m
		O = np.array(np.linalg.eig(C)[1])#,dtype='float64')
		OT = O.transpose()
		AT = OT[:m]
		A = AT.transpose()
		self.blockspinA = A
		C_appr = AT@C@A
		C_max=C_appr.max()
		self.C_now = C_appr/C_max
		self.lnC_max += np.log(C_max)

	def blockspinP(self):
		P = self.P_now
		size = np.shape(P)[0]
		m = self.blockspin_m
		A = self.blockspinA
		AT = A.transpose()
		P_appr_ = AT@P.reshape(size,2*size)
		P_appr__ = P_appr_.reshape(m*2,size)@A
		P_appr = P_appr__.reshape(m,2,m)
		P_max = P_appr.max()
		self.P_now = P_appr/P_max
		self.lnP_max += np.log(P_max)

	def setC_sizeof(self,nC):
		self.setC_init()
		#print('1',self.C_now)
		self.setP_init()
		#print('2',self.P_now)
		self.lnP_max=0
		self.lnC_max=0
		for n in range(1,nC):#order important!
			self.setC_added()
			#print('3',self.C_now)
			self.setP_added()
			#print('4',self.P_now)
			if self.blockspin_m!=0 and self.blockspin_m<np.shape(self.C_now)[0]:
				self.blockspinC()
				#print('5',self.C_now)
				self.lnC_max += 2*self.lnP_max
				self.blockspinP()
				#print('6',self.P_now)
		#print('7',self.C_now)
		#C,P_now becomes(m,m)&(m,2,m)

	def sumC(self,nC=0):
		self.setC_sizeof(nC)
		C=np.array(self.C_now)
		ro=np.trace(C@C@C@C)
		return ro

class data():
	def __init__(self,data):
		self.data=data
		self.numofdata=len(data)
		self.xlabel="T/$K_B$J"
		self.ylabel="C"
		self.xfontsize="16"
		self.yfontsize="16"
		self.tfontsize="16"
		self.category=list(np.arange(self.numofdata))
		self.title="title_undefined"
	def show_given_data(self,save=0,plot=1,mark='.-'):
		for i in range(1,self.numofdata):
			plt.plot(self.data[0],self.data[i],mark)
		plt.xlabel(self.xlabel,fontsize=self.xfontsize)
		plt.ylabel(self.ylabel,fontsize=self.yfontsize)
		plt.legend(self.category)
		plt.title(self.title,fontsize=self.tfontsize)
		if save!=0:
			try:
				plt.savefig('./figure/{}_{}.png'.format(self.title,save))
			except:
				print("save must be type of str->save with str(save)")
				plt.savefig('./figure/{}_{}.png'.format(self.title,str(save)))
		if plot==1:
			plt.show()

class FSS(data):
	def __init__(self,data,where="peak",optfunc="quad"):
		self.data=data
		self.numofdata=len(data)
		self.xlabel="T/$K_B$J"
		self.ylabel="C"
		self.xfontsize="16"
		self.yfontsize="16"
		self.tfontsize="16"
		self.category=list(np.arange(self.numofdata))
		self.given_data=data
		self.alpha=1
		self.aL=1
		self.num_pickup=5
		self.numofL=self.numofdata-1
		self.scaled=0
		self.where=where
		self.optfunc=optfunc
		self.title="title_undefined"

	def get_a(self,L):
		i=self.Llist.index(L)
		b,c= self.get_bc(self.Llist[i])
		x,y=self.processed_data[0][i],self.processed_data[1][i]
		yi=c*np.array(y)-c
		xi=np.array(x)**2
		k=-np.sum(xi*yi)/np.sum(xi**2)
		#print(k)
		return np.sqrt(k)

	def get_bc(self,L):
		Lidx=self.category.index(L)
		ylist=self.given_data[Lidx+1]
		c=max(ylist)
		b=self.given_data[0][ylist.index(c)]
		#print(L,b)
		return b,c
	def set_picked_data(self):
		self.picked_data=np.zeros(4*(self.numofL)*self.num_pickup).reshape(2,self.numofL,2*self.num_pickup)
		for i in range(self.numofL):
			b,c= self.get_bc(self.Llist[i])
			if self.where=="peak":
				Tc_idx=self.given_data[0].index(b)
			sidx,lidx=Tc_idx-self.num_pickup,Tc_idx+self.num_pickup
			self.picked_data[0][i]=self.given_data[0][sidx:lidx]
			self.picked_data[1][i]=self.given_data[i+1][sidx:lidx]
	def set_processed_data(self):
		self.processed_data=np.array(self.picked_data)
		for i in range(self.numofL):
			b,c= self.get_bc(self.Llist[i])
			self.processed_data[0][i]=self.picked_data[0][i]-b
			self.processed_data[1][i]=self.picked_data[1][i]/c
	def set_list_abc(self):
		self.set_picked_data()
		#print(self.picked_data[0])
		self.set_processed_data()
		#print(self.picked_data[0])
		self.list_a=np.zeros(self.numofL)
		self.list_b=np.zeros(self.numofL)
		self.list_c=np.zeros(self.numofL)
		#print(self.picked_data[0])
		for i in range(self.numofL):
			self.list_a[i]=self.get_a(self.Llist[i])
			self.list_b[i],self.list_c[i]=self.get_bc(self.Llist[i])

	def get_coefficient_and_exponent(self,ylist):
		lnL_list=np.log(np.array(self.Llist))
		#print('lnx:',lnL_list)
		lny_list=np.log(np.array(ylist))
		#print('lny:',lny_list)
		lnex,lnco=LSM(lnL_list,lny_list)
		return np.exp(lnco),lnex
		
	def get_coa_and_alpha(self):
		co,ex=self.get_coefficient_and_exponent(self.list_a)
		return co,ex
	def get_cob_and_beta(self):
		co,ex=self.get_coefficient_and_exponent(self.list_b)
		return co,ex
	def get_coc_and_gamma(self):
		co,ex=self.get_coefficient_and_exponent(self.list_c)
		return co,ex
	def set_scaled_data(self,pltrange):
		self.scaled_data=[[],[]]
		coa,alpha=self.get_coa_and_alpha()
		cob,beta=self.get_cob_and_beta()
		coc,gamma=self.get_coc_and_gamma()
		if pltrange == "all":
			self.scaled_data[0]=np.array([self.given_data[0]]*self.numofL)
			self.scaled_data[1]=np.array(self.given_data[1:])
		if pltrange == "peak":
			self.scaled_data[0]=np.array(self.picked_data[0])
			self.scaled_data[1]=np.array(self.picked_data[1])
		print("a=",coa,"*L^",alpha)
		print("b=",cob,"*L^",beta)
		print("c=",coc,"*L^",gamma)
		for i in range(self.numofL):
			L=self.Llist[i]
			self.scaled_data[0][i]-=cob*L**beta
			self.scaled_data[1][i]/=coc*L**gamma
			self.scaled_data[0][i]*=coa*L**alpha/np.sqrt(coc*L**gamma)
	def set_FSS(self,plotrange="all"):
		self.set_list_abc()
		self.set_scaled_data(pltrange=plotrange)

	def show_scaled_data(self,save=0):
		for i in range(self.numofL):
			plt.plot(self.scaled_data[0][i],self.scaled_data[1][i],'.-')
		plt.xlabel(self.xlabel,fontsize=self.xfontsize)
		plt.ylabel(self.ylabel,fontsize=self.yfontsize)
		plt.legend(self.category)
		if self.title!=0:
			plt.title(self.title)
		else:
			plt.title(str("scaled"+self.ylabel[0]+"on"+self.xlabel[0]),fontsize=self.tfontsize)
		if save!=0:
			try:
				plt.savefig('./figure/scaled_{}_{}.png'.format(self.title,save))
			except:
				print("save must be type of char->save with str(save)")
				plt.savefig('./figure/scaled_{}_{}.png'.format(self.title,str(save)))
		plt.show()

def Z(temp,Csize):
	Partial_function_calculater.beta=1/temp
	Z=Partial_function_calculater.sumC(nC=round(Csize))
	return Z

def CriticalHeat(temp,Csize):#L=2*Csize
	beta=1/temp
	def Zbeta(beta):
		Partial_function_calculater.beta=beta
		Z=Partial_function_calculater.sumC(nC=round(Csize))
		return Z
	def lnZ(beta):
		return np.log(Zbeta(beta))+4*Partial_function_calculater.lnC_max
	def diff(f,x):
		e=10**-2
		fxee=f(x+2*e)
		fxe=f(x+e)
		fx=f(x)
		Fxe=(fxee-fxe)/e
		Fx=(fxe-fx)/e
		return (Fxe-Fx)/e
	return diff(lnZ,beta)*beta**2/(2*Csize)**2


Partial_function_calculater=CTM()