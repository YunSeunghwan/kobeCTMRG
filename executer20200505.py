from funcs20200505 import *

def init_pfc(m=2):
	global pfc,sizelist
	pfc=Partial_function_calculater
	pfc.blockspin_m=m
def set_xylist():
	global x_list,y_list,dat
	y_list=[]
	x_list=list(np.arange(100)/1000+2.2)
	for size in sizelist:
		y_list.append(list(map(lambda t:CriticalHeat(t,size),x_list)))
		#y_list.append(list(map(lambda t:Z(t,size),x_list)))
	dat=[x_list]
	dat[1:]=y_list
def plot_fsc(plot=1,mark='.-'):
	global fsc
	fsc=FSS(dat)
	fsc.category=sizelist
	fsc.Llist=sizelist
	fsc.xlabel="T/$K_B$J"
	fsc.ylabel="C$_{heat}$"
	fsc.title="C on T of size1000 on each m"
	fsc.show_given_data(save=str("fixedm{}_".format(pfc.blockspin_m)+today),plot=plot,mark=mark)
	# fsc.show_scaled_data(save="fixedscaled"+today)
def get_peak_list():
	global pfc,sizelist,mlist,t_peak_list,c_peak_list
	t_peak_list={}
	c_peak_list={}
	for m in mlist:
		t_peak_list[m]=[]
		c_peak_list[m]=[]
		for size in sizelist:
			def get_C_heat(t):
				return CriticalHeat(t,size)
			print('m:',m,'size:', size)
			init_pfc(m)
			t_peak,c_peak=peak_find(get_C_heat,1,3)
			t_peak_list[m].append(t_peak)
			c_peak_list[m].append(c_peak)

def plot_Tc_on_size(show=True):
	global pfc,sizelist,mlist,t_peak_list,c_peak_list
	for m in mlist:
		plt.plot(sizelist,t_peak_list[m],'.-')
	
	mlegend=mlist[:]
	for m in mlegend:
		mlegend[mlegend.index(m)]=str('m='+str(m))
	
	plt.xlabel('N',fontsize=16)
	plt.ylabel('$T_c$',fontsize=16)
	plt.legend(mlegend)
	plt.title('Tc on each m',fontsize=16)
	if show:
		plt.show()	
def plot_C_peak_on_size(show=True):
	global pfc,sizelist,mlist,t_peak_list,c_peak_list
	for m in mlist:
		plt.plot(sizelist,c_peak_list[m],'.-')

	mlegend=mlist[:]
	for m in mlegend:
		mlegend[mlegend.index(m)]=str('m='+str(m))
	
	plt.xlabel('N',fontsize=16)
	plt.ylabel('$C_{peak}$',fontsize=16)
	plt.legend(mlegend)
	plt.title('$C_{peak}$ on each m',fontsize=16)
	if show:
		plt.show()

def executer_plot_C_heat():
	global pfc,sizelist,x_list,y_list,dat,fsc
	print("start:",perf())
	for m in mlist[:-1]:
		print(m)
		init_pfc(m)
		set_xylist()
		plot_fsc(plot=0)

	init_pfc(m=mlist[-1])
	set_xylist()
	# m_hundred=[x_list,y_list]
	for m in mlist:
		m=str('m='+str(m))
	sizelist=mlist
	print("calculation end:",perf())
	plot_fsc(mark='.-')
def executer_plot_Tc_and_C_heat_on_size():
	global pfc,sizelist,mlist
	get_peak_list()
	print(perf())
	plot_Tc_on_size()
	#plot_C_peak_on_size()

sizelist=[10,20,50,100]
mlist=[8,16,32,64]
executer_plot_Tc_and_C_heat_on_size()

print(c_peak_list)
print(t_peak_list)
#pickle.dump(c_peak_list, open('c_peak_list.pkl', 'wb'))
#pickle.dump(t_peak_list, open('t_peak_list.pkl', 'wb'))





# sizelist=[50,100]
# mlist=[8,16,32,64,128]
# {8: [4.106878983430332, 3.599195255275718], 16: [3.7488484873092602, 4.751018995231728], 32: [3.9912496897579217, 4.729505045565531], 64: [7.376839215286629, 3.918704026488031], 128: [5.010100685029055, 14.396073161769234]}
# {8: [2.2487679439812718, 2.292185144773452], 16: [2.2734193788717336, 2.258733991386231], 32: [2.261959330358195, 2.261959330358195], 64: [2.221891208135877, 2.2219121665537553], 128: [2.2522091386209953, 2.2524308366798467]}