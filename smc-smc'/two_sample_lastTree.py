import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plot
import copy

def twosim(rho):
	height_list=[]
	breakpoint_list=[0]
	current_breakpoint=0
	recomb_time=0
	current_tree_height=np.random.default_rng().exponential(scale=1)
	height_list.append(current_tree_height)
	current_breakpoint+=np.random.default_rng().exponential(scale=1/(rho*current_tree_height))
	while(current_breakpoint<1):
		recomb_time=np.random.default_rng().uniform(0,current_tree_height)
		temp_recoal=np.random.default_rng().exponential(scale=1/2)
		if temp_recoal<current_tree_height-recomb_time:
			if np.random.default_rng().uniform(0,1)<0.5:
				height_list.append(current_tree_height)
			else:
				height_list.append(recomb_time+temp_recoal)
				current_tree_height=recomb_time+temp_recoal
		else:
			current_tree_height+=np.random.default_rng().exponential(scale=1)
			height_list.append(current_tree_height)
		breakpoint_list.append(current_breakpoint)
		current_breakpoint+=np.random.default_rng().exponential(scale=1/(rho*current_tree_height))
	return(breakpoint_list,height_list)

def smcsim(rho):
	height_list=[]
	breakpoint_list=[0]
	current_breakpoint=0
	recomb_time=0
	current_tree_height=np.random.default_rng().exponential(scale=1)
	height_list.append(current_tree_height)
	current_breakpoint+=np.random.default_rng().exponential(scale=1/(rho*current_tree_height))
	while(current_breakpoint<1):
		recomb_time=np.random.default_rng().uniform(0,current_tree_height)
		temp_recoal=np.random.default_rng().exponential(scale=1)
		current_tree_height=recomb_time+temp_recoal
		height_list.append(current_tree_height)
		breakpoint_list.append(current_breakpoint)
		current_breakpoint+=np.random.default_rng().exponential(scale=1/(rho*current_tree_height))
	return(breakpoint_list,height_list)

def firstlast(bplist,heightlist):
	#if len(bplist)>=2:
	return(heightlist[0],heightlist[-1])
	#	return(heightlist[0]/(bplist[1]),heightlist[-1]/(1-bplist[-1]))
		#return(bplist[1],1-bplist[-1])
	#else:
	#	return(False)
def marjwall(num_reps, rho, max_i):
	storer=[]
	for x in range(max_i):
		storer.append([])
#	print(storer,max_i,len(storer))
	output=[]
	for x in range(0,num_reps):
		currentsim=smcsim(rho)
		for indy,height in enumerate(currentsim[1]):
			if indy<len(storer):
				storer[indy].append(height)
				#print(indy,height, storer)
			else: break
	#print("storer: ",storer)
	for y in range(0,max_i):
		if(storer[y]):		
			output.append(np.mean(storer[y]))
		else:
			output.append([])
	return(output)
			
print(marjwall(1000,100,10))

glomfirst=[]
glomlast=[]
##for x in range(0,10000):
##	currentsim=twosim(10)
	#print(currentsim[1])
##	res=firstlast(currentsim[0],currentsim[1])
##	if res:
##		glomfirst.append(res[0])
##		glomlast.append(res[1])
#plot.hist(glomlast,bins=[25,50,75,100,125,150,175,200,225,250,275,300,325])
#plot.show()
#plot.hist(glomfirst,bins=[25,50,75,100,125,150,175,200,225,250,275,300,325])
#plot.show()
##print(np.mean(glomfirst),np.mean(glomlast))
##print(sp.stats.ttest_ind(glomfirst,glomlast))


#heightmeans=[]
#expos=[]
#numrec=[]
#for x in range(0,10000):
#	expos.append(np.random.default_rng().exponential(3))
##	heightmeans.append(np.mean(twosim(10)[1][-1]))
#	numrec.append(len(twosim(10)[0]))
#print(np.mean(numrec))

#print(twosim(10)[1])
