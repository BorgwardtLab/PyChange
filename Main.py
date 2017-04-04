import numpy as np 
import Exact as Ex

from matplotlib import pyplot as plt 

#Seq 1: Mean 0,+alpha, -alpha
#Seq 2: Std 0, +alpha, -alpha
#Seq 2: Diff 0, +alpha, -alpha

#TP = 2 Changes indicated as significant
#FP = 2 Changes where there are none
#AUC(alpha)
#RL(alpha) = sum of abs diff to gt TP in units of 1./len(seq)

def mean_shift(seed,alpha):
	r = np.random.RandomState(seed)
	return np.concatenate((r.randn(40), r.randn(20)+alpha,r.randn(20)-alpha))


def var_shift(seed,alpha):
	r = np.random.RandomState(seed)
	return np.concatenate((r.randn(40), r.randn(20)*alpha,r.randn(20)*alpha))

def diff_shift(seed,alpha):
	r = np.random.RandomState(seed)
	return np.cumsum(np.concatenate((r.randn(40), r.randn(20)+alpha,r.randn(20)-alpha)))

methods = ['TT', 'KS', 'MW'] #,'AD']
Names_M = ['Welch', 'Kolmogorov-Smirnov', 'Mann-Whitney'] #,'Anderson-Darling'] 

tests = [mean_shift,var_shift] #,diff_shift]
Names_T = ['Mean', 'Variance'] #, 'Trend']

def Run_all():

	trials = 100
	
	#Get the run lengths: Increase strength of signal and watch what happens
	alpha_range = np.linspace(0,4,10)
	np.save('alpha_range.npy',alpha_range)
	for seq_test,n in zip(tests,Names_T):
		for method,name in zip(methods,Names_M):
			RL =  []
			RL_err = []
			for alpha in alpha_range:
				RL_ = []
				for seed in range(trials):
					seq = seq_test(seed,alpha)
					loc,p = Ex.Exact(seq,method,2,3)
					#print loc
					RL_.append(1.*np.abs(loc[0]-40) + np.abs(loc[1]-60)/80)
				print alpha, np.mean(RL_)
				RL.append(np.mean(RL_))
				RL_err.append(np.std(RL_))
			np.save(name+n+'RL.npy',RL )
			np.save(name+n+'RL_err.npy',RL_err )
	

	#Get the ROC: Play with threshold parameter for signal and no signal and watch what happens
	thres_range=  np.logspace(0,-12,10)
	np.save('thres_range.npy',thres_range)
	for seq_test,n in zip(tests,Names_T):
		for method,name in zip(methods,Names_M):
			FPR = []
			TPR = []
			for thres in thres_range:
				FP = 0
				TP = 0
				for seed in range(trials):
					seq = seq_test(seed,3.)
					loc,p = Ex.Exact(seq,method,2,3)
					if p < thres:
						TP = TP + 1

					seq = seq_test(seed,0.)
					loc,p = Ex.Exact(seq,method,2,3)
					if p < thres:
						FP = FP + 1
				print thres, 1.*TP/trials
				TPR.append(1.*TP/trials)
				FPR.append(1.*FP/trials)
			np.save(name+n+'TPR.npy',TPR )
			np.save(name+n+'FPR.npy',FPR )


def Plot_all():
	alpha_range = np.load('alpha_range.npy')
	for seq_test,n in zip(tests,Names_T):
		f, axarr = plt.subplots(2)
		for method,name in zip(methods,Names_M):
			TPR = np.load(name+n+'TPR.npy')
			FPR = np.load(name+n+'FPR.npy')
			RL = np.load(name+n+'RL.npy')
			RL_err = np.load(name+n+'RL_err.npy')
			print name, n, RL, RL_err
			axarr[0].plot(TPR,FPR,label=name+' AUC: '+ str(np.trapz(TPR,FPR)))
			axarr[1].errorbar(range(len(RL)),RL,yerr=RL_err)#, label=name)
		axarr[0].set_xlabel('FPR',fontsize=16)
		axarr[0].set_ylabel('TPR', fontsize=16)
		axarr[0].set_title('ROC for change in '+n,fontsize=16)
		axarr[0].legend(loc='best',fontsize=14)
		axarr[1].set_xlabel('signal to noise ratio',fontsize=16)
		axarr[1].set_ylabel('Run length', fontsize=16)
		axarr[1].set_title('Run length for change in '+ n,fontsize=16)
		axarr[1].legend(loc='best',fontsize=14)
		plt.savefig(n+'.pdf')
		plt.clf()
		plt.close()

#Run_all()
Plot_all()

"""



np.random.seed(2)
seq = np.concatenate((np.random.randn(40), np.random.randn(20)+1.,np.random.randn(20)-1.))
constraint = [[30,45],[55,65]]
#constraint = []
#print seq
offset=5

locKS,pKS = Ex.Exact(seq,'KS',2,2,0.05,constraint)
print "Kolmogorov-Smirnov:", locKS, pKS

locTT,pTT = Ex.Exact(seq,'TT',2,2,0.05,constraint)
print "Student's T:", locTT, pTT


locMW,pMW = Ex.Exact(seq,'MW',2,2,0.05,constraint)
print "Mann-Whitney:", locMW,pMW

locAD,pAD = Ex.Exact(seq,'AD',2,2,0.05,constraint)
print "Anderson-Darling:", locAD, pAD

"""
