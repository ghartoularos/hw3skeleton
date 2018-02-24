from hw3skeleton import smithwaterman, pairsio
import numpy as np
import pandas as pd
from tqdm import tqdm
from itertools import product
from scipy.spatial.distance import squareform
import os
import sys
from copy import copy

negpairs = pairsio.get_pairs('pairs/Negpairs.txt')
pospairs = pairsio.get_pairs('pairs/Pospairs.txt')

def p1q1(opentop=2,exttop=2,do_all=False,pkling=False):
	data = pd.DataFrame(np.zeros((20,4)),columns=[
		'Gap Open', 'Gap Extension',
		'Thresh TP > 70%','FP'])
	if do_all == True:
		opentop = 21
		exttop = 6
	else:
		assert type(opentop) == int and opentop > 0 and \
		type(exttop) == int and exttop > 0, \
		'Range indices accept integers only.'

		assert opentop > 1 and exttop > 1, \
		'Range upper bound must be greater than 1.'
	row = -1
	for i in range(1,opentop):
		for j in range(1,exttop):
			posscores = []
			print('####### Gapopen = %d ########' % i)
			print('Pospairs:')
			for pair in tqdm(pospairs):
				score = smithwaterman.main('submats/BLOSUM50',
					seq1=pair[0], seq2=pair[1],
					gapopen=-i, gapext=-j, js=True)
				posscores.append(score)
			threshold = np.percentile(posscores,30.0)
			# np.percentile uses a linear interpolation when the quantile is
			# between two data points
			negscores = []
			print('Negpairs:')
			for pair in tqdm(negpairs):
				score = smithwaterman.main('submats/BLOSUM50',
					seq1=pair[0], seq2=pair[1],
					gapopen=-i, gapext=-j,
					js=True)
				negscores.append(score)
			FP = float(sum([score > threshold for score in negscores]))/\
				len(negscores)
			row += 1
			data.iloc[row,0] = i
			data.iloc[row,1] = j
			data.iloc[row,2] = threshold
			data.iloc[row,3] = FP
			if pkling == True:
				data.to_pickle('q1data_ext%d.pkl' % j)
	return data

def p1q2a(nummats=1,do_all=False,pkling=False):

	files = ['submats/' + f for f in os.listdir('submats/')]
	data = pd.DataFrame(np.zeros((len(files),3)),columns=[
		'Matrix', 'Thresh TP > 70%','FP'])

	row = -1
	files = ['submats/' + f for f in os.listdir('submats/')]

	if do_all == True:
		nummats = len(files)

	assert type(nummats) == int and nummats > 0, \
	'Invalid number of matrices.'

	assert nummats <= len(files), \
	'More matrices requested than available.'

	files = [files[i] for i in range(nummats)]
	for f in files:
		posscores = []
		print('####### Matrix = %s ########' % f.split('/')[1])
		print('Pospairs:')
		for pair in tqdm(pospairs):
			score = smithwaterman.main(f,
				seq1=pair[0], seq2=pair[1],
				gapopen=-5, gapext=-3, js=True)
			posscores.append(score)
		threshold = np.percentile(posscores,30.0)
		# np.percentile uses a linear interpolation when the quantile is
		# between two data points
		negscores = []
		print('Negpairs:')
		for pair in tqdm(negpairs):
			score = smithwaterman.main(f,
				seq1=pair[0], seq2=pair[1],
				gapopen=-5, gapext=-3,
				js=True)
			negscores.append(score)
		FP = float(sum([score > threshold for score in negscores]))/\
			len(negscores)
		row += 1
		data.iloc[row,0] = f.split('/')[1]
		data.iloc[row,1] = threshold
		data.iloc[row,2] = FP
		if pkling == True:
			data.to_pickle('q2adata.pkl')
	return data

def p1q2b(nummats=1,do_all=False,pkling=False):

	files = ['submats/' + f for f in os.listdir('submats/')]

	ROCdatax = pd.DataFrame(np.zeros((len(files),50)),index=
		[f.split('/')[1] for f in files])

	if do_all == True:
		nummats = len(files)

	assert type(nummats) == int and nummats > 0, \
	'Invalid number of matrices.'

	assert nummats <= len(files), \
	'More matrices requested than available.'

	files = [files[i] for i in range(nummats)]

	for i in range(len(files)):
		posscores = []
		print('####### Matrix = %s ########' % files[i].split('/')[1])
		print('Pospairs:')
		for pair in tqdm(pospairs):
			score = smithwaterman.main(files[i],
				seq1=pair[0], seq2=pair[1],
				gapopen=-5, gapext=-3, js=True)
			posscores.append(score)
		negscores = []
		print('Negpairs:')
		for pair in tqdm(negpairs):
			score = smithwaterman.main(files[i],
				seq1=pair[0], seq2=pair[1],
				gapopen=-5, gapext=-3,
				js=True)
			negscores.append(score)
		for j in range(len(ROCdatax.columns)):
			'''
			Because np.percentile uses a linear interpolation when the 
			quantile is between two data points, it will simply floor
			to the lowest value of posscores at percentile 0, but 
			lowest score in possscores still might be higher than some
			scores in negscores, so ROC curve would not hit (1,1);
			therefore, at last iteration, set threshold to 0.
			'''
			if j != len(ROCdatax.columns) - 1:
				threshold = np.percentile(posscores,100-2*(j+1))
			else:
				threshold = 0
			FP = float(sum([score > threshold for score in negscores]))/\
						len(negscores)
			ROCdatax.iloc[i,j] = FP
		if pkling == True:
			ROCdatax.to_pickle('q2bdata.pkl') # pickling at every iteration, in case it doesn't finish
	return ROCdatax

def p1q3(pkling=False):

	ROCdatax = pd.DataFrame(np.zeros((2,50)),index=
	['PAM100-Raw','PAM100-Normalized'])

	posscores = []
	pospairlens = []
	print('Pospairs:')
	for pair in tqdm(pospairs):
		score = smithwaterman.main('submats/PAM100',
			seq1=pair[0], seq2=pair[1],
			gapopen=-5, gapext=-3, js=True)
		pospairlens.append(min(len(pair[0]),len(pair[1])))
		posscores.append(score)
	negscores = []
	negpairlens = []
	print('Negpairs:')
	for pair in tqdm(negpairs):
		score = smithwaterman.main('submats/PAM100',
			seq1=pair[0], seq2=pair[1],
			gapopen=-5, gapext=-3,
			js=True)
		negpairlens.append(min(len(pair[0]),len(pair[1])))
		negscores.append(score)
	for i in range(len(ROCdatax)):
		if i == 1:
			posscores = [k/m for k, m in product(posscores,pospairlens)]
			negscores = [k/m for k, m in product(negscores,negpairlens)]
		for j in range(len(ROCdatax.columns)):
			if j != len(ROCdatax.columns) - 1:
				threshold = np.percentile(posscores,100-2*(j+1))
			else:
				threshold = 0
			FP = float(sum([score > threshold for score in negscores]))/\
						len(negscores)
			ROCdatax.iloc[i,j] = FP
	if pkling == True:
		ROCdatax.to_pickle('q3data.pkl') # pickling at every iteration, in case it doesn't finish
	return ROCdatax

def p2q1(randompairs=10, q3=False):
	if q3 == True:
		pop = ['submats/' + f for f in os.listdir('submats/')]
		matioind = pop.index('submats/MATIO')
	else:
		pop = ['submats/' + f for f in os.listdir('submats/') if f != 'MATIO']
	subs = list()
	cols = list()
	for i in range(len(pop)):
		submat = smithwaterman.get_submat(pop[i])
		if i == 0:
			subs.append(submat)
			cols.append(list(submat.columns))
		else:
			if list(submat.columns) != cols[i-1]:
				try:
					submat = submat[cols]
				except:
					sys.exit('Matrices not compatible')
			subs.append(submat)
			cols.append(list(submat.columns))
	subs = np.stack(subs)
	cols = cols[0]
	gens = 5
	allfits = list()
	for generation in range(gens):
		print('###########Generation: %d############' % generation)
		np.save('genq3%d.pkl' % generation,subs)
		randints = np.random.choice(len(pospairs),randompairs,replace=False)
		fits = list()
		for i in range(len(subs)):
		# for i in tqdm(range(len(subs))):
			TPs = list()
			print('Sub%d - Getting Fitness:' % i)
			for FP in (0,0.1,0.2,0.3):
				submat = pd.DataFrame(subs[i],index=cols,columns=cols)
				negscores = []
				for pair in tqdm([negpairs[rando] for rando in randints]):
					score = smithwaterman.main(submat,
						seq1=pair[0], seq2=pair[1],
						gapopen=-5, gapext=-3,
						js=True,fm=True)
					negscores.append(score)
				threshold = np.percentile(negscores,((1 - FP)*100))
				# np.percentile uses a linear interpolation when the quantile is
				# between two data points
				posscores = []
				for pair in tqdm([pospairs[rando] for rando in randints]):
					score = smithwaterman.main(submat,
						seq1=pair[0], seq2=pair[1],
						gapopen=-5, gapext=-3, js=True, fm=True)
					posscores.append(score)
				TP = float(sum([score > threshold for score in posscores]))/\
					len(posscores)
				print(TP)
				TPs.append(TP)
			fits.append(sum(TPs))
		allfits.append(copy(fits))
		parents = list()
		if generation == 0 and q3==True:
			fits[matioind] = 4 # to make sure its derived from MATIO matrix
		for _ in range(2):
			ind = np.argmax(fits)
			fits[ind] = -1
			parents.append(subs[ind])
		subs = list()
		conds = list()
		diags = list()
		for parent in parents:
			diag = [parent[i][i] for i in range(len(parent))]
			diags.append(diag)
			conds.append(squareform(parent,checks=False))
		for i in range(len(pop)):
			sizecond = len(conds[0]) # size of condensed matrix
			sizediag = len(diags[0])
			randcond1 = set(np.random.choice(sizecond,sizecond/2,replace=False))
			randcond2 = set(range(sizecond)) - randcond1
			randdiag1 = set(np.random.choice(sizediag,sizediag/2,replace=False))
			randdiag2 = set(range(sizediag)) - randdiag1
			randconds = [list(randcond1), list(randcond2)]
			randdiags = [list(randdiag1), list(randdiag2)]
			newconds = np.full(sizecond,np.inf)
			newdiags = np.full(sizediag, np.inf)
			for n in range(2):
				for j in randconds[n]:
					newconds[j] = conds[n][j]
				for j in randdiags[n]:
					newdiags[j] = diags[n][j]
			sub = squareform(newconds)
			for i in range(len(sub)):
				sub[i][i] = newdiags[i]
			subs.append(sub)
		subs = np.stack(subs)
	np.save('fitnessesq3.pkl',allfits)
	return pd.DataFrame(parents[0],index=cols,columns=cols)

def p2q2n3(pkling=True):
	mats = ['MATIO','MATOPT']
	ROCdatax = pd.DataFrame(np.zeros((2,50)),index=mats)
	matopt = pd.read_pickle('MATOPT.pkl')
	for matnumber in range(len(mats)):
		print('Pospairs:')
		posscores = []
		pospairlens = []
		if mats[matnumber] == 'MATIO':
			for pair in tqdm(pospairs):
					score = smithwaterman.main('submats/MATIO',
						seq1=pair[0], seq2=pair[1],
						gapopen=-5, gapext=-3, js=True)
					posscores.append(score)
		else:
			for pair in tqdm(pospairs):
					score = smithwaterman.main(matopt,
						seq1=pair[0], seq2=pair[1],
						gapopen=-5, gapext=-3, js=True, fm=True)
					posscores.append(score)	
		negscores = []
		negpairlens = []
		print('Negpairs:')
		if mats[matnumber] == 'MATIO':
			for pair in tqdm(negpairs):
					score = smithwaterman.main('submats/MATIO',
						seq1=pair[0], seq2=pair[1],
						gapopen=-5, gapext=-3, js=True)
					negscores.append(score)
		else:
			for pair in tqdm(negpairs):
					score = smithwaterman.main(matopt,
						seq1=pair[0], seq2=pair[1],
						gapopen=-5, gapext=-3, js=True, fm=True)
					negscores.append(score)	
		for j in range(len(ROCdatax.columns)):
			if j != len(ROCdatax.columns) - 1:
				threshold = np.percentile(posscores,100-2*(j+1))
			else:
				threshold = 0
			FP = float(sum([score > threshold for score in negscores]))/\
						len(negscores)
			ROCdatax.iloc[matnumber,j] = FP
	if pkling == True:
		ROCdatax.to_pickle('p2q3data.pkl') # pickling at every iteration, in case it doesn't finish
	return ROCdatax