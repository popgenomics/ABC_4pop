#!/gepv/home2/croux/bin/pypy

# #!/usr/bin/pypy
## #!/usr/bin/python


# ##!/usr/bin/pypy

# assumes four populations with the following species tree:
# ((A,B), (C,D))
import os
import sys


def cr_sqrt(x):
	# returns the square root of a list of variables [x]
	if x == 0.0:
		return 0.0
	else:
		M = 1.0
		xN = x 
		while xN >= 2.0:
			xN = 0.25*xN
			M = 2.0*M
		while xN < 0.5:
			xN = 4.0*xN
			M = 0.5*M
		A = xN
		B = 1.0-xN
		while 1==1:
			A = A*(1.0+0.5*B)
			B = 0.25*(3.0+B)*B*B
			if B < 1.0E-15:
				return A*M


def cr_sum(x):
	# returns the sum of a list by removing 'na' values (-9)
	x2 = [ i for i in x if i!=-9 ]
	if len(x2) != 0:
		return(sum(x2))
	else:
		return(-9)


def cr_mean(x):
	# returns the mean of a list
	# removes the 'na' values from the list (-9)
	x2 = [ i for i in x if i!=-9 ] # -9 = 'na' value
	nElement = len(x2)
	if nElement == 0:
		return(-9)
	else:
		return(sum(x2)/(1.0 * nElement))

def cr_std(x, exp_X):
	# returns the standard variation of a list
	# removes the 'na' values from the list (-9)
	x2 = [ i for i in x if i!=-9 ] # -9 = 'na' value
	nElement = len(x2)
	if nElement == 0:
		return(-9)
	else:
		if sum(x2) == 0:
			return(0)
		else:
			A = sum( [ i**2 for i in x2 ] )
			A = A/(1.0 * nElement)
			return(cr_sqrt(A-exp_X**2))


def cr_pearsonR(x, y):
	# computes correlation between arrays x and y
	## removes 'na' values (-9)
	x2 = [ x[i] for i in range(len(x)) if x[i]!=-9 and y[i]!=-9 ]
	y2 = [ y[i] for i in range(len(y)) if x[i]!=-9 and y[i]!=-9 ]
	
	sumXi = 0.0
	sumYi = 0.0
	sumSquareX = 0.0
	sumSquareY = 0.0
	sumXiYi = 0.0
	
	nX = len(x2)

	if nX!=0:
		for i in range(nX):
			sumXi += x2[i];
			sumYi += y2[i];
			sumSquareX += x2[i]*x2[i];
			sumSquareY += y2[i]*y2[i];
			sumXiYi += x2[i]*y2[i];
		
		numerator = nX*sumXiYi - sumXi * sumYi
		denom1 = cr_sqrt(nX*sumSquareX - sumXi*sumXi)
		denom2 = cr_sqrt(nX*sumSquareY - sumYi*sumYi)
		if denom1 == 0 or denom2 == 0:
			return(0)
		else:
			return(numerator/(denom1*denom2))
	else:
		return(-9)


def compFreq(sequences, segsites):
	# returns derived allele frequency for a number of 'segsites' positions
	nDerAll = []
	nInd = len(sequences)
	nPair = nInd*(nInd-1)/2.0
	for i in range(segsites):
		nDerAll.append(0)
		for j in sequences:
			if j[i] == "1":
				nDerAll[i] += 1.0
	pi = [ i * (nInd-i) / nPair for i in nDerAll ]
	freq = [ i/(1.0 * nInd) for i in nDerAll ]
	res = {}
	res['nDer'] = nDerAll # list of the number of occurence of the derived allele at each position
	res['pi_SNPs'] = pi # list of pi measured at each position
	res['pi'] = sum(pi) # total pi = sum of the pi = total number of differences over all positions + over pairwise comparisons between species
	res['freq'] = freq # frequencies of the derived allele at each position
	res['nSNPs'] = len([ i for i in freq if i>0.0 and i<1.0 ]) # number of SNPs (with frequencies in ]0, 1[) in the alignement
	return(res)


def piTot(nDerA, nDerB, nSamA, nSamB, segsites):
	# returns pi tot for the pooled populations from two vectors of allele count per locus
	piT = []
	nTot = nSamA + nSamB
	for i in range(segsites):
		if nDerA==nSamA and nDerB==nSamB:
			piT.append(-9)
		if nDerA==0 and nDerB==0:
			piT.append(-9)
		else:
			tmp = nDerA[i] + nDerB[i]
			tmp = (nTot - tmp) * tmp / (nTot*(nTot-1)/2.0) # "n ancesral" x "n derived" / C(2,k)
			piT.append(tmp)
	return(piT)


def Fst(piA, piB, piT):
	# returns Fst as: 1 - mean_pi_over_populations / pi_in_the_whole_alignment
	if piT==0:
		res = -9
	else:
		res = 1.0 - (piA+piB)/(2.0*piT)
	return(res)


def sites(freqA, freqB, segsites):
	# test whether the SNP is a polymorphism specific to species A (sxA), species B (sxB), is found in both species (Ss) or differentialy fixed among species (Sf)
	sxA, sxB, ss, sf, same = 0, 0, 0, 0, 0
	successive_ss = [0]
	previous_site = ""
	for i in range(segsites):
		if freqA[i] == 0:
			if freqB[i] == 0:
				same += 1
			if freqB[i] == 1:
				sf += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sf"
			if freqB[i] < 1:
				sxB += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxB"
			continue
		if freqA[i] == 1:
			if freqB[i] == 1:
				same += 1
			if freqB[i] == 0:
				sf += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sf"
			if freqB[i] > 0:
				sxB += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxB"
			continue
		else:
			if freqB[i] == 0 or freqB[i] == 1:
				sxA += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxA"
			else:
				ss += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss[len(successive_ss)-1] += 1
				previous_site = "ss"
			continue
	res = {'sxA':sxA, 'sxB':sxB, 'sf':sf, 'ss':ss, 'same':same, 'successive_ss':max(successive_ss)}
	return(res)


def tajD(pi, theta, n, S, a1, a2):
	# returns Tajima D
	# pi = pi for a locus
	# theta = theta for a locus
	# n = sample size
	# S = number of polymorphic positions
	# a1 = sum(1/i)
	# a2 = sum(1/i**2)
	b1 = (n + 1.0) / (3*(n-1))
	b2 = 2.0 * (n**2 + n + 3) / (9*n*(n-1.0))
	c1 = b1 - 1.0/a1
	c2 = b2 - (n + 2.0)/(a1*n) + a2/(a1**2)
	e1 = c1/a1
	e2 = c2/(a1**2 + a2)
	denom = cr_sqrt(e1*S + e2*S*(S-1))
	if denom != 0:
		D = (pi - theta) / denom
	else:
		D = -9
	return(D)


def compDiv(spA, spB, segsites):
	# spA and spB are list of haplotypes sequenced from species A and from species B
	div = [] # vector of divergence between A and B: length = number of pairwise comparisons. Contains the number of differences between individuals spA_i and spB_j
	nPair = 0	
	for i in spA:
		for j in spB:
			nPair += 1
			div.append(0)
			for k in range(segsites):
				if i[k]!=j[k]:
					div[nPair - 1] += 1
	res = {}
	res['div'] = cr_mean(div) # div = average of the raw divergence between individuals from A and those from B
	res['minDiv'] = min(div) # minimum div among all pairwises between A and B
	res['maxDiv'] = max(div)
	if res['div'] > 0:
		res['Gmin'] = res['minDiv']/res['div'] # Gmin = smallest divergence measured throughout all pairwise comp. / average div
		res['Gmax'] = res['maxDiv']/res['div']
	else:
		res['Gmin'] = -9
		res['Gmax'] = -9
	return(res)	


def nSs_nSf(ss, sf):
	# returns the number of loci with ss and sf, with ss but no sf, etc ...
	nLoci = len(ss)
	ss_sf, noSs_sf, ss_noSf, noSs_noSf = 0, 0, 0, 0
	for i in range(nLoci):
		if ss[i] == 0:
			if sf[i] == 0:
				noSs_noSf += 1
			if sf[i] > 0:
				noSs_sf += 1
		else:
			if sf[i] == 0:
				ss_noSf += 1
			if sf[i] > 0:
				ss_sf += 1
	res = {'ss_sf': ss_sf, 'ss_noSf': ss_noSf, 'noSs_sf': noSs_sf, 'noSs_noSf': noSs_noSf}
	return(res)


def ABBA_BABA_v2(pA, pB, pC):
	# pA, pB and pC are lists of allele frequencies
	nSNPs = len(pA)
	ABBAsum = 0.0
	BABAsum = 0.0
	maxABBAsumHom = 0.0
	maxBABAsumHom = 0.0
	maxABBAsumD = 0.0
	maxBABAsumD = 0.0
	#get derived frequencies for all biallelic siites
	for snp in range(nSNPs):
		p1Freq = pA[snp]
		p2Freq = pB[snp] 
		p3Freq = pC[snp]
		# get weigtings for ABBAs and BABAs
		try: # this was added to ignore crashes when there is missing data for a population at a site - we just ignore these sites
			ABBAsum += (1 - p1Freq) * p2Freq * p3Freq
			BABAsum += p1Freq * (1 - p2Freq) * p3Freq
			maxABBAsumHom += (1 - p1Freq) * p3Freq * p3Freq
			maxBABAsumHom += p1Freq * (1 - p3Freq) * p3Freq
			if p3Freq >= p2Freq:
				maxABBAsumD += (1 - p1Freq) * p3Freq * p3Freq 
				maxBABAsumD += p1Freq * (1 - p3Freq) * p3Freq
			else:
				maxABBAsumD += (1 - p1Freq) * p2Freq * p2Freq
				maxBABAsumD += p1Freq * (1 - p2Freq) * p2Freq
		except:
			continue
	#calculate D, f and fb
	output = {}
	try:
		output["D"] = (ABBAsum - BABAsum) / (ABBAsum + BABAsum)
	except:
		output["D"] = -9 
	try:
		output["fhom"] = (ABBAsum - BABAsum) / (maxABBAsumHom - maxBABAsumHom)
	except:
		output["fhom"] = -9 
	try:
		output["fd"] = (ABBAsum - BABAsum) / (maxABBAsumD - maxBABAsumD)
	except:
		output["fd"] = -9 
	output["ABBA"] = ABBAsum
	output["BABA"] = BABAsum
	return output


def ABBA_BABA(spA, spB, spC, spD):
	# spX is a list of haplotypes for the species X
	nS = len(spA[0]) # number of SNPs
	
	ABBA_sum_AB_C = 0.0
	ABBA_sum_BA_C = 0.0
	ABBA_sum_AB_D = 0.0
	ABBA_sum_BA_D = 0.0
	ABBA_sum_CD_A = 0.0
	ABBA_sum_DC_A = 0.0
	ABBA_sum_CD_B = 0.0
	ABBA_sum_DC_B = 0.0

	BABA_sum_AB_C = 0.0
	BABA_sum_BA_C = 0.0
	BABA_sum_AB_D = 0.0
	BABA_sum_BA_D = 0.0
	BABA_sum_CD_A = 0.0
	BABA_sum_DC_A = 0.0
	BABA_sum_CD_B = 0.0
	BABA_sum_DC_B = 0.0
	
	max_ABBA_sum_AB_C = 0.0
	max_ABBA_sum_BA_C = 0.0
	max_ABBA_sum_AB_D = 0.0
	max_ABBA_sum_BA_D = 0.0
	max_ABBA_sum_CD_A = 0.0
	max_ABBA_sum_DC_A = 0.0
	max_ABBA_sum_CD_B = 0.0
	max_ABBA_sum_DC_B = 0.0

	max_BABA_sum_AB_C = 0.0
	max_BABA_sum_BA_C = 0.0
	max_BABA_sum_AB_D = 0.0
	max_BABA_sum_BA_D = 0.0
	max_BABA_sum_CD_A = 0.0
	max_BABA_sum_DC_A = 0.0
	max_BABA_sum_CD_B = 0.0
	max_BABA_sum_DC_B = 0.0

	if nS >= 1:	
		# loop over SNPs
		for snp in range(nS):
			p1, p2, p3, p4 = 0.0, 0.0, 0.0, 0.0
			
			# loop over individuals
			for ind in range(len(spA)):
		#		print(spA[ind])
				if spA[ind][snp] == '1':
					p1 += 1
			for ind in range(len(spB)):
				if spB[ind][snp] == '1':
					p2 += 1
			for ind in range(len(spC)):
				if spC[ind][snp] == '1':
					p3 += 1
			for ind in range(len(spC)):
				if spD[ind][snp] == '1':
					p4 += 1
				
			# frequencies
			p1 /= float(len(spA))
			p2 /= float(len(spB))
			p3 /= float(len(spC))
			p4 /= float(len(spD))
			
			# ((A, B), C) to test admixture between B and C, or between A and C
			if (p1+p2+p3)>0 and (p1+p2+p3)<3:
				# D
				ABBA_sum_AB_C += (1-p1) * p2 * p3
				ABBA_sum_BA_C += (1-p2) * p1 * p3
				BABA_sum_AB_C += p1 * (1-p2) * p3
				BABA_sum_BA_C += p2 * (1-p1) * p3
				
				# fd
				## AB, C
				if p2 >= p3:
					max_ABBA_sum_AB_C += (1-p1) * p2 * p2
					max_BABA_sum_AB_C += p1 * (1-p2) * p2
				if p3 > p2:
					max_ABBA_sum_AB_C += (1-p1) * p3 * p3
					max_BABA_sum_AB_C += p1 * (1-p3) * p3
				
				## BA, C
				if p1 >= p3:
					max_ABBA_sum_BA_C += (1-p2) * p1 * p1
					max_BABA_sum_BA_C += p2 * (1-p1) * p1
				if p3 > p1:
					max_ABBA_sum_BA_C += (1-p1) * p3 * p3
					max_BABA_sum_BA_C += p1 * (1-p3) * p3
		
			# ((A, B), D) to test admixture between B and D, or between A and D
			if (p1+p2+p4)>0 and (p1+p2+p4)<3:
				# D
				ABBA_sum_AB_D += (1-p1) * p2 * p4
				ABBA_sum_BA_D += (1-p2) * p1 * p4
				BABA_sum_AB_D += p1 * (1-p2) * p4
				BABA_sum_BA_D += p2 * (1-p1) * p4
				
				# fd
				## AB, D
				if p2 >= p4:
					max_ABBA_sum_AB_D += (1-p1) * p2 * p2
					max_BABA_sum_AB_D += p1 * (1-p2) * p2
				if p4 > p2:
					max_ABBA_sum_AB_D += (1-p1) * p4 * p4
					max_BABA_sum_AB_D += p1 * (1-p4) * p4
				
				## BA, D
				if p1 >= p4:
					max_ABBA_sum_BA_D += (1-p2) * p1 * p1
					max_BABA_sum_BA_D += p2 * (1-p1) * p1
				if p4 > p1:
					max_ABBA_sum_BA_D += (1-p1) * p4 * p4
					max_BABA_sum_BA_D += p1 * (1-p4) * p4
		
			# ((C, D), A) to test admixture between D and A, or between C and A
			if (p3+p4+p1)>0 and (p3+p4+p1)<3:
				# D
				ABBA_sum_CD_A += (1-p3) * p4 * p1
				ABBA_sum_DC_A += (1-p4) * p3 * p1
				BABA_sum_CD_A += p3 * (1-p4) * p1
				BABA_sum_DC_A += p4 * (1-p3) * p1
				
				# fd
				## CD, A
				if p4 >= p1:
					max_ABBA_sum_CD_A += (1-p3) * p4 * p4
					max_BABA_sum_CD_A += p3 * (1-p4) * p4
				if p1 > p4:
					max_ABBA_sum_CD_A += (1-p3) * p1 * p1
					max_BABA_sum_CD_A += p3 * (1-p1) * p1
				
				## DC, A
				if p3 >= p1:
					max_ABBA_sum_DC_A += (1-p4) * p3 * p3
					max_BABA_sum_DC_A += p4 * (1-p3) * p3
				if p1 > p3:
					max_ABBA_sum_DC_A += (1-p4) * p1 * p1
					max_BABA_sum_DC_A += p4 * (1-p1) * p1
			
			# ((C, D), B) to test admixture between D and B, or between C and B
			if (p3+p4+p2)>0 and (p3+p4+p2)<3:
				# D
				ABBA_sum_CD_B += (1-p3) * p4 * p2
				ABBA_sum_DC_B += (1-p4) * p3 * p2
				BABA_sum_CD_B += p3 * (1-p4) * p2
				BABA_sum_DC_B += p4 * (1-p3) * p2
				
				# fd
				## CD, B
				if p4 >= p2:
					max_ABBA_sum_CD_B += (1-p3) * p4 * p4
					max_BABA_sum_CD_B += p3 * (1-p4) * p4
				if p2 > p4:
					max_ABBA_sum_CD_B += (1-p3) * p2 * p2
					max_BABA_sum_CD_B += p3 * (1-p2) * p2
				
				## DC, B
				if p3 >= p2:
					max_ABBA_sum_DC_B += (1-p4) * p3 * p3
					max_BABA_sum_DC_B += p4 * (1-p3) * p3
				if p2 > p3:
					max_ABBA_sum_DC_B += (1-p3) * p2 * p2
					max_BABA_sum_DC_B += p3 * (1-p2) * p2
		

		# AB_C
		## D
		if (ABBA_sum_AB_C + BABA_sum_AB_C) != 0:
			D_AB_C = (ABBA_sum_AB_C - BABA_sum_AB_C) / (ABBA_sum_AB_C + BABA_sum_AB_C)
		else:
			D_AB_C = -9
		## fd
		if (max_ABBA_sum_AB_C - max_BABA_sum_AB_C) != 0:
			fd_AB_C = (ABBA_sum_AB_C - BABA_sum_AB_C) / (max_ABBA_sum_AB_C - max_BABA_sum_AB_C)
		else:
			fd_AB_C = -9
		
		# BA_C
		## D
		if (ABBA_sum_BA_C + BABA_sum_BA_C) != 0:
			D_BA_C = (ABBA_sum_BA_C - BABA_sum_BA_C) / (ABBA_sum_BA_C + BABA_sum_BA_C)
		else:
			D_BA_C = -9
		## fd
		if (max_ABBA_sum_BA_C - max_BABA_sum_BA_C) != 0:
			fd_BA_C = (ABBA_sum_BA_C - BABA_sum_BA_C) / (max_ABBA_sum_BA_C - max_BABA_sum_BA_C)
		else:
			fd_BA_C = -9

		# AB_D
		## D
		if (ABBA_sum_AB_D + BABA_sum_AB_D) != 0:
			D_AB_D = (ABBA_sum_AB_D - BABA_sum_AB_D) / (ABBA_sum_AB_D + BABA_sum_AB_D)
		else:
			D_AB_D = -9
		## fd
		if (max_ABBA_sum_AB_D - max_BABA_sum_AB_D) != 0:
			fd_AB_D = (ABBA_sum_AB_D - BABA_sum_AB_D) / (max_ABBA_sum_AB_D - max_BABA_sum_AB_D)
		else:
			fd_AB_D = -9
		
		# BA_D
		## D
		if (ABBA_sum_BA_D + BABA_sum_BA_D) != 0:
			D_BA_D = (ABBA_sum_BA_D - BABA_sum_BA_D) / (ABBA_sum_BA_D + BABA_sum_BA_D)
		else:
			D_BA_D = -9
		## fd
		if (max_ABBA_sum_BA_D - max_BABA_sum_BA_D) != 0:
			fd_BA_D = (ABBA_sum_BA_D - BABA_sum_BA_D) / (max_ABBA_sum_BA_D - max_BABA_sum_BA_D)
		else:
			fd_BA_D = -9
		
		# CD_A
		## D
		if (ABBA_sum_CD_A + BABA_sum_CD_A) != 0:
			D_CD_A = (ABBA_sum_CD_A - BABA_sum_CD_A) / (ABBA_sum_CD_A + BABA_sum_CD_A)
		else:
			D_CD_A = -9
		## fd
		if (max_ABBA_sum_CD_A - max_BABA_sum_CD_A) != 0:
			fd_CD_A = (ABBA_sum_CD_A - BABA_sum_CD_A) / (max_ABBA_sum_CD_A - max_BABA_sum_CD_A)
		else:
			fd_CD_A = -9
		
		# DC_A
		## D
		if (ABBA_sum_DC_A + BABA_sum_DC_A) != 0:
			D_DC_A = (ABBA_sum_DC_A - BABA_sum_DC_A) / (ABBA_sum_DC_A + BABA_sum_DC_A)
		else:
			D_DC_A = -9
		## fd
		if (max_ABBA_sum_DC_A - max_BABA_sum_DC_A) != 0:
			fd_DC_A = (ABBA_sum_DC_A - BABA_sum_DC_A) / (max_ABBA_sum_DC_A - max_BABA_sum_DC_A)
		else:
			fd_DC_A = -9
		
		# CD_B
		## D
		if (ABBA_sum_CD_B + BABA_sum_CD_B) != 0:
			D_CD_B = (ABBA_sum_CD_B - BABA_sum_CD_B) / (ABBA_sum_CD_B + BABA_sum_CD_B)
		else:
			D_CD_B = -9
		## fd
		if (max_ABBA_sum_CD_B - max_BABA_sum_CD_B) != 0:
			fd_CD_B = (ABBA_sum_CD_B - BABA_sum_CD_B) / (max_ABBA_sum_CD_B - max_BABA_sum_CD_B)
		else:
			fd_CD_B = -9
		
		# DC_B
		## D
		if (ABBA_sum_DC_B + BABA_sum_DC_B) != 0:
			D_DC_B = (ABBA_sum_DC_B - BABA_sum_DC_B) / (ABBA_sum_DC_B + BABA_sum_DC_B)
		else:
			D_DC_B = -9
		## fd
		if (max_ABBA_sum_DC_B - max_BABA_sum_DC_B) != 0:
			fd_DC_B = (ABBA_sum_DC_B - BABA_sum_DC_B) / (max_ABBA_sum_DC_B - max_BABA_sum_DC_B)
		else:
			fd_DC_B = -9
		
		res = {}
		res['D_AB_C'] = D_AB_C
		res['fd_AB_C'] = fd_AB_C
		res['D_AB_D'] = D_AB_D
		res['fd_AB_D'] = fd_AB_D
		res['D_BA_C'] = D_BA_C
		res['fd_BA_C'] = fd_BA_C
		res['D_BA_D'] = D_BA_D
		res['fd_BA_D'] = fd_BA_D
		res['D_CD_A'] = D_CD_A
		res['fd_CD_A'] = fd_CD_A
		res['D_CD_B'] = D_CD_B
		res['fd_CD_B'] = fd_CD_B
		res['D_DC_A'] = D_DC_A
		res['fd_DC_A'] = fd_DC_A
		res['D_DC_B'] = D_DC_B
		res['fd_DC_B'] = fd_DC_B
	#	print(fd_DC_B)
	else:
		res = {}
		res['D_AB_C'] = -9
		res['fd_AB_C'] = -9
		res['D_AB_D'] = -9
		res['fd_AB_D'] = -9
		res['D_BA_C'] = -9
		res['fd_BA_C'] = -9
		res['D_BA_D'] = -9
		res['fd_BA_D'] = -9
		res['D_CD_A'] = -9
		res['fd_CD_A'] = -9
		res['D_CD_B'] = -9
		res['fd_CD_B'] = -9
		res['D_DC_A'] = -9
		res['fd_DC_A'] = -9
		res['D_DC_B'] = -9
		res['fd_DC_B'] = -9
	return(res)




def computeStats(locus, outfile):
        # read information about loci from the bpfile
        bpfile = locus + ".bpfile"
        msfile = locus + ".ms"

        contigName = '' # name of the contig where the locus was found
        posStart = ''
        posEnd = ''
        locusName = locus.split('_')
        locusName = locusName[4] + '_' + locusName[5]# + '_' + locusName[4]

        if os.path.isfile(bpfile) == False:
                sys.exit("\n\t\033[1;31;40mERROR: {0} was not found\n\033[0m\n".format(bpfile))


        infile = open(bpfile, "r")

	tmp = infile.readline() # first empty line
	L = [ float(i) for i in infile.readline().strip().split("\t") ]
	# (((A, B), (C, D)), outgroup)
	nSamA = [ int(i) for i in infile.readline().strip().split("\t") ]
	nSamB = [ int(i) for i in infile.readline().strip().split("\t") ]
	nSamC = [ int(i) for i in infile.readline().strip().split("\t") ]
	nSamD = [ int(i) for i in infile.readline().strip().split("\t") ]

	nLoci = len(L) 
	infile.close()
	
	a1_spA, a1_spB, a2_spA, a2_spB= [], [], [], []
	a1_spC, a1_spD, a2_spC, a2_spD= [], [], [], []
	for nsam in nSamA:
		a1_spA.append(sum([ 1.0/i for i in range(1, nsam) ]))
		a2_spA.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))
	for nsam in nSamB:
		a1_spB.append(sum([ 1.0/i for i in range(1, nsam) ]))
		a2_spB.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))
	for nsam in nSamC:
		a1_spC.append(sum([ 1.0/i for i in range(1, nsam) ]))
		a2_spC.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))
	for nsam in nSamD:
		a1_spD.append(sum([ 1.0/i for i in range(1, nsam) ]))
		a2_spD.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))


	#infile = open(msfile, "r")

	test = 0 
	nSim_cnt = 0 # count the number of treated multilocus simulations
	nLoci_cnt = 0 # count the number of treated loci within a simulation


	# READ THE MS's OUTPUTFILE
	##for line in infile:
	#for line in sys.stdin: # read the ms's output from the stdin
	infile = open(msfile, 'r')
	for line in infile:
		line = line.strip()
		if line[0:2] == '//':
			contigName = line.split(' ')[1]
		if "segsites" in line:
			if nLoci_cnt == 0:
				#ss_sf, noSs_sf, ss_noSf, noSs_noSf = 0, 0, 0, 0
				bialsites = []
				sfAB, sfAC, sfAD, sfBC, sfBD, sfCD = [], [], [], [], [], []
				sxA, sxB, sxC, sxD = [], [], [], []
				ssAB, ssAC, ssAD, ssBC, ssBD, ssCD = [], [], [], [], [], []
				#successive_ss = []
				piA, piB, piC, piD = [], [], [], []
				thetaA, thetaB, thetaC, thetaD = [], [], [], []
				DtajA, DtajB, DtajC, DtajD = [], [], [], []
				divAB, divAC, divAD, divBC, divBD, divCD = [], [], [], [], [], []
				netdivAB, netdivAC, netdivAD, netdivBC, netdivBD, netdivCD = [], [], [], [], [], []
				minDivAB, minDivAC, minDivAD, minDivBC, minDivBD, minDivCD = [], [], [], [], [], []
				maxDivAB, maxDivAC, maxDivAD, maxDivBC, maxDivBD, maxDivCD = [], [], [], [], [], []
				GminAB, GminAC, GminAD, GminBC, GminBD, GminCD = [], [], [], [], [], []
				GmaxAB, GmaxAC, GmaxAD, GmaxBC, GmaxBD, GmaxCD = [], [], [], [], [], []
				FST_AB, FST_AC, FST_AD, FST_BC, FST_BD, FST_CD = [], [], [], [], [], []
				D_AB_C, D_AB_D, D_BA_C, D_BA_D, D_CD_A, D_CD_B, D_DC_A, D_DC_B = [], [], [], [], [], [], [], []
				fd_AB_C, fd_AB_D, fd_BA_C, fd_BA_D, fd_CD_A, fd_CD_B, fd_DC_A, fd_DC_B = [], [], [], [], [], [], [], []
				fhom_AB_C, fhom_AB_D, fhom_BA_C, fhom_BA_D, fhom_CD_A, fhom_CD_B, fhom_DC_A, fhom_DC_B = [], [], [], [], [], [], [], []
			
			nLoci_cnt += 1
			nSam_cnt = 0 # count the number of treated individuals within a locus
			test = 1
			segsites = int(line.split(":")[1])
			bialsites.append(segsites)
			spA, spB, spC, spD = [], [], [], []
			continue
		if test == 1:
			if segsites == 0:
				test = 0
				posStart = 0
				posEnd = 0
				sfAB.append(0)
				sfAC.append(0)
				sfAD.append(0)
				sfBC.append(0)
				sfBD.append(0)
				sfCD.append(0)
				
				sxA.append(0)
				sxB.append(0)
				sxC.append(0)
				sxD.append(0)
				
				ssAB.append(0)
				ssAC.append(0)
				ssAD.append(0)
				ssBC.append(0)
				ssBD.append(0)
				ssCD.append(0)
				
				#successive_ss.append(0)
				piA.append(0)
				piB.append(0)
				piC.append(0)
				piD.append(0)
				
				thetaA.append(0)
				thetaB.append(0)
				thetaC.append(0)
				thetaD.append(0)
				
				DtajA.append(-9)
				DtajB.append(-9)
				DtajC.append(-9)
				DtajD.append(-9)
				
				divAB.append(0)
				divAC.append(0)
				divAD.append(0)
				divBC.append(0)
				divBD.append(0)
				divCD.append(0)
				
				netdivAB.append(0)
				netdivAC.append(0)
				netdivAD.append(0)
				netdivBC.append(0)
				netdivBD.append(0)
				netdivCD.append(0)
				
				minDivAB.append(0)
				minDivAC.append(0)
				minDivAD.append(0)
				minDivBC.append(0)
				minDivBD.append(0)
				minDivCD.append(0)
				
				maxDivAB.append(0)
				maxDivAC.append(0)
				maxDivAD.append(0)
				maxDivBC.append(0)
				maxDivBD.append(0)
				maxDivCD.append(0)
				
				GminAB.append(0)
				GminAC.append(0)
				GminAD.append(0)
				GminBC.append(0)
				GminBD.append(0)
				GminCD.append(0)
				
				FST_AB.append(-9)
				FST_AC.append(-9)
				FST_AD.append(-9)
				FST_BC.append(-9)
				FST_BD.append(-9)
				FST_CD.append(-9)
				
				D_AB_C.append(0)
				D_AB_D.append(0)
				D_BA_C.append(0)
				D_BA_D.append(0)
				D_CD_A.append(0)
				D_CD_B.append(0)
				D_DC_A.append(0)
				D_DC_B.append(0)
				
				fd_AB_C.append(0)
				fd_AB_D.append(0)
				fd_BA_C.append(0)
				fd_BA_D.append(0)
				fd_CD_A.append(0)
				fd_CD_B.append(0)
				fd_DC_A.append(0)
				fd_DC_B.append(0)
				
				fhom_AB_C.append(0)
				fhom_AB_D.append(0)
				fhom_BA_C.append(0)
				fhom_BA_D.append(0)
				fhom_CD_A.append(0)
				fhom_CD_B.append(0)
				fhom_DC_A.append(0)
				fhom_DC_B.append(0)
				#noSs_noSf += 1
			if segsites != 0:
				if "positions" in line:
					tmp = line.strip().split(' ')
					posStart = tmp[1]
					posEnd = tmp[-1]
				if "positions" not in line and line!="\n":
					nSam_cnt += 1
					if nSam_cnt <= nSamA[nLoci_cnt - 1]:
						spA.append(line.strip())
					if nSam_cnt > nSamA[nLoci_cnt - 1] and nSam_cnt <= (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1]):
						spB.append(line.strip())
					if nSam_cnt > (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1]) and nSam_cnt <= (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1] + nSamC[nLoci_cnt - 1]):
						spC.append(line.strip())
					if nSam_cnt > (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1] + nSamC[nLoci_cnt - 1]) and nSam_cnt <= (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1] + nSamC[nLoci_cnt - 1] + nSamD[nLoci_cnt - 1]):
						spD.append(line.strip())
					
					# end of the block of sequences -> start the calculations of various statistics
					if nSam_cnt == (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1] + nSamC[nLoci_cnt - 1] + nSamD[nLoci_cnt - 1]):
						#print("locus {0}".format(nLoci_cnt))
						#print("spA\n{0}\n\nspB\n{1}\n\nspC\n{2}\n\nspD\n{3}\n".format("\n".join(spA), "\n".join(spB), "\n".join(spC), "\n".join(spD)))
						tmpA = compFreq(spA, segsites)
						freqA = tmpA['freq']
						piA.append(tmpA['pi']/L[nLoci_cnt - 1])
						
						tmpB = compFreq(spB, segsites)
						freqB = tmpB['freq']
						piB.append(tmpB['pi']/L[nLoci_cnt - 1])
						
						tmpC = compFreq(spC, segsites)
						freqC = tmpC['freq']
						piC.append(tmpC['pi']/L[nLoci_cnt - 1])
						
						tmpD = compFreq(spD, segsites)
						freqD = tmpD['freq']
						piD.append(tmpD['pi']/L[nLoci_cnt - 1])
						
						## Sx
						tmpABC = compFreq(spA+spB+spC, segsites)
						freqABC = tmpABC['freq']
						tmpABD = compFreq(spA+spB+spD, segsites)
						freqABD = tmpABD['freq']
						tmpACD = compFreq(spA+spC+spD, segsites)
						freqACD = tmpACD['freq']
						tmpBCD = compFreq(spB+spC+spD, segsites)
						freqBCD = tmpBCD['freq']
						
						# SxA
						tmp_A_BCD = sites(freqA, freqBCD, segsites)
						sxA.append(tmp_A_BCD['sxA']/(1.0*L[nLoci_cnt - 1]))
						# SxB
						tmp_B_ACD = sites(freqB, freqACD, segsites)
						sxB.append(tmp_B_ACD['sxA']/(1.0*L[nLoci_cnt - 1]))
						# SxC
						tmp_C_ABD = sites(freqC, freqABD, segsites)
						sxC.append(tmp_C_ABD['sxA']/(1.0*L[nLoci_cnt - 1]))
						# SxD
						tmp_D_ABC = sites(freqD, freqABC, segsites)
						sxD.append(tmp_D_ABC['sxA']/(1.0*L[nLoci_cnt - 1]))
						
						# AB
						tmp_AB = sites(freqA, freqB, segsites) # determines the # of sxA, sxB, ss, sf
						sfAB.append(tmp_AB['sf']/(1.0*L[nLoci_cnt - 1]))
						#sxA.append(tmp_AB['sxA']/(1.0*L[nLoci_cnt - 1]))
						#sxB.append(tmp_AB['sxB']/(1.0*L[nLoci_cnt - 1]))
						ssAB.append(tmp_AB['ss']/(1.0*L[nLoci_cnt - 1]))

						# AC
						tmp_AC = sites(freqA, freqC, segsites) # determines the # of sxA, sxC, ss, sf
						sfAC.append(tmp_AC['sf']/(1.0*L[nLoci_cnt - 1]))
						ssAC.append(tmp_AC['ss']/(1.0*L[nLoci_cnt - 1]))
						
						# AD
						tmp_AD = sites(freqA, freqD, segsites) # determines the # of sxA, sxD, ss, sf
						sfAD.append(tmp_AD['sf']/(1.0*L[nLoci_cnt - 1]))
						ssAD.append(tmp_AD['ss']/(1.0*L[nLoci_cnt - 1]))

						# BC
						tmp_BC = sites(freqB, freqC, segsites) # determines the # of sxB, sxC, ss, sf
						sfBC.append(tmp_BC['sf']/(1.0*L[nLoci_cnt - 1]))
						ssBC.append(tmp_BC['ss']/(1.0*L[nLoci_cnt - 1]))
						
						# BD
						tmp_BD = sites(freqB, freqD, segsites) # determines the # of sxB, sxD, ss, sf
						sfBD.append(tmp_BD['sf']/(1.0*L[nLoci_cnt - 1]))
						ssBD.append(tmp_BD['ss']/(1.0*L[nLoci_cnt - 1]))
						
						# CD
						tmp_CD = sites(freqC, freqD, segsites) # determines the # of sxC, sxD, ss, sf
						sfCD.append(tmp_CD['sf']/(1.0*L[nLoci_cnt - 1]))
						ssCD.append(tmp_CD['ss']/(1.0*L[nLoci_cnt - 1]))
						
						#successive_ss.append(tmp['successive_ss'])
						# test if loci has both ss, sf, or only one of two, etc ...
						#if tmp['sf'] > 0:
						#	if tmp['ss'] > 0:
						#		ss_sf += 1
						#	if tmp['ss'] == 0:
						#		noSs_sf += 1
						#if tmp['sf'] == 0:
						#	if tmp['ss'] > 0:
						#		ss_noSf += 1
						#	if tmp['ss'] == 0:
						#		noSs_noSf += 1
						
						# theta = sx + ss
	#					thetaA_locus = (tmp['sxA']+tmp['ss'])/a1_spA[nLoci_cnt - 1]
	#					thetaB_locus = (tmp['sxB']+tmp['ss'])/a1_spB[nLoci_cnt - 1]
						thetaA_locus = (tmpA['nSNPs'])/a1_spA[nLoci_cnt - 1]
						thetaB_locus = (tmpB['nSNPs'])/a1_spB[nLoci_cnt - 1]
						thetaC_locus = (tmpC['nSNPs'])/a1_spC[nLoci_cnt - 1]
						thetaD_locus = (tmpD['nSNPs'])/a1_spD[nLoci_cnt - 1]
						thetaA.append(thetaA_locus/L[nLoci_cnt - 1 ])
						thetaB.append(thetaB_locus/L[nLoci_cnt - 1 ])
						thetaC.append(thetaC_locus/L[nLoci_cnt - 1 ])
						thetaD.append(thetaD_locus/L[nLoci_cnt - 1 ])
						
						# Taj
	#					DtajA.append(tajD(tmpA['pi'], thetaA_locus, nSamA[nLoci_cnt - 1], tmp['sxA']+tmp['ss'], a1_spA[nLoci_cnt-1], a2_spA[nLoci_cnt-1]))
	#					DtajB.append(tajD(tmpB['pi'], thetaB_locus, nSamB[nLoci_cnt - 1], tmp['sxB']+tmp['ss'], a1_spB[nLoci_cnt-1], a2_spB[nLoci_cnt-1]))
	#					DtajC.append(tajD(tmpC['pi'], thetaC_locus, nSamC[nLoci_cnt - 1], tmp['sxC']+tmp['ss'], a1_spC[nLoci_cnt-1], a2_spC[nLoci_cnt-1]))
	#					DtajD.append(tajD(tmpD['pi'], thetaD_locus, nSamD[nLoci_cnt - 1], tmp['sxD']+tmp['ss'], a1_spD[nLoci_cnt-1], a2_spD[nLoci_cnt-1]))
						DtajA.append(tajD(tmpA['pi'], thetaA_locus, nSamA[nLoci_cnt - 1], tmpA['nSNPs'], a1_spA[nLoci_cnt-1], a2_spA[nLoci_cnt-1]))
						DtajB.append(tajD(tmpB['pi'], thetaB_locus, nSamB[nLoci_cnt - 1], tmpB['nSNPs'], a1_spB[nLoci_cnt-1], a2_spB[nLoci_cnt-1]))
						DtajC.append(tajD(tmpC['pi'], thetaC_locus, nSamC[nLoci_cnt - 1], tmpC['nSNPs'], a1_spC[nLoci_cnt-1], a2_spC[nLoci_cnt-1]))
						DtajD.append(tajD(tmpD['pi'], thetaD_locus, nSamD[nLoci_cnt - 1], tmpD['nSNPs'], a1_spD[nLoci_cnt-1], a2_spD[nLoci_cnt-1]))
						## divergences AB, AC, AD, BC, BD, CD					
						# divAB
						divAB_tmp = compDiv(spA, spB, segsites)
						divAB.append(divAB_tmp['div']/L[nLoci_cnt - 1 ])
						netdivAB.append(divAB_tmp['div']/L[nLoci_cnt - 1 ] - (tmpA['pi']/L[nLoci_cnt - 1] + tmpB['pi']/L[nLoci_cnt - 1]) / 2.0)
						minDivAB.append(divAB_tmp['minDiv']/L[nLoci_cnt - 1 ])
						maxDivAB.append(divAB_tmp['maxDiv']/L[nLoci_cnt - 1 ])
						GminAB.append(divAB_tmp['Gmin'])
						GmaxAB.append(divAB_tmp['Gmax'])
						
						# divAC
						divAC_tmp = compDiv(spA, spC, segsites)
						divAC.append(divAC_tmp['div']/L[nLoci_cnt - 1 ])
						netdivAC.append(divAC_tmp['div']/L[nLoci_cnt - 1 ] - (tmpA['pi']/L[nLoci_cnt - 1] + tmpC['pi']/L[nLoci_cnt - 1]) / 2.0)
						minDivAC.append(divAC_tmp['minDiv']/L[nLoci_cnt - 1 ])
						maxDivAC.append(divAC_tmp['maxDiv']/L[nLoci_cnt - 1 ])
						GminAC.append(divAC_tmp['Gmin'])
						GmaxAC.append(divAC_tmp['Gmax'])
						
						# divAD
						divAD_tmp = compDiv(spA, spD, segsites)
						divAD.append(divAD_tmp['div']/L[nLoci_cnt - 1 ])
						netdivAD.append(divAD_tmp['div']/L[nLoci_cnt - 1 ] - (tmpA['pi']/L[nLoci_cnt - 1] + tmpD['pi']/L[nLoci_cnt - 1]) / 2.0)
						minDivAD.append(divAD_tmp['minDiv']/L[nLoci_cnt - 1 ])
						maxDivAD.append(divAD_tmp['maxDiv']/L[nLoci_cnt - 1 ])
						GminAD.append(divAD_tmp['Gmin'])
						GmaxAD.append(divAD_tmp['Gmax'])
						
						# divBC
						divBC_tmp = compDiv(spB, spC, segsites)
						divBC.append(divBC_tmp['div']/L[nLoci_cnt - 1 ])
						netdivBC.append(divBC_tmp['div']/L[nLoci_cnt - 1 ] - (tmpB['pi']/L[nLoci_cnt - 1] + tmpC['pi']/L[nLoci_cnt - 1]) / 2.0)
						minDivBC.append(divBC_tmp['minDiv']/L[nLoci_cnt - 1 ])
						maxDivBC.append(divBC_tmp['maxDiv']/L[nLoci_cnt - 1 ])
						GminBC.append(divBC_tmp['Gmin'])
						GmaxBC.append(divBC_tmp['Gmax'])
						
						# divBD
						divBD_tmp = compDiv(spB, spD, segsites)
						divBD.append(divBD_tmp['div']/L[nLoci_cnt - 1 ])
						netdivBD.append(divBD_tmp['div']/L[nLoci_cnt - 1 ] - (tmpB['pi']/L[nLoci_cnt - 1] + tmpD['pi']/L[nLoci_cnt - 1]) / 2.0)
						minDivBD.append(divBD_tmp['minDiv']/L[nLoci_cnt - 1 ])
						maxDivBD.append(divBD_tmp['maxDiv']/L[nLoci_cnt - 1 ])
						GminBD.append(divBD_tmp['Gmin'])
						GmaxBD.append(divBD_tmp['Gmax'])
						
						# divCD
						divCD_tmp = compDiv(spC, spD, segsites)
						divCD.append(divCD_tmp['div']/L[nLoci_cnt - 1 ])
						netdivCD.append(divCD_tmp['div']/L[nLoci_cnt - 1 ] - (tmpC['pi']/L[nLoci_cnt - 1] + tmpD['pi']/L[nLoci_cnt - 1]) / 2.0)
						minDivCD.append(divCD_tmp['minDiv']/L[nLoci_cnt - 1 ])
						maxDivCD.append(divCD_tmp['maxDiv']/L[nLoci_cnt - 1 ])
						GminCD.append(divCD_tmp['Gmin'])
						GmaxCD.append(divCD_tmp['Gmax'])
						
						# Fst
						# vector of piT over segsites
						piT_AB = piTot(tmpA['nDer'], tmpB['nDer'], nSamA[nLoci_cnt - 1], nSamB[nLoci_cnt - 1], segsites)
						piT_AC = piTot(tmpA['nDer'], tmpC['nDer'], nSamA[nLoci_cnt - 1], nSamC[nLoci_cnt - 1], segsites)
						piT_AD = piTot(tmpA['nDer'], tmpD['nDer'], nSamA[nLoci_cnt - 1], nSamD[nLoci_cnt - 1], segsites)
						piT_BC = piTot(tmpB['nDer'], tmpC['nDer'], nSamB[nLoci_cnt - 1], nSamC[nLoci_cnt - 1], segsites)
						piT_BD = piTot(tmpB['nDer'], tmpD['nDer'], nSamB[nLoci_cnt - 1], nSamD[nLoci_cnt - 1], segsites)
						piT_CD = piTot(tmpC['nDer'], tmpD['nDer'], nSamC[nLoci_cnt - 1], nSamD[nLoci_cnt - 1], segsites)
						
						# mean Fst
						FST_AB.append(Fst(cr_sum(tmpA['pi_SNPs']), cr_sum(tmpB['pi_SNPs']), cr_sum(piT_AB)))
						FST_AC.append(Fst(cr_sum(tmpA['pi_SNPs']), cr_sum(tmpC['pi_SNPs']), cr_sum(piT_AC)))
						FST_AD.append(Fst(cr_sum(tmpA['pi_SNPs']), cr_sum(tmpD['pi_SNPs']), cr_sum(piT_AD)))
						FST_BC.append(Fst(cr_sum(tmpB['pi_SNPs']), cr_sum(tmpC['pi_SNPs']), cr_sum(piT_BC)))
						FST_BD.append(Fst(cr_sum(tmpB['pi_SNPs']), cr_sum(tmpD['pi_SNPs']), cr_sum(piT_BD)))
						FST_CD.append(Fst(cr_sum(tmpC['pi_SNPs']), cr_sum(tmpD['pi_SNPs']), cr_sum(piT_CD)))
						
						# ABBA_BABA
						#ABBA_BABA_res = ABBA_BABA(spA, spB, spC, spD)
						#D_AB_C.append(ABBA_BABA_res['D_AB_C'])
						#D_AB_D.append(ABBA_BABA_res['D_AB_D'])
						#D_BA_C.append(ABBA_BABA_res['D_BA_C'])
						#D_BA_D.append(ABBA_BABA_res['D_BA_D'])
						#D_CD_A.append(ABBA_BABA_res['D_CD_A'])
						#D_CD_B.append(ABBA_BABA_res['D_CD_B'])
						#D_DC_A.append(ABBA_BABA_res['D_DC_A'])
						#D_DC_B.append(ABBA_BABA_res['D_DC_B'])

						#fd_AB_C.append(ABBA_BABA_res['fd_AB_C'])
						#fd_AB_D.append(ABBA_BABA_res['fd_AB_D'])
						#fd_BA_C.append(ABBA_BABA_res['fd_BA_C'])
						#fd_BA_D.append(ABBA_BABA_res['fd_BA_D'])
						#fd_CD_A.append(ABBA_BABA_res['fd_CD_A'])
						#fd_CD_B.append(ABBA_BABA_res['fd_CD_B'])
						#fd_DC_A.append(ABBA_BABA_res['fd_DC_A'])
						#fd_DC_B.append(ABBA_BABA_res['fd_DC_B'])
						
						ABBA_BABA_AB_C = ABBA_BABA_v2(tmpA['freq'], tmpB['freq'], tmpC['freq'])
						ABBA_BABA_AB_D = ABBA_BABA_v2(tmpA['freq'], tmpB['freq'], tmpD['freq'])
						ABBA_BABA_BA_C = ABBA_BABA_v2(tmpB['freq'], tmpA['freq'], tmpC['freq'])
						ABBA_BABA_BA_D = ABBA_BABA_v2(tmpB['freq'], tmpA['freq'], tmpD['freq'])
						ABBA_BABA_CD_A = ABBA_BABA_v2(tmpC['freq'], tmpD['freq'], tmpA['freq'])
						ABBA_BABA_CD_B = ABBA_BABA_v2(tmpC['freq'], tmpD['freq'], tmpB['freq'])
						ABBA_BABA_DC_A = ABBA_BABA_v2(tmpD['freq'], tmpC['freq'], tmpA['freq'])
						ABBA_BABA_DC_B = ABBA_BABA_v2(tmpD['freq'], tmpC['freq'], tmpB['freq'])
						
						D_AB_C.append(ABBA_BABA_AB_C['D'])
						D_AB_D.append(ABBA_BABA_AB_D['D'])
						D_BA_C.append(ABBA_BABA_BA_C['D'])
						D_BA_D.append(ABBA_BABA_BA_D['D'])
						D_CD_A.append(ABBA_BABA_CD_A['D'])
						D_CD_B.append(ABBA_BABA_CD_B['D'])
						D_DC_A.append(ABBA_BABA_DC_A['D'])
						D_DC_B.append(ABBA_BABA_DC_B['D'])

						fd_AB_C.append(ABBA_BABA_AB_C['fd'])
						fd_AB_D.append(ABBA_BABA_AB_D['fd'])
						fd_BA_C.append(ABBA_BABA_BA_C['fd'])
						fd_BA_D.append(ABBA_BABA_BA_D['fd'])
						fd_CD_A.append(ABBA_BABA_CD_A['fd'])
						fd_CD_B.append(ABBA_BABA_CD_B['fd'])
						fd_DC_A.append(ABBA_BABA_DC_A['fd'])
						fd_DC_B.append(ABBA_BABA_DC_B['fd'])

						fhom_AB_C.append(ABBA_BABA_AB_C['fhom'])
						fhom_AB_D.append(ABBA_BABA_AB_D['fhom'])
						fhom_BA_C.append(ABBA_BABA_BA_C['fhom'])
						fhom_BA_D.append(ABBA_BABA_BA_D['fhom'])
						fhom_CD_A.append(ABBA_BABA_CD_A['fhom'])
						fhom_CD_B.append(ABBA_BABA_CD_B['fhom'])
						fhom_DC_A.append(ABBA_BABA_DC_A['fhom'])
						fhom_DC_B.append(ABBA_BABA_DC_B['fhom'])
		
		# compute average and std over of statistics over loci
		if nLoci_cnt != 0 and len(sxA) == nLoci:
			test = 0
			nSim_cnt += 1
			nLoci_cnt = 0
			
			# statistics
			bialsites_avg = cr_mean(bialsites)
			bialsites_std = cr_std(bialsites, bialsites_avg)
			sfAB_avg = cr_mean(sfAB)
			sfAB_std = cr_std(sfAB, sfAB_avg)
			sfAC_avg = cr_mean(sfAC)
			sfAC_std = cr_std(sfAC, sfAC_avg)
			sfAD_avg = cr_mean(sfAD)
			sfAD_std = cr_std(sfAD, sfAD_avg)
			sfBC_avg = cr_mean(sfBC)
			sfBC_std = cr_std(sfBC, sfBC_avg)
			sfBD_avg = cr_mean(sfBD)
			sfBD_std = cr_std(sfBD, sfBD_avg)
			sfCD_avg = cr_mean(sfCD)
			sfCD_std = cr_std(sfCD, sfCD_avg)
			
			sxA_avg = cr_mean(sxA)
			sxA_std = cr_std(sxA, sxA_avg)
			sxB_avg = cr_mean(sxB)
			sxB_std = cr_std(sxB, sxB_avg)
			sxC_avg = cr_mean(sxC)
			sxC_std = cr_std(sxC, sxC_avg)
			sxD_avg = cr_mean(sxD)
			sxD_std = cr_std(sxD, sxD_avg)
		
			ssAB_avg = cr_mean(ssAB)
			ssAB_std = cr_std(ssAB, ssAB_avg)
			ssAC_avg = cr_mean(ssAC)
			ssAC_std = cr_std(ssAC, ssAC_avg)
			ssAD_avg = cr_mean(ssAD)
			ssAD_std = cr_std(ssAD, ssAD_avg)
			ssBC_avg = cr_mean(ssBC)
			ssBC_std = cr_std(ssBC, ssBC_avg)
			ssBD_avg = cr_mean(ssBD)
			ssBD_std = cr_std(ssBD, ssBD_avg)
			ssCD_avg = cr_mean(ssCD)
			ssCD_std = cr_std(ssCD, ssCD_avg)
			
			#successive_ss_avg = cr_mean(successive_ss)
			#successive_ss_std = cr_std(successive_ss, successive_ss_avg)
			
			piA_avg = cr_mean(piA)
			piA_std = cr_std(piA, piA_avg)
			piB_avg = cr_mean(piB)
			piB_std = cr_std(piB, piB_avg)
			piC_avg = cr_mean(piC)
			piC_std = cr_std(piC, piC_avg)
			piD_avg = cr_mean(piD)
			piD_std = cr_std(piD, piD_avg)
			
			pearson_r_pi_AB = cr_pearsonR(piA, piB)
			pearson_r_pi_AC = cr_pearsonR(piA, piC)
			pearson_r_pi_AD = cr_pearsonR(piA, piD)
			pearson_r_pi_BC = cr_pearsonR(piB, piC)
			pearson_r_pi_BD = cr_pearsonR(piB, piD)
			pearson_r_pi_CD = cr_pearsonR(piC, piD)
			
			thetaA_avg = cr_mean(thetaA)
			thetaA_std = cr_std(thetaA, thetaA_avg)
			thetaB_avg = cr_mean(thetaB)
			thetaB_std = cr_std(thetaB, thetaB_avg)
			thetaC_avg = cr_mean(thetaC)
			thetaC_std = cr_std(thetaC, thetaC_avg)
			thetaD_avg = cr_mean(thetaD)
			thetaD_std = cr_std(thetaD, thetaD_avg)
			
			pearson_r_theta_AB = cr_pearsonR(thetaA, thetaB)
			pearson_r_theta_AC = cr_pearsonR(thetaA, thetaC)
			pearson_r_theta_AD = cr_pearsonR(thetaA, thetaD)
			pearson_r_theta_BC = cr_pearsonR(thetaB, thetaC)
			pearson_r_theta_BD = cr_pearsonR(thetaB, thetaD)
			pearson_r_theta_CD = cr_pearsonR(thetaC, thetaD)
			
			DtajA_avg = cr_mean(DtajA)
			DtajA_std = cr_std(DtajA, DtajA_avg)
			DtajB_avg = cr_mean(DtajB)
			DtajB_std = cr_std(DtajB, DtajB_avg)
			DtajC_avg = cr_mean(DtajC)
			DtajC_std = cr_std(DtajC, DtajC_avg)
			DtajD_avg = cr_mean(DtajD)
			DtajD_std = cr_std(DtajD, DtajD_avg)
			
			divAB_avg = cr_mean(divAB)
			divAB_std = cr_std(divAB, divAB_avg)
			divAC_avg = cr_mean(divAC)
			divAC_std = cr_std(divAC, divAC_avg)
			divAD_avg = cr_mean(divAD)
			divAD_std = cr_std(divAD, divAD_avg)
			divBC_avg = cr_mean(divBC)
			divBC_std = cr_std(divBC, divBC_avg)
			divBD_avg = cr_mean(divBD)
			divBD_std = cr_std(divBD, divBD_avg)
			divCD_avg = cr_mean(divCD)
			divCD_std = cr_std(divCD, divCD_avg)
			
			netdivAB_avg = cr_mean(netdivAB)
			netdivAB_std = cr_std(netdivAB, netdivAB_avg)
			netdivAC_avg = cr_mean(netdivAC)
			netdivAC_std = cr_std(netdivAC, netdivAC_avg)
			netdivAD_avg = cr_mean(netdivAD)
			netdivAD_std = cr_std(netdivAD, netdivAD_avg)
			netdivBC_avg = cr_mean(netdivBC)
			netdivBC_std = cr_std(netdivBC, netdivBC_avg)
			netdivBD_avg = cr_mean(netdivBD)
			netdivBD_std = cr_std(netdivBD, netdivBD_avg)
			netdivCD_avg = cr_mean(netdivCD)
			netdivCD_std = cr_std(netdivCD, netdivCD_avg)
			
			minDivAB_avg = cr_mean(minDivAB)
			minDivAB_std = cr_std(minDivAB, minDivAB_avg)
			minDivAC_avg = cr_mean(minDivAC)
			minDivAC_std = cr_std(minDivAC, minDivAC_avg)
			minDivAD_avg = cr_mean(minDivAD)
			minDivAD_std = cr_std(minDivAD, minDivAD_avg)
			minDivBC_avg = cr_mean(minDivBC)
			minDivBC_std = cr_std(minDivBC, minDivBC_avg)
			minDivBD_avg = cr_mean(minDivBD)
			minDivBD_std = cr_std(minDivBD, minDivBD_avg)
			minDivCD_avg = cr_mean(minDivCD)
			minDivCD_std = cr_std(minDivCD, minDivCD_avg)

			maxDivAB_avg = cr_mean(maxDivAB)
			maxDivAB_std = cr_std(maxDivAB, maxDivAB_avg)
			maxDivAC_avg = cr_mean(maxDivAC)
			maxDivAC_std = cr_std(maxDivAC, maxDivAC_avg)
			maxDivAD_avg = cr_mean(maxDivAD)
			maxDivAD_std = cr_std(maxDivAD, maxDivAD_avg)
			maxDivBC_avg = cr_mean(maxDivBC)
			maxDivBC_std = cr_std(maxDivBC, maxDivBC_avg)
			maxDivBD_avg = cr_mean(maxDivBD)
			maxDivBD_std = cr_std(maxDivBD, maxDivBD_avg)
			maxDivCD_avg = cr_mean(maxDivCD)
			maxDivCD_std = cr_std(maxDivCD, maxDivCD_avg)
			
			GminAB_avg = cr_mean(GminAB)
			GminAB_std = cr_std(GminAB, GminAB_avg)
			GminAC_avg = cr_mean(GminAC)
			GminAC_std = cr_std(GminAC, GminAC_avg)
			GminAD_avg = cr_mean(GminAD)
			GminAD_std = cr_std(GminAD, GminAD_avg)
			GminBC_avg = cr_mean(GminBC)
			GminBC_std = cr_std(GminBC, GminBC_avg)
			GminBD_avg = cr_mean(GminBD)
			GminBD_std = cr_std(GminBD, GminBD_avg)
			GminCD_avg = cr_mean(GminCD)
			GminCD_std = cr_std(GminCD, GminCD_avg)

			GmaxAB_avg = cr_mean(GmaxAB)
			GmaxAB_std = cr_std(GmaxAB, GmaxAB_avg)
			GmaxAC_avg = cr_mean(GmaxAC)
			GmaxAC_std = cr_std(GmaxAC, GmaxAC_avg)
			GmaxAD_avg = cr_mean(GmaxAD)
			GmaxAD_std = cr_std(GmaxAD, GmaxAD_avg)
			GmaxBC_avg = cr_mean(GmaxBC)
			GmaxBC_std = cr_std(GmaxBC, GmaxBC_avg)
			GmaxBD_avg = cr_mean(GmaxBD)
			GmaxBD_std = cr_std(GmaxBD, GmaxBD_avg)
			GmaxCD_avg = cr_mean(GmaxCD)
			GmaxCD_std = cr_std(GmaxCD, GmaxCD_avg)
			
			FST_AB_avg = cr_mean(FST_AB)
			FST_AB_std = cr_std(FST_AB, FST_AB_avg)
			FST_AC_avg = cr_mean(FST_AC)
			FST_AC_std = cr_std(FST_AC, FST_AC_avg)
			FST_AD_avg = cr_mean(FST_AD)
			FST_AD_std = cr_std(FST_AD, FST_AD_avg)
			FST_BC_avg = cr_mean(FST_BC)
			FST_BC_std = cr_std(FST_BC, FST_BC_avg)
			FST_BD_avg = cr_mean(FST_BD)
			FST_BD_std = cr_std(FST_BD, FST_BD_avg)
			FST_CD_avg = cr_mean(FST_CD)
			FST_CD_std = cr_std(FST_CD, FST_CD_avg)
			
			D_AB_C_avg = cr_mean(D_AB_C)
			D_AB_C_std = cr_std(D_AB_C, D_AB_C_avg)
			D_AB_D_avg = cr_mean(D_AB_D)
			D_AB_D_std = cr_std(D_AB_D, D_AB_D_avg)
			D_BA_C_avg = cr_mean(D_BA_C)
			D_BA_C_std = cr_std(D_BA_C, D_BA_C_avg)
			D_BA_D_avg = cr_mean(D_BA_D)
			D_BA_D_std = cr_std(D_BA_D, D_BA_D_avg)
			D_CD_A_avg = cr_mean(D_CD_A)
			D_CD_A_std = cr_std(D_CD_A, D_CD_A_avg)
			D_CD_B_avg = cr_mean(D_CD_B)
			D_CD_B_std = cr_std(D_CD_B, D_CD_B_avg)
			D_DC_A_avg = cr_mean(D_DC_A)
			D_DC_A_std = cr_std(D_DC_A, D_DC_A_avg)
			D_DC_B_avg = cr_mean(D_DC_B)
			D_DC_B_std = cr_std(D_DC_B, D_DC_B_avg)
			
			fd_AB_C_avg = cr_mean(fd_AB_C)
			fd_AB_C_std = cr_std(fd_AB_C, fd_AB_C_avg)
			fd_AB_D_avg = cr_mean(fd_AB_D)
			fd_AB_D_std = cr_std(fd_AB_D, fd_AB_D_avg)
			fd_BA_C_avg = cr_mean(fd_BA_C)
			fd_BA_C_std = cr_std(fd_BA_C, fd_BA_C_avg)
			fd_BA_D_avg = cr_mean(fd_BA_D)
			fd_BA_D_std = cr_std(fd_BA_D, fd_BA_D_avg)
			fd_CD_A_avg = cr_mean(fd_CD_A)
			fd_CD_A_std = cr_std(fd_CD_A, fd_CD_A_avg)
			fd_CD_B_avg = cr_mean(fd_CD_B)
			fd_CD_B_std = cr_std(fd_CD_B, fd_CD_B_avg)
			fd_DC_A_avg = cr_mean(fd_DC_A)
			fd_DC_A_std = cr_std(fd_DC_A, fd_DC_A_avg)
			fd_DC_B_avg = cr_mean(fd_DC_B)
			fd_DC_B_std = cr_std(fd_DC_B, fd_DC_B_avg)
			
			fhom_AB_C_avg = cr_mean(fhom_AB_C)
			fhom_AB_C_std = cr_std(fhom_AB_C, fhom_AB_C_avg)
			fhom_AB_D_avg = cr_mean(fhom_AB_D)
			fhom_AB_D_std = cr_std(fhom_AB_D, fhom_AB_D_avg)
			fhom_BA_C_avg = cr_mean(fhom_BA_C)
			fhom_BA_C_std = cr_std(fhom_BA_C, fhom_BA_C_avg)
			fhom_BA_D_avg = cr_mean(fhom_BA_D)
			fhom_BA_D_std = cr_std(fhom_BA_D, fhom_BA_D_avg)
			fhom_CD_A_avg = cr_mean(fhom_CD_A)
			fhom_CD_A_std = cr_std(fhom_CD_A, fhom_CD_A_avg)
			fhom_CD_B_avg = cr_mean(fhom_CD_B)
			fhom_CD_B_std = cr_std(fhom_CD_B, fhom_CD_B_avg)
			fhom_DC_A_avg = cr_mean(fhom_DC_A)
			fhom_DC_A_std = cr_std(fhom_DC_A, fhom_DC_A_avg)
			fhom_DC_B_avg = cr_mean(fhom_DC_B)
			fhom_DC_B_std = cr_std(fhom_DC_B, fhom_DC_B_avg)
			
			#pearson_r_div_netDiv = cr_pearsonR(divAB, netdivAB)
			#pearson_r_div_FST = cr_pearsonR(divAB, FST)
			#pearson_r_netDiv_FST = cr_pearsonR(netdivAB, FST)
			
			#print("dataset {0}: {1} loci".format(nSim_cnt-1, len(ss)))
			#res += "{0}\t{1:.5f}\t{2:.5f}\t".format(nSim_cnt-1, bialsites_avg, bialsites_std)
			res = "{0}\t{1}\t{2}\t{3}\t{4}\t".format(contigName, posStart, posEnd, locusName, L[0]) # contig\tposStart\tposEnd\tname\tlength
			res += "{0:.5f}\t{1:.5f}\t".format(bialsites_avg, bialsites_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sfAB_avg, sfAB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sfAC_avg, sfAC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sfAD_avg, sfAD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sfBC_avg, sfBC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sfBD_avg, sfBD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sfCD_avg, sfCD_std)
			
			res += "{0:.5f}\t{1:.5f}\t".format(sxA_avg, sxA_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sxB_avg, sxB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sxC_avg, sxC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(sxD_avg, sxD_std)
		
			res += "{0:.5f}\t{1:.5f}\t".format(ssAB_avg, ssAB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(ssAC_avg, ssAC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(ssAD_avg, ssAD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(ssBC_avg, ssBC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(ssBD_avg, ssBD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(ssCD_avg, ssCD_std)
			#res += "{0:.5f}\t{1:.5f}\t".format(successive_ss_avg, successive_ss_std)

			res += "{0:.5f}\t{1:.5f}\t".format(piA_avg, piA_std)
			res += "{0:.5f}\t{1:.5f}\t".format(piB_avg, piB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(piC_avg, piC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(piD_avg, piD_std)

			res += "{0:.5f}\t".format(pearson_r_pi_AB)
			res += "{0:.5f}\t".format(pearson_r_pi_AC)
			res += "{0:.5f}\t".format(pearson_r_pi_AD)
			res += "{0:.5f}\t".format(pearson_r_pi_BC)
			res += "{0:.5f}\t".format(pearson_r_pi_BD)
			res += "{0:.5f}\t".format(pearson_r_pi_CD)

			res += "{0:.5f}\t{1:.5f}\t".format(thetaA_avg, thetaA_std)
			res += "{0:.5f}\t{1:.5f}\t".format(thetaB_avg, thetaB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(thetaC_avg, thetaC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(thetaD_avg, thetaD_std)

			res += "{0:.5f}\t".format(pearson_r_theta_AB)
			res += "{0:.5f}\t".format(pearson_r_theta_AC)
			res += "{0:.5f}\t".format(pearson_r_theta_AD)
			res += "{0:.5f}\t".format(pearson_r_theta_BC)
			res += "{0:.5f}\t".format(pearson_r_theta_BD)
			res += "{0:.5f}\t".format(pearson_r_theta_CD)

			res += "{0:.5f}\t{1:.5f}\t".format(DtajA_avg, DtajA_std)
			res += "{0:.5f}\t{1:.5f}\t".format(DtajB_avg, DtajB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(DtajC_avg, DtajC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(DtajD_avg, DtajD_std)

			res += "{0:.5f}\t{1:.5f}\t".format(divAB_avg, divAB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(divAC_avg, divAC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(divAD_avg, divAD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(divBC_avg, divBC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(divBD_avg, divBD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(divCD_avg, divCD_std)
			
			res += "{0:.5f}\t{1:.5f}\t".format(netdivAB_avg, netdivAB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(netdivAC_avg, netdivAC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(netdivAD_avg, netdivAD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(netdivBC_avg, netdivBC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(netdivBD_avg, netdivBD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(netdivCD_avg, netdivCD_std)

			res += "{0:.5f}\t{1:.5f}\t".format(minDivAB_avg, minDivAB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(minDivAC_avg, minDivAC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(minDivAD_avg, minDivAD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(minDivBC_avg, minDivBC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(minDivBD_avg, minDivBD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(minDivCD_avg, minDivCD_std)

			res += "{0:.5f}\t{1:.5f}\t".format(maxDivAB_avg, maxDivAB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(maxDivAC_avg, maxDivAC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(maxDivAD_avg, maxDivAD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(maxDivBC_avg, maxDivBC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(maxDivBD_avg, maxDivBD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(maxDivCD_avg, maxDivCD_std)

			res += "{0:.5f}\t{1:.5f}\t".format(GminAB_avg, GminAB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GminAC_avg, GminAC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GminAD_avg, GminAD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GminBC_avg, GminBC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GminBD_avg, GminBD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GminCD_avg, GminCD_std)

			res += "{0:.5f}\t{1:.5f}\t".format(GmaxAB_avg, GmaxAB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GmaxAC_avg, GmaxAC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GmaxAD_avg, GmaxAD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GmaxBC_avg, GmaxBC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GmaxBD_avg, GmaxBD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(GmaxCD_avg, GmaxCD_std)
			
			res += "{0:.5f}\t{1:.5f}\t".format(FST_AB_avg, FST_AB_std)
			res += "{0:.5f}\t{1:.5f}\t".format(FST_AC_avg, FST_AC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(FST_AD_avg, FST_AD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(FST_BC_avg, FST_BC_std)
			res += "{0:.5f}\t{1:.5f}\t".format(FST_BD_avg, FST_BD_std)
			res += "{0:.5f}\t{1:.5f}\t".format(FST_CD_avg, FST_CD_std)
			
			res += "{0:.5f}\t{1:.5f}\t".format(D_AB_C_avg, D_AB_C_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_AB_C.count(-1), D_AB_C.count(1))
			res += "{0:.5f}\t{1:.5f}\t".format(fd_AB_C_avg, fd_AB_C_std)
			res += "{0:.5f}\t{1:.5f}\t".format(fhom_AB_C_avg, fhom_AB_C_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_AB_D_avg, D_AB_D_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_AB_D.count(-1), D_AB_D.count(1))
			res += "{0:.5f}\t{1:.5f}\t".format(fd_AB_D_avg, fd_AB_D_std)
			res += "{0:.5f}\t{1:.5f}\t".format(fhom_AB_D_avg, fhom_AB_D_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_BA_C_avg, D_BA_C_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_BA_C.count(-1), D_BA_C.count(1))
			res += "{0:.5f}\t{1:.5f}\t".format(fd_BA_C_avg, fd_BA_C_std)
			res += "{0:.5f}\t{1:.5f}\t".format(fhom_BA_C_avg, fhom_BA_C_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_BA_D_avg, D_BA_D_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_BA_D.count(-1), D_BA_D.count(1))
			res += "{0:.5f}\t{1:.5f}\t".format(fd_BA_D_avg, fd_BA_D_std)
			res += "{0:.5f}\t{1:.5f}\t".format(fhom_BA_D_avg, fhom_BA_D_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_CD_A_avg, D_CD_A_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_CD_A.count(-1), D_CD_A.count(1))
			res += "{0:.5f}\t{1:.5f}\t".format(fd_CD_A_avg, fd_CD_A_std)
			res += "{0:.5f}\t{1:.5f}\t".format(fhom_CD_A_avg, fhom_CD_A_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_CD_B_avg, D_CD_B_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_CD_B.count(-1), D_CD_B.count(1))
			res += "{0:.5f}\t{1:.5f}\t".format(fd_CD_B_avg, fd_CD_B_std)
			res += "{0:.5f}\t{1:.5f}\t".format(fhom_CD_B_avg, fhom_CD_B_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_DC_A_avg, D_DC_A_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_DC_A.count(-1), D_DC_A.count(1))
			res += "{0:.5f}\t{1:.5f}\t".format(fd_DC_A_avg, fd_DC_A_std)
			res += "{0:.5f}\t{1:.5f}\t".format(fhom_DC_A_avg, fhom_DC_A_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_DC_B_avg, D_DC_B_std)
			res += "{0:.5f}\t{1:.5f}\t".format(D_DC_B.count(-1), D_DC_B.count(1))
			res += "{0:.5f}\t{1:.5f}\t".format(fd_DC_B_avg, fd_DC_B_std)
			res += "{0:.5f}\t{1:.5f}".format(fhom_DC_B_avg, fhom_DC_B_std)
			#res += "{0:.5f}\t".format(pearson_r_div_netDiv)
			#res += "{0:.5f}\t".format(pearson_r_div_FST)
			#res += "{0:.5f}\t".format(pearson_r_netDiv_FST)
			#res += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}".format(ss_sf, ss_noSf, noSs_sf, noSs_noSf) # total number of ss_sf, ss_noSf and noSs_sf noSs_noSf loci
			#res += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}".format(ss_sf/(1.0*nLoci), ss_noSf/(1.0*nLoci), noSs_sf/(1.0*nLoci), noSs_noSf/(1.0*nLoci)) # proportion of ss_sf, ss_noSf, etc ... loci
			res += "\n"
			outfile.write(res)
	infile.close()




## ms' output file
outfile = open("ABCstat_allLoci.txt", "w")
# header
res = "contig\tpos_first_SNP\tpos_last_SNP\tname\tlength\tbialsites_avg\tbialsites_std\t"

res += "sfAB_avg\tsfAB_std\t"
res += "sfAC_avg\tsfAC_std\t"
res += "sfAD_avg\tsfAD_std\t"
res += "sfBC_avg\tsfBC_std\t"
res += "sfBD_avg\tsfBD_std\t"
res += "sfCD_avg\tsfCD_std\t"

res += "sxA_avg\tsxA_std\t"
res += "sxB_avg\tsxB_std\t"
res += "sxC_avg\tsxC_std\t"
res += "sxD_avg\tsxD_std\t"

res += "ssAB_avg\tssAB_std\t"
res += "ssAC_avg\tssAC_std\t"
res += "ssAD_avg\tssAD_std\t"
res += "ssBC_avg\tssBC_std\t"
res += "ssBD_avg\tssBD_std\t"
res += "ssCD_avg\tssCD_std\t"

res += "piA_avg\tpiA_std\t"
res += "piB_avg\tpiB_std\t"
res += "piC_avg\tpiC_std\t"
res += "piD_avg\tpiD_std\t"

res += "pearson_r_pi_AB\t"
res += "pearson_r_pi_AC\t"
res += "pearson_r_pi_AD\t"
res += "pearson_r_pi_BC\t"
res += "pearson_r_pi_BD\t"
res += "pearson_r_pi_CD\t"

res += "thetaA_avg\tthetaA_std\t"
res += "thetaB_avg\tthetaB_std\t"
res += "thetaC_avg\tthetaC_std\t"
res += "thetaD_avg\tthetaD_std\t"

res += "pearson_r_theta_AB\t"
res += "pearson_r_theta_AC\t"
res += "pearson_r_theta_AD\t"
res += "pearson_r_theta_BC\t"
res += "pearson_r_theta_BD\t"
res += "pearson_r_theta_CD\t"

res += "DtajA_avg\tDtajA_std\t"
res += "DtajB_avg\tDtajB_std\t"
res += "DtajC_avg\tDtajC_std\t"
res += "DtajD_avg\tDtajD_std\t"

res += "divAB_avg\tdivAB_std\t"
res += "divAC_avg\tdivAC_std\t"
res += "divAD_avg\tdivAD_std\t"
res += "divBC_avg\tdivBC_std\t"
res += "divBD_avg\tdivBD_std\t"
res += "divCD_avg\tdivCD_std\t"

res += "netdivAB_avg\tnetdivAB_std\t"
res += "netdivAC_avg\tnetdivAC_std\t"
res += "netdivAD_avg\tnetdivAD_std\t"
res += "netdivBC_avg\tnetdivBC_std\t"
res += "netdivBD_avg\tnetdivBD_std\t"
res += "netdivCD_avg\tnetdivCD_std\t"

res += "minDivAB_avg\tminDivAB_std\t"
res += "minDivAC_avg\tminDivAC_std\t"
res += "minDivAD_avg\tminDivAD_std\t"
res += "minDivBC_avg\tminDivBC_std\t"
res += "minDivBD_avg\tminDivBD_std\t"
res += "minDivCD_avg\tminDivCD_std\t"

res += "maxDivAB_avg\tmaxDivAB_std\t"
res += "maxDivAC_avg\tmaxDivAC_std\t"
res += "maxDivAD_avg\tmaxDivAD_std\t"
res += "maxDivBC_avg\tmaxDivBC_std\t"
res += "maxDivBD_avg\tmaxDivBD_std\t"
res += "maxDivCD_avg\tmaxDivCD_std\t"

res += "GminAB_avg\tGminAB_std\t"
res += "GminAC_avg\tGminAC_std\t"
res += "GminAD_avg\tGminAD_std\t"
res += "GminBC_avg\tGminBC_std\t"
res += "GminBD_avg\tGminBD_std\t"
res += "GminCD_avg\tGminCD_std\t"

res += "GmaxAB_avg\tGmaxAB_std\t"
res += "GmaxAC_avg\tGmaxAC_std\t"
res += "GmaxAD_avg\tGmaxAD_std\t"
res += "GmaxBC_avg\tGmaxBC_std\t"
res += "GmaxBD_avg\tGmaxBD_std\t"
res += "GmaxCD_avg\tGmaxCD_std\t"

res += "FST_AB_avg\tFST_AB_std\t"
res += "FST_AC_avg\tFST_AC_std\t"
res += "FST_AD_avg\tFST_AD_std\t"
res += "FST_BC_avg\tFST_BC_std\t"
res += "FST_BD_avg\tFST_BD_std\t"
res += "FST_CD_avg\tFST_CD_std\t"

res += 'D_AB_C_avg\tD_AB_C_std\t'
res += 'D_negOne_AB_C\tD_posOne_AB_C\t'
res += 'fd_AB_C_avg\tfd_AB_C_std\t'
res += 'fhom_AB_C_avg\tfhom_AB_C_std\t'
res += 'D_AB_D_avg\tD_AB_D_std\t'
res += 'D_negOne_AB_D\tD_posOne_AB_D\t'
res += 'fd_AB_D_avg\tfd_AB_D_std\t'
res += 'fhom_AB_D_avg\tfhom_AB_D_std\t'
res += 'D_BA_C_avg\tD_BA_C_std\t'
res += 'D_negOne_BA_C\tD_posOne_BA_C\t'
res += 'fd_BA_C_avg\tfd_BA_C_std\t'
res += 'fhom_BA_C_avg\tfhom_BA_C_std\t'
res += 'D_BA_D_avg\tD_BA_D_std\t'
res += 'D_negOne_BA_D\tD_posOne_BA_D\t'
res += 'fd_BA_D_avg\tfd_BA_D_std\t'
res += 'fhom_BA_D_avg\tfhom_BA_D_std\t'
res += 'D_CD_A_avg\tD_CD_A_std\t'
res += 'D_negOne_CD_A\tD_posOne_CD_A\t'
res += 'fd_CD_A_avg\tfd_CD_A_std\t'
res += 'fhom_CD_A_avg\tfhom_CD_A_std\t'
res += 'D_CD_B_avg\tD_CD_B_std\t'
res += 'D_negOne_CD_B\tD_posOne_CD_B\t'
res += 'fd_CD_B_avg\tfd_CD_B_std\t'
res += 'fhom_CD_B_avg\tfhom_CD_B_std\t'
res += 'D_DC_A_avg\tD_DC_A_std\t'
res += 'D_negOne_DC_A\tD_posOne_DC_A\t'
res += 'fd_DC_A_avg\tfd_DC_A_std\t'
res += 'fhom_DC_A_avg\tfhom_DC_A_std\t'
res += 'D_DC_B_avg\tD_DC_B_std\t'
res += 'D_negOne_DC_B\tD_posOne_DC_B\t'
res += 'fd_DC_B_avg\tfd_DC_B_std\t'
res += 'fhom_DC_B_avg\tfhom_DC_B_std\n'
outfile.write(res)


locus_list = [ i.split('.')[0] for i in os.listdir('./') if '.ms' in i and 'data' not in i ]
for locus in locus_list:
        computeStats(locus, outfile)
outfile.close()


