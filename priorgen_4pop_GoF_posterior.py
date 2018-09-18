#!/usr/bin/python
# #!/home/roux/python/Python-2.7.14/python
# -*- coding: utf-8 -*-

import sys
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from numpy.random import seed

#seed(100)

help = "\n\t\033[1;31;40mTakes one indicator of migration, one model specifier, a posterior file and a number of replicates as arguments:\033[0m\n\t"
help += "\033[1;32;40m#Command line for ms: \033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 1 2 0 -em tbs 2 1 0 -ej tbs 2 1 -en tbs 1 tbs -em tbs 3 4 0 -em tbs 4 3 0 -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\033[0m\n\n"
help += "\t\033[1;31;40mAccepted model specifiers:\033[0m\n\t"
help += ", ".join(["1M_1N", "1M_2N", "2M_1N", "2M_2N"])
help += "\n"
help += "\t\033[1;32;40m#1M\033[0m\t\tAll loci share the same M13, M31, M24 and M42\033[0m\n"
help += "\t\033[1;32;40m#1N\033[0m\t\tAll loci share the same N1, N2, N3, N4, Na_12, Na_34 and Na\033[0m\n"
help += "\t\033[1;32;40m#2M\033[0m\t\tAll loci have different values for M13, M31, M24 and M42 (Beta distributed)\033[0m\n"
help += "\t\033[1;32;40m#2N\033[0m\t\tAll loci have different N1, N2, N3, N4, Na_12, Na_34 and Na (Beta distributed)\033[0m\n"
help += "\n"
help += "\t\033[1;31;40mAccepted indicators of migration:\033[0m\n\t"
help += ", ".join(["AB AC BD CD"])
help += "\n\t\033[1;32;40m#ac\033[0m\t\tNo migration between pop_A and pop_C, A | C\033[0m\n"
help += "\t\033[1;32;40m#Ac\033[0m\t\tSecondary contact between pop_A and pop_C, unidirectional A<--C\033[0m\n"
help += "\t\033[1;32;40m#aC\033[0m\t\tSecondary contact between pop_A and pop_C, unidirectional A-->C\033[0m\n"
help += "\t\033[1;32;40m#AC\033[0m\t\tSecondary contact between pop_A and pop_C, bidirectional A<->C\033[0m\n\n"
help += "\t\033[1;32;40mExample: ./ABC_4pop_GoF_posterior.py 2M_2N AB_AC_BD_CD post_AB_AC_BD_CD_2datasets_250_TajD_1 4\033[0m\n"

if len(sys.argv) != 8:
	print(help)
	sys.exit()

# Configuration of the prior distribution
posterior_file = sys.argv[6]
nReplicates = int(sys.argv[7])

# read the posterior file
# values for parameters are stored in the dic param[parameter][line]
infile = open(posterior_file, "r")
header = infile.readline().strip()
header = header.replace('"', '')
header = header.replace(" ", "\t")
header = header.split("\t")

param = {}
colonne = {}
cnt = -1
for i in header:
	cnt += 1
	colonne[cnt] = i
	param[i] = []

nMultilocus = 0
for line in infile:
	nMultilocus += 1
	line = line.strip().replace(" ", "\t")
	line = line.split("\t")
	for i in range(len(line)):
		param[colonne[i]].append(float(line[i]))
infile.close()


# read bpfile
infile = open("bpfile", "r")
tmp = infile.readline()
L = [ float(i) for i in infile.readline().strip().split("\t") ]
nsamA = [ int(i) for i in infile.readline().strip().split("\t") ]
nsamB = [ int(i) for i in infile.readline().strip().split("\t") ]
nsamC = [ int(i) for i in infile.readline().strip().split("\t") ]
nsamD = [ int(i) for i in infile.readline().strip().split("\t") ]
theta = [ float(i) for i in infile.readline().strip().split("\t") ]
rho = [ float(i) for i in infile.readline().strip().split("\t") ]
infile.close()

# number of loci
nLoci = len(L)

# sum of nsamA + nsamB
nsam_tot = [ nsamA[i] + nsamB[i] + nsamC[i] + nsamD[i] for i in range(nLoci) ]

# param multilocus: values that will be printed in priorfile.txt
## N = N_pop_i / Nref
# N1      N2      N3      N4      Na_12   Na_34   Na      shape_N_a       shape_N_b       Tsplit_12       Tsplit_34       Tsplit  Tsc_13  Tsc_24  M12  M21      M13     shape_M13_a     shape_M13_b     M31     shape_M31_a     shape_M31_b     M14     shape_M14_a     shape_M14_b     M41     shape_M41_a  shape_M41_b      M23     shape_M23_a     shape_M23_b     M32     shape_M32_a     shape_M32_b     M24     shape_M24_a      shape_M24_b    M42     shape_M42_a   shape_M42_b     M34     M43
N1 = param['N1']
N2 = param['N2']
N3 = param['N3']
N4 = param['N4']

Na_12 = param['Na_12']
Na_34 = param['Na_34']
Na = param['Na']
	
## factor of local reduction in Ne. Model of "background selection"
shape_N_a = param['shape_N_a']
shape_N_b = param['shape_N_b']

## Miration rates
### between A and B
if 'a' in sys.argv[2]:
	M12 = [0] * nMultilocus
if 'A' in sys.argv[2]:
	M12 = param['M12']
if 'b' in sys.argv[2]:
	M21 = [0] * nMultilocus
if 'B' in sys.argv[2]:
	M21 = param['M21']

### between A and C
if 'a' in sys.argv[3]:
	M13 = [0] * nMultilocus
if 'A' in sys.argv[3]:
	M13 = param['M13']
if 'c' in sys.argv[3]:
	M31 = [0] * nMultilocus
if 'C' in sys.argv[3]:
	M31 = param['M31']
if 'a' in sys.argv[3] and 'c' in sys.argv[3]:
	Tsc_13 = [0] * nMultilocus

### between B and D
if 'b' in sys.argv[4]:
	M24 = [0] * nMultilocus
if 'B' in sys.argv[4]:
	M24 = param['M24']
if 'd' in sys.argv[4]:
	M42 = [0] * nMultilocus
if 'D' in sys.argv[4]:
	M42 = param['M42']
if 'b' in sys.argv[4] and 'd' in sys.argv[4]:
	Tsc_24 = [0] * nMultilocus

### between C and D
if 'c' in sys.argv[5]:
	M34 = [0] * nMultilocus
if 'C' in sys.argv[5]:
	M34 = param['M34']
if 'd' in sys.argv[5]:
	M43 = [0] * nMultilocus
if 'D' in sys.argv[5]:
	M43 = param['M43']


## factor of local reduction in Me. One Beta distribution for each of the migration rates 
try:
	shape_M13_a = param['shape_M13_a']
except:
	shape_M13_a = [999]*nMultilocus
	sys.stderr.write("no shape_M13_a\n")

try:
	shape_M13_b = param['shape_M13_b']
except:
	shape_M13_b = [999]*nMultilocus
	sys.stderr.write("no shape_M13_b\n")

try:
	shape_M31_a = param['shape_M31_a']
except:
	shape_M31_a = [999]*nMultilocus
	sys.stderr.write("no shape_M31_a\n")

try:
	shape_M31_b = param['shape_M31_b']
except:
	shape_M31_b = [999]*nMultilocus
	sys.stderr.write("no shape_M31_b\n")

try:
	shape_M24_a = param['shape_M24_a']
except:
	shape_M24_a = [999]*nMultilocus
	sys.stderr.write("no shape_M24_a\n")

try:
	shape_M24_b = param['shape_M24_b']
except:
	shape_M24_b = [999]*nMultilocus
	sys.stderr.write("no shape_M24_b\n")

try:
	shape_M42_a = param['shape_M42_a']
except:
	shape_M42_a = [999]*nMultilocus
	sys.stderr.write("no shape_M42_a\n")

try:
	shape_M42_b = param['shape_M42_b']
except:
	shape_M42_b = [999]*nMultilocus
	sys.stderr.write("no shape_M42_b\n")

## times
Tsplit = param['Tsplit']
Tsplit_12 = param['Tsplit_12']
Tsplit_34 = param['Tsplit_34']

for i in range(len(Tsplit)):
	if Tsplit_34[i]>Tsplit[i]:
		Tsplit_34[i]=0.99999*Tsplit[i]
	if Tsplit_12[i]>Tsplit[i]:
		Tsplit_12[i]=0.99999*Tsplit[i]

try:
	Tsc_13 = param['Tsc_13']
except:
	sys.stderr.write("no Tsc_13\n")
try:
	Tsc_24 = param['Tsc_24'] 
except:
	sys.stderr.write("no Tsc_24\n")

# Four pop version
# param monolocus: values that will be read by ms
priorfile = "N1\tN2\tN3\tN4\tNa_12\tNa_34\tNa\tshape_N_a\tshape_N_b\tTsplit_12\tTsplit_34\tTsplit\tTsc_13\tTsc_24\tM12\tM21\tM13\tshape_M13_a\tshape_M13_b\tM31\tshape_M31_a\tshape_M31_b\tM14\tshape_M14_a\tshape_M14_b\tM41\tshape_M41_a\tshape_M41_b\tM23\tshape_M23_a\tshape_M23_b\tM32\tshape_M32_a\tshape_M32_b\tM24\tshape_M24_a\t shape_M24_b\tM42\tshape_M42_a\tshape_M42_b\tM34\tM43\n"

for sim in range(nMultilocus):
	#priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}\t{29:.5f}\t{30:.5f}\t{31:.5f}\t{32:.5f}\t{33:.5f}\t{34:.5f}\t{35:.5f}\t{36:.5f}\t{37:.5f}\t{38:.5f}\t{39:.5f}\t{40:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim], Tsc_13[sim], Tsc_24[sim], M12[sim], M21[sim], M13[sim], shape_M13_a[sim], shape_M13_b[sim], M31[sim], shape_M31_a[sim], shape_M31_b[sim], M14[sim], shape_M14_a[sim], shape_M14_b[sim], M41[sim], shape_M41_a[sim], shape_M41_b[sim], M23[sim], shape_M23_a[sim], shape_M23_b[sim], M32[sim], shape_M32_a[sim], shape_M32_b[sim], M24[sim], shape_M24_a[sim], shape_M24_b[sim], M42[sim], shape_M42_a[sim], shape_M42_b[sim], M34[sim], M43[sim])
	priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}\t{29:.5f}\t{30:.5f}\t{31:.5f}\t{32:.5f}\t{33:.5f}\t{34:.5f}\t{35:.5f}\t{36:.5f}\t{37:.5f}\t{38:.5f}\t{39:.5f}\t{40:.5f}\t{41:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim], Tsc_13[sim], Tsc_24[sim], M12[sim], M21[sim], M13[sim], shape_M13_a[sim], shape_M13_b[sim], M31[sim], shape_M31_a[sim], shape_M31_b[sim], 0, -9, -9, 0, -9, -9, 0, -9, -9, 0, -9, -9, M24[sim], shape_M24_a[sim], shape_M24_b[sim], M42[sim], shape_M42_a[sim], shape_M42_b[sim], M34[sim], M43[sim])
	# vectors of size 'nLoci' containing parameters
	## effective sizes
	if "1N" in sys.argv[1]:
		N1_vec = [N1[sim]] * nLoci
		N2_vec = [N2[sim]] * nLoci
		N3_vec = [N3[sim]] * nLoci
		N4_vec = [N4[sim]] * nLoci
		Na_12_vec = [Na_12[sim]] * nLoci
		Na_34_vec = [Na_34[sim]] * nLoci
		Na_vec = [Na[sim]] * nLoci
	
	if "2N" in sys.argv[1]:
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		N1_vec = [ N1[sim]*i for i in scalar_N ]
		N2_vec = [ N2[sim]*i for i in scalar_N ]
		N3_vec = [ N3[sim]*i for i in scalar_N ]
		N4_vec = [ N4[sim]*i for i in scalar_N ]
		Na_12_vec = [ Na_12[sim]*i for i in scalar_N ]
		Na_34_vec = [ Na_34[sim]*i for i in scalar_N ]
		Na_vec = [ Na[sim]*i for i in scalar_N ]
	
	## migration rates
	M12_vec = [M12[sim]] * nLoci
	M21_vec = [M21[sim]] * nLoci
	M34_vec = [M34[sim]] * nLoci
	M43_vec = [M43[sim]] * nLoci
	
	if "1M" in sys.argv[1]:
		M13_vec = [M13[sim]] * nLoci
		M31_vec = [M31[sim]] * nLoci
		M14_vec = [M14[sim]] * nLoci
		M41_vec = [M41[sim]] * nLoci
		M23_vec = [M23[sim]] * nLoci
		M32_vec = [M32[sim]] * nLoci
		M24_vec = [M24[sim]] * nLoci
		M42_vec = [M42[sim]] * nLoci
	
	if "2M" in sys.argv[1]:	
		scalar_M13 = beta(shape_M13_a[sim], shape_M13_b[sim], size = nLoci)
		M13_vec = [ M13[sim] * i for i in scalar_M13 ]
		scalar_M31 = beta(shape_M31_a[sim], shape_M31_b[sim], size = nLoci)
		M31_vec = [ M31[sim] * i for i in scalar_M31 ]
		
		scalar_M24 = beta(shape_M24_a[sim], shape_M24_b[sim], size = nLoci)
		M24_vec = [ M24[sim] * i for i in scalar_M24 ]
		scalar_M42 = beta(shape_M42_a[sim], shape_M42_b[sim], size = nLoci)
		M42_vec = [ M42[sim] * i for i in scalar_M42 ]

	for replicates in range(nReplicates):	
		for locus in range(nLoci):
			# msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 1 2 0 -em tbs 2 1 0 -ej tbs 2 1 -en tbs 1 tbs -em tbs 3 4 0 -em tbs 4 3 0 -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs
			print("{0}\t{1:.5f}\t{2:.5f}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}\t{29:.5f}\t{30:.5f}\t{31:.5f}\t{32:.5f}\t{33:.5f}\t{34:.5f}\t{35:.5f}\t{36:.5f}\t{37:.5f}\t{38:.5f}\t{39:.5f}\t{40:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], nsamC[locus], nsamD[locus], N1_vec[locus], N2_vec[locus], N3_vec[locus], N4_vec[locus], M12_vec[locus], M21_vec[locus], M13_vec[locus], M31_vec[locus], 0, 0, 0, 0, M24_vec[locus], M42_vec[locus], M34_vec[locus], M43_vec[locus], Tsc_13[sim], Tsc_13[sim], Tsc_24[sim], Tsc_24[sim], Tsplit_12[sim], Tsplit_12[sim], Tsplit_12[sim], Tsplit_12[sim], Na_12_vec[locus], Tsplit_34[sim], Tsplit_34[sim], Tsplit_34[sim], Tsplit_34[sim], Na_34_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

outfile = open("priorfile.txt", "w")
outfile.write(priorfile)
outfile.close()

# 1 nsam_tot[locus]
# 2 theta[locus]
# 3 rho[locus]
# 4 L[locus]
# 5 nsamA[locus]
# 6 nsamB[locus]
# 7 nsamC[locus]
# 8 nsamD[locus]
# 9 N1_vec[locus]
# 10 N2_vec[locus]
# 11 N3_vec[locus]
# 12 N4_vec[locus]
# 13 M12_vec[locus]
# 14 M21_vec[locus]
# 15 M13_vec[locus]
# 16 M31_vec[locus]
# 17 0
# 18 0
# 19 0
# 20 0
# 21 M24_vec[locus]
# 22 M42_vec[locus]
# 23 M34_vec[locus]
# 24 M43_vec[locus]
# 25 Tsc_13[sim]
# 26 Tsc_13[sim]
# 27 Tsc_24[sim]
# 28 Tsc_24[sim]
# 29 Tsplit_12[sim]
# 30 Tsplit_12[sim]
# 31 Tsplit_12[sim]
# 32 Tsplit_12[sim]
# 33 Na_12_vec[locus]
# 34 Tsplit_34[sim]
# 35 Tsplit_34[sim]
# 36 Tsplit_34[sim]
# 37 Tsplit_34[sim]
# 38 Na_34_vec[locus]
# 39 Tsplit[sim]
# 40 Tsplit[sim]
# 41 Na_vec[locus]))

