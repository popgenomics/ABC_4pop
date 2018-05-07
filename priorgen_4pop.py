#!/home/roux/python/Python-2.7.14/python
# #!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from numpy.random import seed

#seed(100)

help = "\t\033[1;31;40mTakes one model specifier, one indicator of migration and a number of multilocus simulations as arguments:\033[0m\n\t"
help += "Accepted model specifiers:\n\t\t"
help += "\n\t\t".join(["SI_1N", "SI_2N", "SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N"])
help += "\n\n"
help += "\tAccepted indicators of migration:\n\t\t"
help += "\n\t\t".join(["none (for SI)", "A (for SC)", "B (for SC)", "C (for SC)", "D (for SC)", "AC (for SC)", "BD (for SC)", "ACBD (for SC)"])
help += "\n\n"
help += "\t\033[1;32;40m#SI\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\033[0m\n\n" # no migration
help += "\t\033[1;32;40m#SC\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\033[0m\n\n" # secondary contact between pop_A and pop_C
help += "\t\033[1;32;40m#1M\033[0m\t\tAll loci share the same M13, M31, M24 and M42\033[0m\n"
help += "\t\033[1;32;40m#1N\033[0m\t\tAll loci share the same N1, N2, N3, N4, Na_12, Na_34 and Na\033[0m\n"
help += "\t\033[1;32;40m#2M\033[0m\t\tAll loci have different values for M13, M31, M24 and M42 (Beta distributed)\033[0m\n"
help += "\t\033[1;32;40m#2N\033[0m\t\tAll loci have different N1, N2, N3, N4, Na_12, Na_34 and Na (Beta distributed)\033[0m\n\n"
help += "\t\033[1;32;40m#A\033[0m\t\tSecondary contact between pop_A and pop_C, unidirectional A<--C\033[0m\n"
help += "\t\033[1;32;40m#B\033[0m\t\tSecondary contact between pop_B and pop_D, unidirectional B<--D\033[0m\n"
help += "\t\033[1;32;40m#C\033[0m\t\tSecondary contact between pop_A and pop_C, unidirectional A-->C\033[0m\n"
help += "\t\033[1;32;40m#D\033[0m\t\tSecondary contact between pop_B and pop_D, unidirectional B-->D\033[0m\n"
help += "\t\033[1;32;40m#AC\033[0m\t\tSecondary contact between pop_A and pop_C, bidirectional A<->C\033[0m\n"
help += "\t\033[1;32;40m#BD\033[0m\t\tSecondary contact between pop_B and pop_D, bidirectional B<->D\033[0m\n"
help += "\t\033[1;32;40m#ACBD\033[0m\t\tSecondary contact between A<->C and B<->D\033[0m\n\n"
help += "\t\033[1;32;40mExample: ./priorgen_4pop_beta.py SC_2M_2N AC 1000\033[0m\n"

if len(sys.argv) != 4:
	print(help)
	sys.exit()

# Configuration of the prior distribution
nMultilocus = int(sys.argv[3])
N_bound = [0, 80]
T_bound = [0, 10]
M_bound = [0, 40]
shape_bound = [0.01, 50]

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
N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
N3 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
N4 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)

Na_12 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) # ancestor between N1 and N2
Na_34 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) # ancestor between N3 and N4
Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) # N1N2-N3N4 
	
## factor of local reduction in Ne. Model of "background selection"
shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

## Miration rates
if sys.argv[2] == 'A':
	M13 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M31 = [0]*nMultilocus
	M24 = [0]*nMultilocus
	M42 = [0]*nMultilocus

if sys.argv[2] == 'B':
	M13 = [0]*nMultilocus
	M31 = [0]*nMultilocus
	M24 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M42 = [0]*nMultilocus

if sys.argv[2] == 'C':
	M13 = [0]*nMultilocus
	M31 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M24 = [0]*nMultilocus
	M42 = [0]*nMultilocus

if sys.argv[2] == 'D':
	M13 = [0]*nMultilocus
	M31 = [0]*nMultilocus
	M24 = [0]*nMultilocus
	M42 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

if sys.argv[2] == 'AC':
	M13 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M31 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M24 = [0]*nMultilocus
	M42 = [0]*nMultilocus

if sys.argv[2] == 'BD':
	M13 = [0]*nMultilocus
	M31 = [0]*nMultilocus
	M24 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M42 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

if sys.argv[2] == 'ACBD':
	M13 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M31 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M24 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M42 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

## factor of local reduction in Me. One Beta distribution for each of the migration rates 
shape_M13_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
shape_M13_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
shape_M31_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
shape_M31_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

shape_M24_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
shape_M24_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
shape_M42_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
shape_M42_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	
## times
Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
Tsplit_12 = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]
Tsplit_34 = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]
Tsc_13 = [ uniform(low = T_bound[0], high = min([Tsplit_12[i], Tsplit_34[i]]), size = 1)[0] for i in range(nMultilocus) ]
Tsc_24 = [ uniform(low = T_bound[0], high = min([Tsplit_12[i], Tsplit_34[i]]), size = 1)[0] for i in range(nMultilocus) ]


# Four pop version
if sys.argv[1] == "SI_1N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\n" # no migration
	# param monolocus: values that will be read by ms
	priorfile = "N_popA\tN_popB\tN_popC\tN_popD\tNa_AB\tNa_CD\tNa\tTsplit_AB\tTsplit_CD\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1:.5f}\t{2:.5f}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], nsamC[locus], nsamD[locus], N1[sim], N2[sim], N3[sim], N4[sim], 0, 0, 0, 0, Tsc_13[sim], Tsc_13[sim], Tsc_24[sim], Tsc_24[sim], Tsplit_12[sim], Tsplit_12[sim], Na_12[sim], Tsplit_34[sim], Tsplit_34[sim], Na_34[sim], Tsplit[sim], Tsplit[sim], Na[sim]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SI_2N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\n" # no migration
	# param monolocus: values that will be read by ms
	priorfile = "N_popA\tN_popB\tN_popC\tN_popD\tNa_AB\tNa_CD\tNa\tshape_a_N\tshape_b_N\tTsplit_AB\tTsplit_CD\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		N1_vec = [ N1[sim]*i for i in scalar_N ]
		N2_vec = [ N2[sim]*i for i in scalar_N ]
		N3_vec = [ N3[sim]*i for i in scalar_N ]
		N4_vec = [ N4[sim]*i for i in scalar_N ]
		Na_12_vec = [ Na_12[sim]*i for i in scalar_N ]
		Na_34_vec = [ Na_34[sim]*i for i in scalar_N ]
		Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1:.5f}\t{2:.5f}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], nsamC[locus], nsamD[locus], N1_vec[locus], N2_vec[locus], N3_vec[locus], N4_vec[locus], 0, 0, 0, 0, Tsc_13[sim], Tsc_13[sim], Tsc_24[sim], Tsc_24[sim], Tsplit_12[sim], Tsplit_12[sim], Na_12_vec[locus], Tsplit_34[sim], Tsplit_34[sim], Na_34_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_1M_1N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\n" # no migration
	# param monolocus: values that will be read by ms
	priorfile = "N_popA\tN_popB\tN_popC\tN_popD\tNa_AB\tNa_CD\tNa\tTsplit_AB\tTsplit_CD\tTsplit\tT_SC_AC\tT_SC_BD\tM13\tM31\tM24\tM42\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim], Tsc_13[sim], Tsc_24[sim], M13[sim], M31[sim], M24[sim], M42[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1:.5f}\t{2:.5f}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], nsamC[locus], nsamD[locus], N1[sim], N2[sim], N3[sim], N4[sim], M13[sim], M31[sim], M24[sim], M42[sim], Tsc_13[sim], Tsc_13[sim], Tsc_24[sim], Tsc_24[sim], Tsplit_12[sim], Tsplit_12[sim], Na_12[sim], Tsplit_34[sim], Tsplit_34[sim], Na_34[sim], Tsplit[sim], Tsplit[sim], Na[sim]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SC_1M_2N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\n" # no migration
	# param monolocus: values that will be read by ms
	priorfile = "N_popA\tN_popB\tN_popC\tN_popD\tNa_AB\tNa_CD\tNa\tshape_a_N\tshape_b_N\tTsplit_AB\tTsplit_CD\tTsplit\tT_SC_AC\tT_SC_BD\tM13\tM31\tM24\tM42\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim], Tsc_13[sim], Tsc_24[sim], M13[sim], M31[sim], M24[sim], M42[sim])

		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		N1_vec = [ N1[sim]*i for i in scalar_N ]
		N2_vec = [ N2[sim]*i for i in scalar_N ]
		N3_vec = [ N3[sim]*i for i in scalar_N ]
		N4_vec = [ N4[sim]*i for i in scalar_N ]
		Na_12_vec = [ Na_12[sim]*i for i in scalar_N ]
		Na_34_vec = [ Na_34[sim]*i for i in scalar_N ]
		Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1:.5f}\t{2:.5f}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], nsamC[locus], nsamD[locus], N1_vec[locus], N2_vec[locus], N3_vec[locus], N4_vec[locus], M13[sim], M31[sim], M24[sim], M42[sim], Tsc_13[sim], Tsc_13[sim], Tsc_24[sim], Tsc_24[sim], Tsplit_12[sim], Tsplit_12[sim], Na_12_vec[locus], Tsplit_34[sim], Tsplit_34[sim], Na_34_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_2M_1N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\n" # no migration
	# param monolocus: values that will be read by ms
	priorfile = "N_popA\tN_popB\tN_popC\tN_popD\tNa_AB\tNa_CD\tNa\tTsplit_AB\tTsplit_CD\tTsplit\tT_SC_AC\tT_SC_BD\tM13\tshape_M13_a\tshape_M13_b\tM31\tshape_M31_a\tshape_M31_b\tM24\tshape_M24_a\tshape_M24_b\tM42\tshape_M42_a\tshape_M42_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim], Tsc_13[sim], Tsc_24[sim], M13[sim], shape_M13_a[sim], shape_M13_b[sim], M31[sim], shape_M31_a[sim], shape_M31_b[sim], M24[sim], shape_M24_a[sim], shape_M24_b[sim], M42[sim], shape_M42_a[sim], shape_M42_b[sim])
		
		# vectors of size 'nLoci' containing parameters
		scalar_M13 = beta(shape_M13_a[sim], shape_M13_b[sim], size = nLoci)
		scalar_M31 = beta(shape_M31_a[sim], shape_M31_b[sim], size = nLoci)
		scalar_M24 = beta(shape_M24_a[sim], shape_M24_b[sim], size = nLoci)
		scalar_M42 = beta(shape_M42_a[sim], shape_M42_b[sim], size = nLoci)
		M13_vec = [ M13[sim] * i for i in scalar_M13 ]
		M31_vec = [ M31[sim] * i for i in scalar_M31 ]
		M24_vec = [ M24[sim] * i for i in scalar_M24 ]
		M42_vec = [ M42[sim] * i for i in scalar_M42 ]
		
		for locus in range(nLoci):
			print("{0}\t{1:.5f}\t{2:.5f}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], nsamC[locus], nsamD[locus], N1[sim], N2[sim], N3[sim], N4[sim], M13_vec[locus], M31_vec[locus], M24_vec[locus], M42_vec[locus], Tsc_13[sim], Tsc_13[sim], Tsc_24[sim], Tsc_24[sim], Tsplit_12[sim], Tsplit_12[sim], Na_12[sim], Tsplit_34[sim], Tsplit_34[sim], Na_34[sim], Tsplit[sim], Tsplit[sim], Na[sim]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_2M_2N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\n" # no migration
	# param monolocus: values that will be read by ms
	priorfile = "N_popA\tN_popB\tN_popC\tN_popD\tNa_AB\tNa_CD\tNa\tshape_N_a\tshape_N_b\tTsplit_AB\tTsplit_CD\tTsplit\tT_SC_AC\tT_SC_BD\tM13\tshape_M13_a\tshape_M13_b\tM31\tshape_M31_a\tshape_M31_b\tM24\tshape_M24_a\tshape_M24_b\tM42\tshape_M42_a\tshape_M42_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim], Tsc_13[sim], Tsc_24[sim], M13[sim], shape_M13_a[sim], shape_M13_b[sim], M31[sim], shape_M31_a[sim], shape_M31_b[sim], M24[sim], shape_M24_a[sim], shape_M24_b[sim], M42[sim], shape_M42_a[sim], shape_M42_b[sim])
		
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		N1_vec = [ N1[sim]*i for i in scalar_N ]
		N2_vec = [ N2[sim]*i for i in scalar_N ]
		N3_vec = [ N3[sim]*i for i in scalar_N ]
		N4_vec = [ N4[sim]*i for i in scalar_N ]
		Na_12_vec = [ Na_12[sim]*i for i in scalar_N ]
		Na_34_vec = [ Na_34[sim]*i for i in scalar_N ]
		Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		scalar_M13 = beta(shape_M13_a[sim], shape_M13_b[sim], size = nLoci)
		scalar_M31 = beta(shape_M31_a[sim], shape_M31_b[sim], size = nLoci)
		scalar_M24 = beta(shape_M24_a[sim], shape_M24_b[sim], size = nLoci)
		scalar_M42 = beta(shape_M42_a[sim], shape_M42_b[sim], size = nLoci)
		M13_vec = [ M13[sim] * i for i in scalar_M13 ]
		M31_vec = [ M31[sim] * i for i in scalar_M31 ]
		M24_vec = [ M24[sim] * i for i in scalar_M24 ]
		M42_vec = [ M42[sim] * i for i in scalar_M42 ]
		
		for locus in range(nLoci):
			print("{0}\t{1:.5f}\t{2:.5f}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], nsamC[locus], nsamD[locus], N1_vec[locus], N2_vec[locus], N3_vec[locus], N4_vec[locus], M13_vec[locus], M31_vec[locus], M24_vec[locus], M42_vec[locus], Tsc_13[sim], Tsc_13[sim], Tsc_24[sim], Tsc_24[sim], Tsplit_12[sim], Tsplit_12[sim], Na_12_vec[locus], Tsplit_34[sim], Tsplit_34[sim], Na_34_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



