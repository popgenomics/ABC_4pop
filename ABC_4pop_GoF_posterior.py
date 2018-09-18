#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys

accepted_models = ["1M_1N", "1M_2N", "2M_1N", "2M_2N"]
accepted_migrations = ['ab_ac_bd_cd', 'ab_ac_bd_cD', 'ab_ac_bd_Cd', 'ab_ac_bd_CD', 'ab_ac_bD_cd', 'ab_ac_bD_cD', 'ab_ac_bD_Cd', 'ab_ac_bD_CD', 'ab_ac_Bd_cd', 'ab_ac_Bd_cD', 'ab_ac_Bd_Cd', 'ab_ac_Bd_CD', 'ab_ac_BD_cd', 'ab_ac_BD_cD', 'ab_ac_BD_Cd', 'ab_ac_BD_CD', 'ab_aC_bd_cd', 'ab_aC_bd_cD', 'ab_aC_bd_Cd', 'ab_aC_bd_CD', 'ab_aC_bD_cd', 'ab_aC_bD_cD', 'ab_aC_bD_Cd', 'ab_aC_bD_CD', 'ab_aC_Bd_cd', 'ab_aC_Bd_cD', 'ab_aC_Bd_Cd', 'ab_aC_Bd_CD', 'ab_aC_BD_cd', 'ab_aC_BD_cD', 'ab_aC_BD_Cd', 'ab_aC_BD_CD', 'ab_Ac_bd_cd', 'ab_Ac_bd_cD', 'ab_Ac_bd_Cd', 'ab_Ac_bd_CD', 'ab_Ac_bD_cd', 'ab_Ac_bD_cD', 'ab_Ac_bD_Cd', 'ab_Ac_bD_CD', 'ab_Ac_Bd_cd', 'ab_Ac_Bd_cD', 'ab_Ac_Bd_Cd', 'ab_Ac_Bd_CD', 'ab_Ac_BD_cd', 'ab_Ac_BD_cD', 'ab_Ac_BD_Cd', 'ab_Ac_BD_CD', 'ab_AC_bd_cd', 'ab_AC_bd_cD', 'ab_AC_bd_Cd', 'ab_AC_bd_CD', 'ab_AC_bD_cd', 'ab_AC_bD_cD', 'ab_AC_bD_Cd', 'ab_AC_bD_CD', 'ab_AC_Bd_cd', 'ab_AC_Bd_cD', 'ab_AC_Bd_Cd', 'ab_AC_Bd_CD', 'ab_AC_BD_cd', 'ab_AC_BD_cD', 'ab_AC_BD_Cd', 'ab_AC_BD_CD', 'aB_ac_bd_cd', 'aB_ac_bd_cD', 'aB_ac_bd_Cd', 'aB_ac_bd_CD', 'aB_ac_bD_cd', 'aB_ac_bD_cD', 'aB_ac_bD_Cd', 'aB_ac_bD_CD', 'aB_ac_Bd_cd', 'aB_ac_Bd_cD', 'aB_ac_Bd_Cd', 'aB_ac_Bd_CD', 'aB_ac_BD_cd', 'aB_ac_BD_cD', 'aB_ac_BD_Cd', 'aB_ac_BD_CD', 'aB_aC_bd_cd', 'aB_aC_bd_cD', 'aB_aC_bd_Cd', 'aB_aC_bd_CD', 'aB_aC_bD_cd', 'aB_aC_bD_cD', 'aB_aC_bD_Cd', 'aB_aC_bD_CD', 'aB_aC_Bd_cd', 'aB_aC_Bd_cD', 'aB_aC_Bd_Cd', 'aB_aC_Bd_CD', 'aB_aC_BD_cd', 'aB_aC_BD_cD', 'aB_aC_BD_Cd', 'aB_aC_BD_CD', 'aB_Ac_bd_cd', 'aB_Ac_bd_cD', 'aB_Ac_bd_Cd', 'aB_Ac_bd_CD', 'aB_Ac_bD_cd', 'aB_Ac_bD_cD', 'aB_Ac_bD_Cd', 'aB_Ac_bD_CD', 'aB_Ac_Bd_cd', 'aB_Ac_Bd_cD', 'aB_Ac_Bd_Cd', 'aB_Ac_Bd_CD', 'aB_Ac_BD_cd', 'aB_Ac_BD_cD', 'aB_Ac_BD_Cd', 'aB_Ac_BD_CD', 'aB_AC_bd_cd', 'aB_AC_bd_cD', 'aB_AC_bd_Cd', 'aB_AC_bd_CD', 'aB_AC_bD_cd', 'aB_AC_bD_cD', 'aB_AC_bD_Cd', 'aB_AC_bD_CD', 'aB_AC_Bd_cd', 'aB_AC_Bd_cD', 'aB_AC_Bd_Cd', 'aB_AC_Bd_CD', 'aB_AC_BD_cd', 'aB_AC_BD_cD', 'aB_AC_BD_Cd', 'aB_AC_BD_CD', 'Ab_ac_bd_cd', 'Ab_ac_bd_cD', 'Ab_ac_bd_Cd', 'Ab_ac_bd_CD', 'Ab_ac_bD_cd', 'Ab_ac_bD_cD', 'Ab_ac_bD_Cd', 'Ab_ac_bD_CD', 'Ab_ac_Bd_cd', 'Ab_ac_Bd_cD', 'Ab_ac_Bd_Cd', 'Ab_ac_Bd_CD', 'Ab_ac_BD_cd', 'Ab_ac_BD_cD', 'Ab_ac_BD_Cd', 'Ab_ac_BD_CD', 'Ab_aC_bd_cd', 'Ab_aC_bd_cD', 'Ab_aC_bd_Cd', 'Ab_aC_bd_CD', 'Ab_aC_bD_cd', 'Ab_aC_bD_cD', 'Ab_aC_bD_Cd', 'Ab_aC_bD_CD', 'Ab_aC_Bd_cd', 'Ab_aC_Bd_cD', 'Ab_aC_Bd_Cd', 'Ab_aC_Bd_CD', 'Ab_aC_BD_cd', 'Ab_aC_BD_cD', 'Ab_aC_BD_Cd', 'Ab_aC_BD_CD', 'Ab_Ac_bd_cd', 'Ab_Ac_bd_cD', 'Ab_Ac_bd_Cd', 'Ab_Ac_bd_CD', 'Ab_Ac_bD_cd', 'Ab_Ac_bD_cD', 'Ab_Ac_bD_Cd', 'Ab_Ac_bD_CD', 'Ab_Ac_Bd_cd', 'Ab_Ac_Bd_cD', 'Ab_Ac_Bd_Cd', 'Ab_Ac_Bd_CD', 'Ab_Ac_BD_cd', 'Ab_Ac_BD_cD', 'Ab_Ac_BD_Cd', 'Ab_Ac_BD_CD', 'Ab_AC_bd_cd', 'Ab_AC_bd_cD', 'Ab_AC_bd_Cd', 'Ab_AC_bd_CD', 'Ab_AC_bD_cd', 'Ab_AC_bD_cD', 'Ab_AC_bD_Cd', 'Ab_AC_bD_CD', 'Ab_AC_Bd_cd', 'Ab_AC_Bd_cD', 'Ab_AC_Bd_Cd', 'Ab_AC_Bd_CD', 'Ab_AC_BD_cd', 'Ab_AC_BD_cD', 'Ab_AC_BD_Cd', 'Ab_AC_BD_CD', 'AB_ac_bd_cd', 'AB_ac_bd_cD', 'AB_ac_bd_Cd', 'AB_ac_bd_CD', 'AB_ac_bD_cd', 'AB_ac_bD_cD', 'AB_ac_bD_Cd', 'AB_ac_bD_CD', 'AB_ac_Bd_cd', 'AB_ac_Bd_cD', 'AB_ac_Bd_Cd', 'AB_ac_Bd_CD', 'AB_ac_BD_cd', 'AB_ac_BD_cD', 'AB_ac_BD_Cd', 'AB_ac_BD_CD', 'AB_aC_bd_cd', 'AB_aC_bd_cD', 'AB_aC_bd_Cd', 'AB_aC_bd_CD', 'AB_aC_bD_cd', 'AB_aC_bD_cD', 'AB_aC_bD_Cd', 'AB_aC_bD_CD', 'AB_aC_Bd_cd', 'AB_aC_Bd_cD', 'AB_aC_Bd_Cd', 'AB_aC_Bd_CD', 'AB_aC_BD_cd', 'AB_aC_BD_cD', 'AB_aC_BD_Cd', 'AB_aC_BD_CD', 'AB_Ac_bd_cd', 'AB_Ac_bd_cD', 'AB_Ac_bd_Cd', 'AB_Ac_bd_CD', 'AB_Ac_bD_cd', 'AB_Ac_bD_cD', 'AB_Ac_bD_Cd', 'AB_Ac_bD_CD', 'AB_Ac_Bd_cd', 'AB_Ac_Bd_cD', 'AB_Ac_Bd_Cd', 'AB_Ac_Bd_CD', 'AB_Ac_BD_cd', 'AB_Ac_BD_cD', 'AB_Ac_BD_Cd', 'AB_Ac_BD_CD', 'AB_AC_bd_cd', 'AB_AC_bd_cD', 'AB_AC_bd_Cd', 'AB_AC_bd_CD', 'AB_AC_bD_cd', 'AB_AC_bD_cD', 'AB_AC_bD_Cd', 'AB_AC_bD_CD', 'AB_AC_Bd_cd', 'AB_AC_Bd_cD', 'AB_AC_Bd_Cd', 'AB_AC_Bd_CD', 'AB_AC_BD_cd', 'AB_AC_BD_cD', 'AB_AC_BD_Cd', 'AB_AC_BD_CD']

help = "\t\033[1;31;40mTakes 4 arguments:\n\t\ti) one model specifier\n\t\tii) one indicator of migration\n\t\tiii) a file containing the posterior distribution\n\t\tiv) a number of replicates for each combination of parameters\033[0m\n\n\t"
help += "\033[1;32;40mAccepted model specifiers:\033[0m\n\t\t"
help += "\n\t\t".join(accepted_models)
help += "\n\t\t\t\033[1;32;40m1M\033[0m\tAll loci share the same M13, M31, M24 and M42\033[0m\n"
help += "\t\t\t\033[1;32;40m1N\033[0m\tAll loci share the same N1, N2, N3, N4, Na_12, Na_34 and Na\033[0m\n"
help += "\t\t\t\033[1;32;40m2M\033[0m\tAll loci have different values for M13, M31, M24 and M42 (Beta distributed)\033[0m\n"
help += "\t\t\t\033[1;32;40m2N\033[0m\tAll loci have different N1, N2, N3, N4, Na_12, Na_34 and Na (Beta distributed)\033[0m\n\n"
help += "\t\033[1;32;40mAccepted indicators of migration:\033[0m\n\t\t"
help += "\n\t\t".join(accepted_migrations[0:5]) + "\tetc ..."
help += "\n"
help += "\t\t\t\033[1;32;40mac\033[0m\tIsolation between pop_A and pop_C, A\033[0m" + ' |' + " C\n"
help += "\t\t\t\033[1;32;40mAc\033[0m\tSecondary contact between pop_A and pop_C, unidirectional A\033[0m" + u'\u2190' + " C\n"
help += "\t\t\t\033[1;32;40maC\033[0m\tSecondary contact between pop_A and pop_C, unidirectional A\033[0m" + u'\u2192' + " C\n"
help += "\t\t\t\033[1;32;40mAC\033[0m\tSecondary contact between pop_A and pop_C, bidirectional  A\033[0m" + u'\u2194' + " C\n"
help += "\t\t\033[1;32;40mExample: ./ABC_4pop.py 2M_2N ab_AC_bd_Cd posterior_file.txt 5\033[0m\n"
help += "\t\t\033[1;32;40mThis will run multilocus simulations with a secondary contact between population A and C (M hetero + N hetero), and unidirectional gene flow from D to C by using parameters values sampled from the posterior distribution into posterior_file.txt, replicated 5 times.\033[0m\n"


help_bpfile = "\t\033[1;31;40mA file named bpfile is needed\033[0m\n\n\t"
help_bpfile += "This file has to contain:"
help_bpfile += "\n\t\tline 1: a useless header with usefull informations for the investigator\n\t\t\t(name of the populations, assumed mutation or recombination rates), just in case..."
help_bpfile += "\n\t\tline 2: length of the loci in nucleotides. One locus per column, one tab between columns."
help_bpfile += "\n\t\tline 3: number of sequences sampled from species A for each locus. One locus per column."
help_bpfile += "\n\t\tline 4: number of sequences sampled from species B for each locus. One locus per column."
help_bpfile += "\n\t\tline 5: number of sequences sampled from species C for each locus. One locus per column."
help_bpfile += "\n\t\tline 6: number of sequences sampled from species D for each locus. One locus per column."
help_bpfile += "\n\t\tline 7: one theta value per locus, with theta = 4.N.L.µ where N is an arbitrary value used as a reference\n\t\t\tto convert coalescent units in demographic units (i.e: N=100,000).\n\t\t\tL is the locus length attributed line_2. µ is the per nucleotide mutation rate."
help_bpfile += "\n\t\tline 8: one rho value per locus, with rho = 4.N.L.r where is the per nucleotide recombination rate.\n"


help_pypy = "\n\t\033[1;31;40mpypy is needed by mscalc_4pop.py\033[0m\n\n\t"
help_pypy += "can be found here: https://pypy.org/\n"


help_msnsam = "\n\t\033[1;31;40mmsnsam is needed. Used as coalescent simulator.\033[0m"
help_msnsam += "\n\n\tCan be found here: https://github.com/rossibarra/msnsam"
help_msnsam += "\n\n\tCan be simply installed with: ./clms"
help_msnsam += "\n\n\tThen you can link the msnsam binary to your bin (ex: sudo ln -s /home/Toto/softwares/msnsam/msnsam /usr/bin/)\n"


help_priorgen = "\n\t\033[1;31;40mpriorgen_4pop.py is needed.\033[0m"
help_priorgen += "\n\n\tYou can link priorgen to your bin as following (for ex): sudo ln -s /home/Toto/softwares/ABC_4pop/priorgen_4pop.py /usr/bin/)\n"


help_mscalc = "\n\t\033[1;31;40mmscalc_4pop.py is needed.\033[0m"
help_mscalc += "\n\n\tYou can link mscalc to your bin as following (for ex): sudo ln -s /home/Toto/softwares/ABC_4pop/mscalc_4pop.py /usr/bin/)\n"


if os.path.isfile('bpfile') == False:
	print(help_bpfile)
	sys.exit()


my_command = "pypy"
test_pypy = any(os.access(os.path.join(path, my_command), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))
if test_pypy == False:
	print(help_pypy)
	sys.exit()


my_command = "msnsam"
test_msnsam = any(os.access(os.path.join(path, my_command), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))
if test_msnsam == False:
	print(help_msnsam)
	sys.exit()


my_command = "priorgen_4pop_GoF_posterior.py"
test_priorgen = any(os.access(os.path.join(path, my_command), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))
if test_priorgen == False:
	print(help_priorgen)
	sys.exit()


my_command = "mscalc_4pop.py"
test_mscalc = any(os.access(os.path.join(path, my_command), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))
if test_mscalc == False:
	print(help_mscalc)
	sys.exit()


if len(sys.argv) != 5:
	print(help)
	sys.exit()


model = sys.argv[1]
migration = sys.argv[2]
posteriorFile = sys.argv[3]
nReplicates = int(sys.argv[4])

if model not in accepted_models:
	print("\n\t\033[1;31;40mThe model '{0}' is not in the list of accepted models\n\n\t".format(model))
	print(help)
	sys.exit()

if migration not in accepted_migrations:
	print("\n\t\033[1;31;40mThe migration '{0}' is not in the list of accepted patterns of migration\n\n\t".format(migration))
	print(help)
	sys.exit()
else:
	migration = migration.replace("_", " ")

infile = open("bpfile", "r")
tmp = infile.readline()
tmp = infile.readline().strip().split("\t")
infile.close()
nLoci = len(tmp)

# nMultilocus = number of combination of joint-parameters in the posterior distribution
nMultilocus = 0
infile = open(posteriorFile, "r")
line = infile.readline()
for i in infile:
	nMultilocus += 1
infile.close()

print(nLoci*nMultilocus)

commandLine = "priorgen_4pop_GoF_posterior.py {0} {1} {3} {4} | msnsam tbs {2} -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 1 2 0 -em tbs 2 1 0 -ej tbs 2 1 -en tbs 1 tbs -em tbs 3 4 0 -em tbs 4 3 0 -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs | mscalc_4pop.py".format(model, migration, nLoci*nMultilocus*nReplicates, posteriorFile, nReplicates)
print(commandLine)
os.system(commandLine)

# ../priorgen_4pop_GoF_posterior.py 2M_2N ab ac BD cd posterior_flo_txn_mal_ama.txt >prior_from_posterior.txt

