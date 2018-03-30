#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys

accepted_models = ["SI_1N", "SI_2N", "SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N"]
accepted_migrations = ["none", "A", "B", "C", "D", "AC", "BD", "ACBD"]

help = "\t\033[1;31;40mTakes 3 arguments:\n\t\ti) one model specifier\n\t\tii) one indicator of migration\n\t\tiii) a number of multilocus simulations\033[0m\n\n\t"
help += "\033[1;32;40mAccepted model specifiers:\033[0m\n\t\t"
help += "\n\t\t".join(accepted_models)
help += "\n\t\t\t\033[1;32;40m1M\033[0m\tAll loci share the same M13, M31, M24 and M42\033[0m\n"
help += "\t\t\t\033[1;32;40m1N\033[0m\tAll loci share the same N1, N2, N3, N4, Na_12, Na_34 and Na\033[0m\n"
help += "\t\t\t\033[1;32;40m2M\033[0m\tAll loci have different values for M13, M31, M24 and M42 (Beta distributed)\033[0m\n"
help += "\t\t\t\033[1;32;40m2N\033[0m\tAll loci have different N1, N2, N3, N4, Na_12, Na_34 and Na (Beta distributed)\033[0m\n\n"
help += "\t\033[1;32;40mAccepted indicators of migration:\033[0m\n\t\t"
help += "\n\t\t".join(accepted_migrations)
help += "\n"
help += "\t\t\t\033[1;32;40mA\033[0m\t\tSecondary contact between pop_A and pop_C, unidirectional A\033[0m" + u'\u2190' + " C\n"
help += "\t\t\t\033[1;32;40mB\033[0m\t\tSecondary contact between pop_B and pop_D, unidirectional B\033[0m" + u'\u2190' + " D\n"
help += "\t\t\t\033[1;32;40mC\033[0m\t\tSecondary contact between pop_A and pop_C, unidirectional A\033[0m" + u'\u2192' + " C\n"
help += "\t\t\t\033[1;32;40mD\033[0m\t\tSecondary contact between pop_B and pop_D, unidirectional B\033[0m" + u'\u2192' + " D\n"
help += "\t\t\t\033[1;32;40mAC\033[0m\t\tSecondary contact between pop_A and pop_C, bidirectional  A\033[0m" + u'\u2194' + " C\n"
help += "\t\t\t\033[1;32;40mBD\033[0m\t\tSecondary contact between pop_B and pop_D, bidirectional  B\033[0m" + u'\u2194' + " D\n"
help += "\t\t\t\033[1;32;40mACBD\033[0m\t\tSecondary contact between A\033[0m" + u'\u2194' + " C and B" + u'\u2194' + " D\n\n"
help += "\t\t\033[1;32;40mExample: ./ABC_4pop.py SC_2M_2N AC 10\033[0m\n"
help += "\t\t\033[1;32;40mThis will run 10 multilocus simulations with a secondary contact between population A and C (M hetero + N hetero)\033[0m\n"


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


my_command = "priorgen_4pop.py"
test_priorgen = any(os.access(os.path.join(path, my_command), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))
if test_priorgen == False:
	print(help_priorgen)
	sys.exit()


my_command = "mscalc_4pop.py"
test_mscalc = any(os.access(os.path.join(path, my_command), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))
if test_mscalc == False:
	print(help_mscalc)
	sys.exit()


if len(sys.argv) != 4:
	print(help)
	sys.exit()


model = sys.argv[1]
migration = sys.argv[2]
nMultilocus = int(sys.argv[3])

if model not in accepted_models:
	print("\n\t\033[1;31;40mThe model '{0}' is not in the list of accepted models\n\n\t".format(model))
	print(help)
	sys.exit()

if migration not in accepted_migrations:
	print("\n\t\033[1;31;40mThe migration '{0}' is not in the list of accepted patterns of migration\n\n\t".format(migration))
	print(help)
	sys.exit()


infile = open("bpfile", "r")
tmp = infile.readline()
tmp = infile.readline().strip().split("\t")
infile.close()
nLoci = len(tmp)

commandLine = "priorgen_4pop.py {0} {1} {2} | msnsam tbs {3} -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs | mscalc_4pop.py".format(model, migration, nMultilocus, nLoci*nMultilocus)
os.system(commandLine)


