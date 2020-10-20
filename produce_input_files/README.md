# converts the fasta file in input files for ABC
./fasta2bpfile.py chr_10.MSA.fasta ADR INT INV MAC 5000 0.2 4 0.00000001 0.0000000001 100000  
_argument 1_: name of the fasta file (ex: chr_10.MSA.fasta)  
_argument 2_: name for species A (ex: ADR)  
_argument 3_: name for species B (ex: INT)  
_argument 4_: name for species C (ex: INV)  
_argument 5_: name for species D (ex: MAC)  
_argument 6_: size of the bin in bp (ex: 5000)  
_argument 7_: max proportion of missing data in the alignement (ex: 0.2)  
_argument 8_: minimum number of individuals in a given species at a given locus. (ex: 4)  
_argument 9_: mutation rate mu per nucleotide and per generation (ex: 0.000000001)  
_argument 10_: recombinaration rate (ex: 0.0000000001 --> keep it small for big bins, i.e, 0.05mu  
_argument 11_: size of the reference population (Nref) arbitrary fixed for ABC (ex: 100000)  
produces 3 output files:  
         . general_infos.txt : informations about each bins (position, length, number of SNPs, etc...)  
         . observed_data.ms : the observed sequenced data in ms's format  
         . bpfile : used with observed_data.ms by mscalc to compute statistics  
  
# computes summary statistics from the converted files  
cat observed_data.ms | pypy mscalc_4pop_observed_data.py  
produces 2 output files:  
     . ABCstat.txt : avg and std statistics over all bins  
     . ABCstat_bins.txt : statistics for each bin  
  
