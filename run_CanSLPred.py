'''import subprocess
test = subprocess.Popen(["python3","iFeature/iFeature.py"], stdout=subprocess.PIPE)
output = test.communicate()[0]
print(output)
import os
import numpy as np
filename = "example"
iFeature_descriptor_list = ["AAC", "PAAC", "APAAC", "DPC", "GAAC", "CKSAAP", "CKSAAGP", "GDPC", "Moran", "Geary",
                            "NMBroto", "CTDC", "CTDD", "CTDT", "CTriad", "KSCTriad", "SOCNumber", "QSOrder"]
for descr in iFeature_descriptor_list:
    os.system("python3 iFeature/iFeature.py --file fasta_files/" +filename+".fasta --type " + descr+\
              " --out iFeature_descriptors_results/"+filename+"_"+descr+".txt")
POSSUM_descriptor_list = ["aac_pssm", "d_fpssm", "smoothed_pssm", "ab_pssm", "pssm_composition", "rpm_pssm",
                          "s_fpssm", "dpc_pssm", "k_separated_bigrams_pssm", "tri_gram_pssm", "eedp", "tpc",
                          "edp", "rpssm", "pse_pssm", "dp_pssm", "pssm_ac", "pssm_cc", "aadp_pssm", "aatp", "medp"]
for possum_descr in POSSUM_descriptor_list:
    os.system("python2 POSSUM_Standalone_Toolkit/src/possum.py -i fasta_files/"+filename+\
              ".fasta -o POSSUM_descriptors_results/"+filename+"_"+possum_descr+".txt -t "+ possum_descr +
              " -p fasta_files/pssm_files")'''
import sys
print(sys.maxsize)
