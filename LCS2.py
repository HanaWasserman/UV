#Loop command:
#for i in L18HD3*.txt;do python3 LCS2.py $i;done

#LCS2
import sys
import pandas as pd
import pylcs

input_txt = sys.argv[1]
print(input_txt)
f = pd.read_table(input_txt, header = None)
#f = pd.read_table("/Users/hanawasserman/Downloads/L18HD3.txt", header = None)

#candidate primer sequences
p1 = "CACATACGCATACACATACGCACATACGCACACGC"
p2 = "GTATGTATGCGTGTGTGTATGCGTGCGTGTATGCG"
p3 = "TACGCATACACGCACACGCACATACACATACATAC"

#reversed candidate primer sequences
p1_r = "CGCACACGCATACACGCATACACATACGCATACAC"
p2_r = "GCGTATGTGCGTGCGTATGTGTGTGCGTATGTATG"
p3_r = "CATACATACACATACACGCACACGCACATACGCAT"

def do_align1(x):
	return pylcs.lcs2(x, p1)

def do_align2(x):
	return pylcs.lcs2(x, p2)

def do_align3(x):
	return pylcs.lcs2(x, p3)


def do_align1_r(x):
	return pylcs.lcs2(x, p1_r)

def do_align2_r(x):
	return pylcs.lcs2(x, p2_r)

def do_align3_r(x):
	return pylcs.lcs2(x, p3_r)


scores1 = f.applymap(do_align1)
scores2 = f.applymap(do_align2)
scores3 = f.applymap(do_align3)

scores1_r = f.applymap(do_align1_r)
scores2_r = f.applymap(do_align2_r)
scores3_r = f.applymap(do_align3_r)

f = f.rename(columns = {0:"barcode_seq"})

scores1 = scores1.rename(columns = {0:"p1_CACATACGCATACACATACGCACATACGCACACGC"})
scores2 = scores2.rename(columns = {0:"p2_GTATGTATGCGTGTGTGTATGCGTGCGTGTATGCG"})
scores3 = scores3.rename(columns = {0:"p3_TACGCATACACGCACACGCACATACACATACATAC"})

scores1_r = scores1_r.rename(columns = {0:"p1_CACATACGCATACACATACGCACATACGCACACGC"})
scores2_r = scores2_r.rename(columns = {0:"p2_GTATGTATGCGTGTGTGTATGCGTGCGTGTATGCG"})
scores3_r = scores3_r.rename(columns = {0:"p3_TACGCATACACGCACACGCACATACACATACATAC"})

res_f = pd.concat([f,scores1,scores2,scores3],axis=1)
res_f["orient"] = "f"

res_r = pd.concat([f,scores1_r,scores2_r,scores3_r],axis=1)
res_r["orient"] = "r"

res = pd.concat([res_f,res_r], axis=0)

input_arr = input_txt.split(".")
res.to_csv(input_arr[0]+"_fr_lcs2.csv",sep=',',index=False)

