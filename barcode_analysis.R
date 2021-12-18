#Hana Wasserman - 12/17/2021
#File containing longest common substring scores between 3 candidate primers and a set of randomly generated barcodes (L18HD3. Mielko, Z.) in forward and reverse orientation
df = read.csv("/Users/hanawasserman/Desktop/GordanLab/Solutions_15/L18HD3_combined_fr_lcs2.csv", header = TRUE, sep = ",", col.names = c("bar_seq", "p1", "p2", "p3","orient"))

d = unique(df)
row.names(d) = NULL
#After dropping duplicates ~8.9% entries are left

dp1 = d[c(1,2,5)]
dp2 = d[c(1,3,5)]
dp3 = d[c(1,4,5)]

dp1["primer"] = "p1"
dp2["primer"] = "p2"
dp3["primer"] = "p3"

names(dp1)[names(dp1) == "p1"] <- "score"
names(dp2)[names(dp2) == "p2"] <- "score"
names(dp3)[names(dp3) == "p3"] <- "score"
d_all = rbind(dp1,dp2,dp3)

#Consistent distribution of LCS2 scores across the three primer candidates in both forward and reverse orientations 
ggplot(data = d_all, aes(x = primer, y = score)) + geom_boxplot(aes(colour = orient))

#Find the primer pair that yields a sufficient number of unique barcodes (~7k) that meet the LCS2 threshold requirement in both the forward and reverse primer directions 
d12 = d[which(d$p1<8 & d$p2<8),]
length(unique(d12$bar_seq))

d13 = d[which(d$p1<8 & d$p3<8),]
length(unique(d13$bar_seq))

d23 = d[which(d$p2<8 & d$p3<8),]
length(unique(d23$bar_seq))

#The lowest LCS2 threshold that yields >7k unique barcodes is 8
#We get the most unique barcodes by pairing p1 and p3 but p1 and p3 also are fairly similar (LCS=13 vs. LCS=2 for the other primer pairings)
#Recommend going with barcode set yielded from pairing primer1 and primer2 