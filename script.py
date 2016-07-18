import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

#read data file from the clinvar ftp
data = pd.read_table("ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/summary_of_conflicting_interpretations.txt")
print data.head()

# Replace NaN with an empty string
data.fillna("")

# Replace NaN with "not provided
data.fillna("not provided", inplace=True)

# how many confilct total?
num_record = data.shape[0]
print (num_record, "conflict")

# how many genes have conflicting interpretaions?
unique_conflict_gene = set(data["#Gene_Symbol"])
num_conflict_gene = len(unique_conflict_gene)
print (num_conflict_gene, " genes have conflicting interpretations")

# how many labs/groups have submitted with at least 1 confict
print (
len(set(list(data["Submitter1"]) + list(data["Submitter2"]))), "labs/groups have submitted with at least 1 conflict")

# how many conflicts of each type?
conflict_type_counts = data["Rank_diff"].value_counts()
print ("the number of conflict of each type")
print (conflict_type_counts)

# Which submitter has the most conflict?
submitter = list(data["Submitter1"]) + list(data["Submitter2"])
individual_submitter_count = Counter(submitter)
majoy_submitter = max(individual_submitter_count, key=individual_submitter_count.get)
print (majoy_submitter, "has the most conflict")


# Which gene is the most conflicted, if we look at those for with rank_diff in 1-4 range
print ("the genes with the most confilcted:")
print (data["#Gene_Symbol"].value_counts()[0:10])
gene_conflict_count = data["#Gene_Symbol"].value_counts()

data_new = data["#Gene_Symbol"][(data["Rank_diff"] >= 1) & (data["Rank_diff"] <= 4)]
print ("the genes with the most conflicted when looking at rank in 1-4 range:")
print (data_new.value_counts()[0:10])

# Which gene is most conflicted if we look at 3-4 range?
data_new = data["#Gene_Symbol"][(data["Rank_diff"] >= 3) & (data["Rank_diff"] <= 4)]
print ("the genes with the most conflicted when looking at rank in 3-4 range:")
print (data_new.value_counts()[0:10])

# What gene has the most egregious conflicts?
data_new = data["#Gene_Symbol"][data["Rank_diff"] == 4]
print ("These genes has the most egregious confilict:")
print (data_new.value_counts()[0:10])

#get the gene table from NCBI
data_gene_table = pd.read_table("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene")
print (data_gene_table.head())

data_new = pd.merge(data, data_gene_table, how='left', left_on='#Gene_Symbol', right_on='Symbol')

data_no_null = data_new[data_new["RSG"].notnull()]
print ("the shape of the new data with gene ID:", data_no_null.shape)
print (data_no_null.shape[0])

# create a new coloum for the length of the gene with missing values
#data_no_null['Gene_length'] = np.nan

print (data_no_null.shape)
print (data_no_null.head())


# get gene length from NCBI reference sequence
import re

data_gene_ref = pd.read_table(
    "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/GCF_000001405.28_refseqgene_alignments.gff3", skiprows=3,
    header=None)

data_gene_ref_shape = data_gene_ref.shape[0]
print ("there is %d genes in the reference table" % (data_gene_ref_shape))
print (type(data_gene_ref[8]))

#get the gene ID using regular expression
# Target=NG_027734.1

gene_id_length = {}
for i in range(0, data_gene_ref_shape):
    gene_length = data_gene_ref[4][i] - data_gene_ref[3][i]
    match = re.search("""Target=(NG_\d{6}\.\d{1,2})""", str(data_gene_ref[8][i]))
    if match:
        genome_id = match.group(1)
    if genome_id in gene_id_length:
        if gene_id_length[genome_id] < gene_length:
            gene_id_length[genome_id] = gene_length
    else:
        gene_id_length[genome_id] = gene_length


df_data_id_length=pd.DataFrame(gene_id_length.items(), columns=['Gene_ID', 'Gene_length'])
print (df_data_id_length.head())

#merge data_not_null with df_data_id_length
data_new = pd.merge(data_no_null, df_data_id_length, how='left', left_on='RSG', right_on='Gene_ID')

print (data_new.head())

length_conflict = data_new["Gene_length"].value_counts()
#print length_conflict
print ("data type:", type(length_conflict))
length_conflict_1=Counter(data_new["Gene_length"])
#for g in sorted(length_conflict_1, key=length_conflict_1.get, reverse=True):
#    print g, length_conflict_1[g]

print ("number of conflict max, median, min:",np.max(length_conflict_1.values()),np.median(length_conflict_1.values()),np.min(length_conflict_1.values()))
print ("gene length max, median, min:", np.max(length_conflict_1.keys()), np.median(length_conflict_1.keys()), np.min(length_conflict_1.keys()))


# scatter plot of the gene length vs number of conflict
plt.scatter(length_conflict_1.keys(),length_conflict_1.values())
plt.xlim(0,2500000)
plt.ylim(0,4000)
plt.xlabel('gene length')
plt.ylabel('number of conflict')
plt.show()


#2.	What gene has the most conflict given length?
gene_id_conflict_count_len={}
gene_id_conflict_count = Counter(data_new["RSG"])
for id in gene_id_conflict_count:
    gene_id_conflict_count_len[id]=float(gene_id_conflict_count[id])/(gene_id_length[id])

for id in sorted(gene_id_conflict_count_len,key=gene_id_conflict_count_len.get, reverse=True):
    print id,gene_id_conflict_count_len[id]
    if id == "NG_011618.3":
        print "OKOKOKOKTTN"


print data_gene_table["Symbol"][data_gene_table["RSG"] == 'NG_007111.1']
print gene_id_conflict_count_len['NG_011618.3']

# explore the total submission and total alleles vs conflict
data_summary=pd.read_table("ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_specific_summary.txt",skiprows=1)
print (data_summary.head())


df_gene_conflict_count=pd.DataFrame(Counter(data_new["#Gene_Symbol"]).items(),columns=['Gene_Symbol','Conflict_Count'])
data_summary_new = pd.merge(df_gene_conflict_count,data_summary, how='left', left_on='Gene_Symbol', right_on='#Symbol')

print ("Total submission max, median, and min", np.max(data_summary_new['Total_submissions']),np.median(data_summary_new["Total_submissions"]),np.min(data_summary_new["Total_submissions"]))
print ("Total Alleles max, median, and min", np.max(data_summary_new['Total_alleles']),np.median(data_summary_new["Total_alleles"]),np.min(data_summary_new["Total_alleles"]))

print "Largest total submission:"
print data_summary_new["Gene_Symbol"][data_summary_new["Total_submissions"] == np.max(data_summary_new['Total_submissions'])]
print "Smallest total submission:"
print data_summary_new["Gene_Symbol"][data_summary_new["Total_submissions"] == np.min(data_summary_new['Total_submissions'])]
print "Top total submission:"
print data_summary_new["Gene_Symbol"][data_summary_new["Total_submissions"] > 2000]

print "Largest total Alleles:"
print data_summary_new["Gene_Symbol"][data_summary_new["Total_alleles"] == np.max(data_summary_new['Total_alleles'])]
print "Smallest total Alleles:"
print data_summary_new["Gene_Symbol"][data_summary_new["Total_alleles"] == np.min(data_summary_new['Total_alleles'])]


plt.hist(data_summary_new["Total_submissions"], bins='auto')  # plt.hist passes it's arguments to np.histogram
plt.title("Number of total submissions")
#plt.xlim(0,4000)
plt.show()

plt.hist(data_summary_new["Total_alleles"], bins='auto')  # plt.hist passes it's arguments to np.histogram
plt.title("Number of total alleles")
plt.xlim(0,3000)
plt.show()


plt.scatter(data_summary_new['Total_submissions'],data_summary_new["Conflict_Count"])
plt.xlim(0,12000)
plt.ylim(0,4000)
plt.xlabel('Total Submission')
plt.ylabel('number of conflict')
plt.show()

plt.scatter(data_summary_new['Total_alleles'],data_summary_new["Conflict_Count"])
plt.xlim(0,2000)
plt.ylim(0,1000)
plt.xlabel('Total Alleles')
plt.ylabel('number of conflict')
plt.show()
