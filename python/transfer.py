import pandas
import numpy
# data types are preserved... all nums are numpy.float64
df1 = pandas.read_csv("../all_pep_seqs/All_peptides_Set1.csv")
df2 = pandas.read_csv("../all_pep_seqs/All_peptides_Set2.csv")
df3 = pandas.read_csv("../all_pep_seqs/All_peptides_Set3.csv")

# drop irrelevant data columns
df1.drop(['S0', 'S1', 'S2', 'S3'], axis=1, inplace=True)
df2.drop(['S0', 'S1', 'S2', 'S3'], axis=1, inplace=True)
df3.drop(['S0', 'S1', 'S2', 'S3'], axis=1, inplace=True)

# series of inner merges (df123 retains only the peptides in all sets)
df12 = pandas.merge(df1, df2, on='AA_seq', how='inner')
df13 = pandas.merge(df1, df3, on='AA_seq', how='inner')
# in df123, _x_x and _x_y are from set1, _y_x is from set2, _y_y is from set3
df123 = pandas.merge(df12, df13, on='AA_seq', how='inner')

# returns the total number of copies of the given peptide that were found in
# the detergent washes for the given NGS set
def sumrow(inputrow, ngsset):
    if ngsset == 'one':
        return inputrow.CE_x_x + inputrow.CP1_x_x + inputrow.CP2_x_x + inputrow.CP3_x_x
    elif ngsset == 'two':
        return inputrow.CE_y_x + inputrow.CP1_y_x + inputrow.CP2_y_x + inputrow.CP3_y_x
    else:
        return inputrow.CE_y_y + inputrow.CP1_y_y + inputrow.CP2_y_y + inputrow.CP3_y_y

# returns an output row if data is consistent (at least 2/3 have binding affinity 
# within 10%)
def getOutputRow(inputrow):
    one_two = inputrow.binding_affinity_x_x / inputrow.binding_affinity_y_x
    two_three = inputrow.binding_affinity_y_x / inputrow.binding_affinity_y_y
    one_three = inputrow.binding_affinity_x_x / inputrow.binding_affinity_y_y
    # binding affinity sum
    total_ba = 0
    # sum of copies of peptide found in washes
    total_pres = 0
    # number of datapoints contributing to the sums
    n = 0
    
    if one_two < 1.1 and one_two > 0.9:
        total_ba += inputrow.binding_affinity_x_x + inputrow.binding_affinity_y_x
        total_pres += sumrow(inputrow, 'one') + sumrow(inputrow, 'two')
        n += 2
        if two_three < 1.1 and two_three > 0.9 and one_three < 1.1 and one_three > 0.9:
            total_ba += inputrow.binding_affinity_y_y
            total_pres += sumrow(inputrow, 'three')
            n += 1
    elif two_three < 1.1 and two_three > 0.9:
        total_ba += inputrow.binding_affinity_y_x + inputrow.binding_affinity_y_y
        total_pres += sumrow(inputrow, 'two') + sumrow(inputrow, 'three')
        n += 2
    elif one_three < 1.1 and one_three > 0.9:
        total_ba += inputrow.binding_affinity_x_x + inputrow.binding_affinity_y_y
        total_pres += sumrow(inputrow, 'one') + sumrow(inputrow, 'three')
        n += 2
    
    if n == 0:
        return None
    else:
        return {'AA_seq': inputrow.AA_seq,
                'binding_affinity': total_ba / n,
                'frequency': total_pres / n
                }

output = []
for row in df123.itertuples():
    res = getOutputRow(row)
    if res is not None:
        output.append(res)
        
output = pandas.DataFrame(output)
output.to_csv(path_or_buf='merged.csv',index=False)