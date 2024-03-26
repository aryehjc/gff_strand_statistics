import pandas as pd
import numpy as np
# remove characters after and including ##FASTA, header and all ##sequence-region.
#run the two bash scripts below on the gff file
# run these on sap1 and sap2 resistant, feb 28 2023
#run sed '/##FASTA/,$d' test.gff > out.gff and sed '/^##/d' out.gff > out.csv in bash
#this is an example, change it later to run on gffs for RN4220, S1 S2..etc
#with open('out.csv', 'r') as f:
#    data=f.read().replace('\n', ' ')
#Reads the CSV and counts number of + strands.. (lines with + and are going to be in updated_plus_strands_list later)..CHANGE THIS LATER
#data.split(' ').count("+")
#print('Number of forward strands:')
#print(data.count("+"))
# If there is an issue with tab separation and you obtain the wrong value or a blank CSV, fix it by copypasting columns directly into input csv file below.
# /home/aryeh/Staphylococcus_Aureus_Assembly_Amended_Params/ASSEMBLIES_DEC_2023/SAP1_Resistant/Final_Contigs/spades_output_directory/Final_SPADES_SA_Annotation_Contigs_2023/out.csv

Annot = pd.read_csv("out.csv", sep='\t', header=None) # Follow above steps to convert .gff to .csv, and insert as input file here. (out.csv)
Annot.columns = ["Node", "Predictor", "GeneType", "BasePairStart", "BasePairEnd", ".", "Strand", "0", "Annotation"]
print(Annot)
#(1): After I run sed '/##FASTA/,$d' test.gff > out.gff and sed '/^##/d' out.gff > out.csv in bash, make a dataframe Annot out of the .csv file with columns according to the descriptors.

#print(Annot[['Node', 'BasePairStart', 'BasePairEnd', 'Strand', 'Annotation']].head(10))
Annot_New = Annot.filter(['Node', 'BasePairStart', 'BasePairEnd', 'Strand', 'Annotation'], axis=1)
#(2): Filtering the above dataframe Annot according to headers I will use for statistics, the node type, base pair start and end (for length) +/- strand and annotation.

#print(Annot_New)

Annot_New['Annotation'] = Annot_New['Annotation'].str.replace('^.*?(?=product)', '', regex=True)
#(3): Filtering extra characters or strings off the annotation header just so I can get the product information.
   
print(Annot_New)
Strand_List_Sorted_By_Ascending = Annot_New.sort_values('BasePairStart')
#(4): By sorting the dataframe Annot_New according to base pair length from smallest to largest, I can see similar strands + and - that can potentially overlap to have TEs.


#Annot_New.to_csv('Complete_Strand_List.csv', sep='\t')
Strand_List_Sorted_By_Ascending.to_csv('Complete_Strand_List.csv', sep='\t')
#(5): I then convert the Steand List Ascending dataframe to a csv to see these changes.


#list_A = Annot_New['Node'].str.split(',\s*').explode().unique().tolist()
#print(list_A)
# count dataframe 'unique' rows, finding total occurence of each
# see https://stackoverflow.com/questions/19960077/how-to-filter-pandas-dataframe-using-in-and-not-in-like-in-sql
# see https://stackoverflow.com/questions/41665659/count-of-unique-value-in-column-pandas

total_strands_per_chromosome = (Annot_New['Node'].value_counts())
#(6): Printing a list of value counts of each Node, i.e. gives how many total strands for each node, both +/-

#print(total_strands_per_chromosome)
total_strands_per_chromosome.to_csv('total_strands_list.csv', sep='\t')
TotalStrands_List = pd.read_csv("total_strands_list.csv", skiprows = 1, sep='\t', header=None)
#(7): Converting (6) to a csv file to visualise total strands in each node.

TotalStrands_List.columns = ["Node", "Total_Strands"]
#(8): Converting (7) back to a dataframe for operations

#print(Annot_New['Node'].loc[Annot_New['Strand']!='-'].nunique())
#print(Annot_New[Annot_New["Strand"] == "\+"]["Node"].nunique())
# unique values of each containing + 
plus_strands = (Annot_New.loc[Annot_New['Strand'] == "+", 'Node'].value_counts())
#(9): Now I am filtering Annot_New dataframe by WHICH Node values contain + Strands and total values altogether.

print(plus_strands)
plus_strands.to_csv('plus_strands_list.csv', sep='\t')
#(10): Converting above to a CSV to visualise the total count of Plus Strands..
PlusStrands_List = pd.read_csv("plus_strands_list.csv", skiprows = 1, sep='\t', header=None) #correct
PlusStrands_List.columns = ["Node", "Plus_Strands"]
PlusStrands_List.to_csv('updated_plus_strands_list.csv', sep='\t') #correct
#(11): Converting (10) back to a dataframe for operations, now I have a Plus Strand list with counts.

#total strands list --> updated plus strands list, indexes different for left column because arranged in descending order..
# arranging total strands list in order of plus strand.

#==============feb 21=============

#TotalStrands_List['Plus_Strands'] = TotalStrands_List.index.map(dict(zip(PlusStrands_List.index,PlusStrands_List['Node']))).astype(float)
#https://stackoverflow.com/questions/70460010/map-index-of-one-dataframe-to-column-of-another-dataframe
#TotalStrands_List = TotalStrands_List.set_index('Node').join(PlusStrands_List.set_index('Node'))


# https://stackoverflow.com/questions/70460010/map-index-of-one-dataframe-to-column-of-another-dataframe
TotalStrands_List.set_index("Node", inplace=True)
PlusStrands_List.set_index("Node", inplace=True)
#(12): Since I am performing calculations to get % of Plus Strands in each node, the index of the Total Strands and Plus Strands csv must be the same once I perform the operations on each row.
TotalStrands_List['Plus_Strands'] = PlusStrands_List['Plus_Strands']
TotalStrands_List['Plus_Strands'] = TotalStrands_List['Plus_Strands'].astype(float)
TotalStrands_List['Total_Strands'] = TotalStrands_List['Total_Strands'].astype(float)
#(13): Sets the index, adds Plus Strands as a column to Total Strand list, and converts Plus Strands and Total Strands columns to float so I can subtract Plus from Total and get Minus Strands.

TotalStrands_List.to_csv('TestmeMar6.csv', sep='\t')
#This is just for testing

#seems to have worked, read this into new df and drop the index column, and use it to replace TotalStrands_List below?*****

# 'add' missing values to updated..from unsorted.. convert unsorted and updated to a df each and concat?, wanting to remove the top empty row.
# https://stackoverflow.com/questions/72946006/pandas-fill-in-missing-index-from-other-dataframe
#if value in column PlusStrandsList not in TotalStrandsList
#remove it from TotalStrandsList
#================feb21==============================
#if column value in updated_plus_strands_list.csv not in last file, remove it

TotalStrands_List['Minus_Strands'] = TotalStrands_List['Total_Strands'] - TotalStrands_List['Plus_Strands']
#(14): Calculate Minus Strands

# need to account for ONLY 1 plus strand and let it also be 100%
print(TotalStrands_List)
TotalStrands_List.to_csv('TestmeMar62.csv', sep='\t')
TotalStrands_List_Amended = pd.read_csv("TestmeMar62.csv", skiprows = 1, sep='\t', header=None) #correct
TotalStrands_List_Amended.columns = ["Node", "Total_Strands", "Plus_Strands", "Minus_Strands"]
print(TotalStrands_List_Amended)
# (15): Now visualising the Total, Plus and Minus Strands for each node.

#TotalStrands_List.loc[TotalStrands_List.Plus_Strands == '', 'Percentage_Plus_Strands'] = '100'
#TotalStrands_List.where(TotalStrands_List.Plus_Strands == "0", 'Plus_Strands')
#====
# https://stackoverflow.com/questions/50685107/pandas-dataframe-nan-values-not-replacing
TotalStrands_List_Amended.Plus_Strands.replace('nan', np.nan, inplace=True)
TotalStrands_List_Amended['Plus_Strands'].fillna(1, inplace=True)
# (16): Some nodes don't contain Minus Strands and only contain 1 Plus Strand, but in the csv it appears as blank for Plus Strands that have only 1 STRAND THAT IS A PLUS STRAND
# (17): Therefore I convert the blank NaN value to 1 (so it becomes 1 plus strand only: i.e. out of 1 total strand, it is a plus strand.. so then it calculates corectly.

#TotalStrands_List_Amended = TotalStrands_List_Amended.set_index('Node',inplace=True)
# this makes sure that if it only contains 1 pos strand it will still calculate properly
TotalStrands_List_Amended['Percentage_Plus_Strands'] = (TotalStrands_List_Amended['Plus_Strands'] / TotalStrands_List_Amended['Total_Strands']) * 100
# (18): Calculating Plus Strands out of Total Strands as a percentage.
print(TotalStrands_List_Amended.columns.tolist())
print(PlusStrands_List.columns.to_list())
print(PlusStrands_List)
# (19): Now I am converting the columns of Total Strands and Plus Strands dataframes to a list, because I want to filter the Total Strands by those that contain Plus Strands only.
PlusStrands_List_Amended = pd.read_csv("updated_plus_strands_list.csv", skiprows = 1, sep='\t', header=None, index_col=0) #correct
# (20): When I use the original Plus Strands List dataframe, it has no headers, so I'm reading the CSV in again as a dataframe and adding headers so that I can filter by Node.
PlusStrands_List_Amended.columns = ["Node", "Plus_Strands"]
print(PlusStrands_List_Amended)
#Had to put 'Node' column back before this worked properly
List_Only_Positive_Strands = TotalStrands_List_Amended[TotalStrands_List_Amended['Node'].isin(PlusStrands_List_Amended['Node'])]
# (21): Now I filter the TotalStrands_List_Amended dataframe (with Total, Plus, Minus and % Plus Strands, by the rows of Nodes only containing + Strands)

#List_Only_Positive_Strands['Plus_Strands'] = List_Only_Positive_Strands['Plus_Strands'].replace(np.nan, 1)
#List_Only_Positive_Strands = List_Only_Positive_Strands.replace(np.nan, 0)
#TotalStrands_List.to_csv('percentage_total_strands_list.csv', sep='\t')

List_Only_Positive_Strands.to_csv('Positive_Strands.csv', sep='\t')
# (22): Converting (21) to a CSV again.

Base_Pairs_Of_Positive_Strands = Strand_List_Sorted_By_Ascending[Strand_List_Sorted_By_Ascending['Node'].isin(List_Only_Positive_Strands['Node'])]
print(Base_Pairs_Of_Positive_Strands)
Base_Pairs_Of_Positive_Strands.to_csv('Base_Pairs_Of_Positive_Strands.csv', sep='\t')
# (23): For ease of visualisation of base pairs of all positive strands in ascending order, after being filtered from the complete strand list
Filtered_Base_Pairs_Of_Positive_Strands = Base_Pairs_Of_Positive_Strands.loc[Base_Pairs_Of_Positive_Strands['Strand'] == "+"]
# (24): + Strands in (23)
print(Filtered_Base_Pairs_Of_Positive_Strands) #next can start to use for graphing
Filtered_Base_Pairs_Of_Positive_Strands.to_csv('Filtered_Base_Pairs_Of_Positive_Strands.csv', sep='\t')
# (25): For subsequent graphing
Filtered_Base_Pairs_Of_Positive_Strands['BasePairStart'] = Filtered_Base_Pairs_Of_Positive_Strands['BasePairStart'].astype(int)
Filtered_Base_Pairs_Of_Positive_Strands['BasePairEnd'] = Filtered_Base_Pairs_Of_Positive_Strands['BasePairEnd'].astype(int)
Filtered_Base_Pairs_Of_Positive_Strands['Range'] = [len((range(i, j+1))) for i, j in Filtered_Base_Pairs_Of_Positive_Strands[['BasePairStart','BasePairEnd']].values]
# (26): Converting Start and End of Base Pair range to int and looping to find the range for each one, to find spanning number of base pairs.
Positive_Strand_Range = Filtered_Base_Pairs_Of_Positive_Strands.assign(freq=Filtered_Base_Pairs_Of_Positive_Strands.groupby('Node')['Node'].transform('count'))\
  .sort_values(by=['freq','Node'],ascending=[False,True])
Positive_Strand_Range.to_csv('1stRangeOfBasePairs_Positive.csv', sep='\t')
# https://stackoverflow.com/questions/56360610/sum-column-based-on-another-column-in-pandas-dataframe
# https://stackoverflow.com/questions/27018622/pandas-groupby-sort-descending-order
dff = Positive_Strand_Range.groupby('Node')['Range'].sum().sort_values(ascending=False)
# (27): Above 2 columns are created, Range and freq, where freq is occurence of each node (frequency) and Range is the corresponding number of base pairs in each node row.
# (28): Following this, we group Positive_Strand_Range by counting the amount of times each row appears, 
# and adding all the 'range' column values in each row, for each node respectively, to obtain spanning base pairs.



#dff2 = dff.sort_values('1')
#Positive_Strand_Range['TotalFreq'] = Positive_Strand_Range.groupby('Node')['Range'].sum()
dff.to_csv('FinalRangeOfBasePairs_Positive.csv', sep='\t') # GIVING RANGE OF FORWARD STRANDS OVER BASE PAIRS
print(dff)
concatenate_columns = pd.read_csv('FinalRangeOfBasePairs_Positive.csv', sep='\t')
# (29): Converting (28) above to CSV.

concatenate_columns.Node = pd.Categorical(concatenate_columns.Node,
                                  categories=List_Only_Positive_Strands.Node.values,
                                  ordered=True)

concatenate_columns = concatenate_columns.sort_values('Node')
print(concatenate_columns)
concatenate_columns.to_csv('Final_SortedRangeOfBasePairs_Positive.csv', sep='\t')
# (30): Making sure the indexes are correct according to Plus_Strand List previously
List_Only_Positive_Strands['Spanning No. of Base Pairs']= concatenate_columns['Range'].set_axis(List_Only_Positive_Strands.index)
# (31): PRESERVE NUMBER ORDER in index as I add it to the final dataframe:
#https://stackoverflow.com/questions/46396257/adding-a-new-column-in-pandas-dataframe-from-another-dataframe-with-differing-in
List_Only_Positive_Strands['Spanning No. of Base Pairs'] = List_Only_Positive_Strands['Spanning No. of Base Pairs'].astype(float)
# (32): ensuring values are in float for future number calculations
print(List_Only_Positive_Strands)
#print(List_Only_Positive_Strands.dtypes)
#print(concatenate_columns.dtypes)
List_Only_Positive_Strands.to_csv('Statistics_Positive_Strands.csv', sep='\t')
# (33): Final statistical output for nodes that contain AT LEAST one Plus Strand.
# This is the file I use
#SEEMED TO WORK MAR 6, SEE BELOW MAYBE IF - STRANDS, perhaps filter minus list with final CSV, etc. Figuring this out on Monday
#===============================================================================================================================================

total_strands_per_chromosome_2 = (Annot_New['Node'].value_counts())
#(34): Printing a list of value counts of each Node again, i.e. gives how many total strands for each node, both +/-, now calculating for - strands.

#print(total_strands_per_chromosome)
total_strands_per_chromosome_2.to_csv('total_strands_list_2.csv', sep='\t')
TotalStrands_List_2 = pd.read_csv("total_strands_list_2.csv", skiprows = 1, sep='\t', header=None)
#(35): Converting (6) to a csv file to visualise total strands in each node.

TotalStrands_List_2.columns = ["Node", "Total_Strands"]


minus_strands = (Annot_New.loc[Annot_New['Strand'] == "-", 'Node'].value_counts())
print(minus_strands)
minus_strands.to_csv('minus_strands_list.csv', sep='\t')
MinusStrands_List = pd.read_csv("minus_strands_list.csv", skiprows = 1, sep='\t', header=None) #correct
MinusStrands_List.columns = ["Node", "Minus_Strands"]
MinusStrands_List.to_csv('updated_minus_strands_list.csv', sep='\t') #correct

TotalStrands_List_2.set_index("Node", inplace=True)
MinusStrands_List.set_index("Node", inplace=True)
#(36): Since I am performing calculations to get % of Minus Strands in each node, the index of the Total Strands and Minus Strands csv must be the same once I perform the operations on each row.
TotalStrands_List_2['Minus_Strands'] = MinusStrands_List['Minus_Strands']
TotalStrands_List_2['Minus_Strands'] = TotalStrands_List_2['Minus_Strands'].astype(float)
TotalStrands_List_2['Total_Strands'] = TotalStrands_List_2['Total_Strands'].astype(float)

TotalStrands_List_2.to_csv('TestmeMar11.csv', sep='\t')




TotalStrands_List_2['Plus_Strands'] = TotalStrands_List_2['Total_Strands'] - TotalStrands_List_2['Minus_Strands']
#(37): Calculate Plus Strands

# need to account for ONLY 1 plus strand and let it also be 100%
print(TotalStrands_List_2)
TotalStrands_List_2.to_csv('TestmeMar112.csv', sep='\t')
TotalStrands_List_2_Amended = pd.read_csv("TestmeMar112.csv", skiprows = 1, sep='\t', header=None) #correct
TotalStrands_List_2_Amended.columns = ["Node", "Total_Strands", "Minus_Strands", "Plus_Strands"]
print(TotalStrands_List_2_Amended)
TotalStrands_List_2_Amended.to_csv("TestmeMar113.csv", sep='\t')

TotalStrands_List_2_Amended.Minus_Strands.replace('nan', np.nan, inplace=True)
TotalStrands_List_2_Amended['Minus_Strands'].fillna(0, inplace=True) #this worked for the previous case if fillna was 1, so anything -1 worked, but now result is NaN when -0.
# I need it to perform the arithmetic operation properly.
# Cast it again?
# (16): Some nodes don't contain Minus Strands and only contain 1 Plus Strand, but in the csv it appears as blank for Plus Strands that have only 1 STRAND THAT IS A PLUS STRAND
# (17): Therefore I convert the blank NaN value to 1 (so it becomes 1 plus strand only: i.e. out of 1 total strand, it is a plus strand.. so then it calculates corectly.

#TotalStrands_List_Amended = TotalStrands_List_Amended.set_index('Node',inplace=True)
# this makes sure that if it only contains 1 pos strand it will still calculate properly
TotalStrands_List_2_Amended['Percentage_Minus_Strands'] = (TotalStrands_List_2_Amended['Minus_Strands'] / TotalStrands_List_2_Amended['Total_Strands']) * 100
print(TotalStrands_List_2_Amended)

print(TotalStrands_List_2_Amended.columns.tolist())
print(MinusStrands_List.columns.to_list())
print(MinusStrands_List)
# (19): Now I am converting the columns of Total Strands and Plus Strands dataframes to a list, because I want to filter the Total Strands by those that contain Plus Strands only.
MinusStrands_List_Amended = pd.read_csv("updated_minus_strands_list.csv", skiprows = 1, sep='\t', header=None, index_col=0) #correct
# (20): When I use the original Plus Strands List dataframe, it has no headers, so I'm reading the CSV in again as a dataframe and adding headers so that I can filter by Node.
MinusStrands_List_Amended.columns = ["Node", "Minus_Strands"]
print(MinusStrands_List_Amended)
#Had to put 'Node' column back before this worked properly
List_Only_Negative_Strands = TotalStrands_List_2_Amended[TotalStrands_List_2_Amended['Node'].isin(MinusStrands_List_Amended['Node'])]
# (21): Now I filter the TotalStrands_List_Amended dataframe (with Total, Plus, Minus and % Plus Strands, by the rows of Nodes only containing + Strands)

#List_Only_Positive_Strands['Plus_Strands'] = List_Only_Positive_Strands['Plus_Strands'].replace(np.nan, 1)
#List_Only_Positive_Strands = List_Only_Positive_Strands.replace(np.nan, 0)
#TotalStrands_List.to_csv('percentage_total_strands_list.csv', sep='\t')

List_Only_Negative_Strands.to_csv('Negative_Strands.csv', sep='\t')

Base_Pairs_Of_Negative_Strands = Strand_List_Sorted_By_Ascending[Strand_List_Sorted_By_Ascending['Node'].isin(List_Only_Negative_Strands['Node'])]
print(Base_Pairs_Of_Negative_Strands)
Base_Pairs_Of_Negative_Strands.to_csv('Base_Pairs_Of_Negative_Strands.csv', sep='\t')

Filtered_Base_Pairs_Of_Negative_Strands = Base_Pairs_Of_Negative_Strands.loc[Base_Pairs_Of_Negative_Strands['Strand'] == "-"]
# (24): + Strands in (23)
print(Filtered_Base_Pairs_Of_Negative_Strands) #next can start to use for graphing
Filtered_Base_Pairs_Of_Negative_Strands.to_csv('Filtered_Base_Pairs_Of_Negative_Strands.csv', sep='\t')


# (25): For subsequent graphing

Filtered_Base_Pairs_Of_Negative_Strands['BasePairStart'] = Filtered_Base_Pairs_Of_Negative_Strands['BasePairStart'].astype(int)
Filtered_Base_Pairs_Of_Negative_Strands['BasePairEnd'] = Filtered_Base_Pairs_Of_Negative_Strands['BasePairEnd'].astype(int)
Filtered_Base_Pairs_Of_Negative_Strands['Range'] = [len((range(i, j+1))) for i, j in Filtered_Base_Pairs_Of_Negative_Strands[['BasePairStart','BasePairEnd']].values]
# (26): Converting Start and End of Base Pair range to int and looping to find the range for each one, to find spanning number of base pairs.
Negative_Strand_Range = Filtered_Base_Pairs_Of_Negative_Strands.assign(freq=Filtered_Base_Pairs_Of_Negative_Strands.groupby('Node')['Node'].transform('count'))\
  .sort_values(by=['freq','Node'],ascending=[False,True])
Negative_Strand_Range.to_csv('1stRangeOfBasePairs_Negative.csv', sep='\t')
# https://stackoverflow.com/questions/56360610/sum-column-based-on-another-column-in-pandas-dataframe
# https://stackoverflow.com/questions/27018622/pandas-groupby-sort-descending-order
dff2 = Negative_Strand_Range.groupby('Node')['Range'].sum().sort_values(ascending=False)

dff2.to_csv('FinalRangeOfBasePairs_Negative.csv', sep='\t') # GIVING RANGE OF FORWARD STRANDS OVER BASE PAIRS
print(dff2)
concatenate_columns_2 = pd.read_csv('FinalRangeOfBasePairs_Negative.csv', sep='\t')
# (29): Converting (28) above to CSV.

concatenate_columns_2.Node = pd.Categorical(concatenate_columns_2.Node,
                                  categories=List_Only_Negative_Strands.Node.values,
                                  ordered=True)

concatenate_columns_2 = concatenate_columns_2.sort_values('Node')
print(concatenate_columns_2)
concatenate_columns_2.to_csv('Final_SortedRangeOfBasePairs_Negative.csv', sep='\t')
# (30): Making sure the indexes are correct according to Plus_Strand List previously
List_Only_Negative_Strands['Spanning No. of Base Pairs']= concatenate_columns_2['Range'].set_axis(List_Only_Negative_Strands.index)
# (31): PRESERVE NUMBER ORDER in index as I add it to the final dataframe:
#https://stackoverflow.com/questions/46396257/adding-a-new-column-in-pandas-dataframe-from-another-dataframe-with-differing-in
List_Only_Negative_Strands['Spanning No. of Base Pairs'] = List_Only_Negative_Strands['Spanning No. of Base Pairs'].astype(float)
# (32): ensuring values are in float for future number calculations
print(List_Only_Negative_Strands)
#print(List_Only_Positive_Strands.dtypes)
#print(concatenate_columns.dtypes)
List_Only_Negative_Strands.to_csv('Statistics_Negative_Strands.csv', sep='\t')
Filtered_Statistics_Negative_Strand_List = List_Only_Negative_Strands.loc[List_Only_Negative_Strands['Plus_Strands'] == 0]
print(Filtered_Statistics_Negative_Strand_List)
Filtered_Statistics_Negative_Strand_List.to_csv('Statistics_Remaining_Negative_Strands.csv', sep='\t')

# Overall, PROKKA shows nodes with only negative strands, but PGAP is more accurate, and there will only be nodes with both +/- strands with PGAP.
# Therefore Statistics_Remaining_Negative_Strands.csv is not produced with PGAP, only the Positive file.
# Done, 11/3/2024, now I have to filter List_Only_Negative_Strands by Plus_Strands = 0, 
# And then I'll have stats for both + and - strands. Test on Mar 12, run on SAP1, SAP2 res, and then 

# now join this Range column to List_Only_Positive_Strands..
# Print above with Positive_Strands.csv (attach there)
# add_back_to_column_positive = pd.read_csv('FinalRangeOfBasePairs_Positive.csv', sep='\t')
# add_back_to_column_positive.rename(columns={'Range': 'Total Base Count'}, inplace=True)
# print(add_back_to_column_positive)
# print(add_back_to_column_positive.dtypes)

# Positive_Strand_Range = Positive_Strand_Range.drop('freq', axis=1)
# Positive_Strand_Range['Total Base Count'] = add_back_to_column_positive['Total Base Count']
# print(Positive_Strand_Range)
# Positive_Strand_Range.to_csv('Positive_Strands_Spanning_Base_Pairs.csv', sep='\t')
#Filtered_Base_Pairs_Of_Positive_Strands.sort_values('Node').to_csv('RangeOfBasePairs_Positive.csv', sep='\t')
#https://stackoverflow.com/questions/28236305/how-do-i-sum-values-in-a-column-that-match-a-given-condition-using-pandas
#Filtered_Base_Pairs_Of_Positive_Strands.to_csv('RangeOfBasePairs_Positive.csv', sep='\t')
# https://stackoverflow.com/questions/49619842/creating-a-range-based-on-two-pandas-columns
# filter this for + and then subtraction operation on basepairend - basepairstart
# https://stackoverflow.com/questions/49619842/creating-a-range-based-on-two-pandas-columns
# order them according to the column order in Positive_Strands - sort via frequency
# https://stackoverflow.com/questions/18265935/how-do-i-create-a-list-with-numbers-between-two-values
# and https://stackoverflow.com/questions/49619842/creating-a-range-based-on-two-pandas-columns
# add column values for rows containing X node name
# separate column, reorder accoridng to first df?
# then join to first df and get range of basepairs.

# Filter Complete_Strand_List according to pos strand values, then print that out to find the pos strands bases
# then add these columns to positive strands csv
# repeat above 2 steps for negative strands.
# 21 FEB WORKS, MAKE OTHERS 100.
# then add % negative strands to Positive_Strands
#% of + per row. , also showing for - strands, meaning 0% + for all those rows
# column for 'spanning bases'
# Above complete..

# 01/03/2024..
# Remaining negative strands..

# make one for list only negative strands, then same filter
# output in text+graphical format
# arrange node start point in ascending order(.. done 'Complete Strand List..)
# cluster them +/- possible TEs
# try on another file..
# overlap graphics..
#TotalStrands_List = TotalStrands_List[~TotalStrands_List['Node'].isin(PlusStrands_List['Node'].values)] #REMOVE VALUES NOT IN PLUS STRANDS LIST
# OR TRY https://stackoverflow.com/questions/56587750/how-do-i-delete-rows-in-one-csv-based-on-another-csv

#total_strands_per_chromosome['PlusStrands'] = plus_strands['Node']
#res = Annot_New.sort_values('Node').groupby('ID')['Color'].apply((list))\
    #.value_counts()
#res = Annot_New.groupby('Node').apply(lambda g: pd.Series(g['Strand'].unique()).str.contains("\+").sum()).reset_index()
#res.to_csv('res.csv', sep='\t')


# count number of times character appears  in each 'unique' row, iterate over all and output to CSV or alike, then..
# see https://stackoverflow.com/questions/57387613/how-to-do-value-counts-based-on-value-in-another-dataframe-column-python
# find % of + per row. # insert the new column 2 from here into the first df, calculate: ***% of forward strand in each chromosome****
# Present in text/graph format, and then add columns for basepair, minus to find across no. base

# on chromosome node (no.) get amount of + strand, - strand.



# + strand covers how many bases..
# dataframe filter row to +, subtract 2 values of base pair difference and add total.







# overlap btw. strands..
#extract fasta seqs at bottom of first file and make new dataframe with header and sequence 2 cols
#match the +/- headers to the 2nd dataframe
# bring the column back with the sequence.
# extract base pair range as coordinates for line plot
# do for all +ve
# 'flip' - to -ve
# draw them
# overlap regions refer back to df and note down
# filter them in original test.gff and then can find the regions
