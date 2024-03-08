#Data is given here is look like sample sheet for demultiplex to get revese complement of index2 column

#we have two ways

#frist one is we can use module from Biopython called Seq ike below

#if there is no biopython then first we need to install biopython

pip install Biopython

#import pandas and Seq pakcage

import pandas as pd

from Bio.Seq import Seq

# Read csv file here i have used colab so path is look like this user need to change path according to file in

#respective folder

df = pd.read_csv('/content/problem_sheet.csv', header=None)

#To check the position of index2

df.head(30)

# Reverse complement of index 2 which is present in column 9 after row 19

df.iloc[20:, 9] = df.iloc[20:, 9].apply(lambda seq: str(Seq(seq).reverse_complement()))

# Save the modified DataFrame to a new CSV file

df.to_csv('/content/problem_sheet_reverse.csv', index=False)

#The second way is we can define our on funcation to do the reverse complement.

#import pandas package

import pandas as pd

# Read csv file here i have used colab so path is look like this user need to change path according to file in

#respective folder

df = pd.read_csv('/content/problem_sheet.csv', header=None)

# user define a function to get the complement

def get_complement(base):

return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get(base, base)

# Apply the function to each value in column 9 for rows after row 19

df.iloc[20:, 9] = df.iloc[20:, 9].apply(lambda sequence: ''.join(map(get_complement,

reversed(sequence))))

# Save the modified DataFrame to a new CSV file

df.to_csv('/content/problem_sheet_modified.csv', index=False)
