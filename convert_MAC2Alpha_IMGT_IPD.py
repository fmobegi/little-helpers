#!/usr/bin/env python3
import pyard
import pandas as pd
import numpy as np

"""
This script takes a single column list of ambiguosly typed alleles in MAC format 
and searches the IPD-IMGT/HLA Database
https://www.ebi.ac.uk/ipd/imgt/hla/ for the MAC=Alpha-numeric conversions
"""
# read the file and extract the desired column as a vector
data = pd.read_csv('search_list_AFDN.txt', header=None, squeeze=True)
list_search = data.values

pyard.max_cache_size = 1_000_000
ard = pyard.ARD()


def process_list(lst):
    results = []
    for item in lst:
        if item[-2:].isupper() and item[-2:] != 'XX' and item[-1:] != 'N':
            print(item)
            result = {"Item": item, "Result": ard.redux_gl(item, 'exon')}
            results.append(result)
    return results


output = process_list(list_search)
df = pd.DataFrame(output)
df = df[['Item', 'Result']]

# write the DataFrame to a file named 'output_file.csv'
df.to_csv('output_file.csv', index=False)
