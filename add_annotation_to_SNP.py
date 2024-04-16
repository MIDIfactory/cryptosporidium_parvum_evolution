"""
Add annotation to SNP
"""

import os
import pandas as pd
import sys



GFF = sys.argv[1].strip()
ChrFile = sys.argv[2].strip()

ReadBin = pd.read_excel(ChrFile)
ToBed = ReadBin[['#CHROM','POS', 'POS']]
ToBed.to_csv("bedded.snp", sep='\t', index=None, header=None)
os.system("cat bedded.snp | bedtools intersect -wb -a IowaII.gff -b stdin | grep exon > resulting_match")
#os.system("sed -i 's/CP044422.1    /CP044422.1     /' resulting_match") #check formatting in the gff and add fix if needed
MatchingGFF = pd.read_csv("resulting_match", sep='\t', names = ['GB', 'type', 'start', 'end', 'w1', 'strand', 'w2', 'feature', 'chrmatch', 'POS', 'fakepos'])
MatchingGFF.reset_index(inplace=True)

splitFeature=MatchingGFF['feature'].str.split(';', expand=True)
MatchingGFF['product'] = splitFeature[6]
MatchingGFF['product'] = MatchingGFF['product'] .str.replace(r'product=', '')
annotated_positions = MatchingGFF[['POS', 'product']]
snpAnno = ReadBin.merge(annotated_positions, on='POS', how='left')
snpAnno.drop_duplicates(subset=['POS'], keep='first', inplace=True)

name = "Annotated" + str(ChrFile)

snpAnno.to_excel(name, engine='openpyxl')

