import sys
import os
import pandas as pd
import itertools

positions = sys.argv[1].strip() #used_sites
fasta = sys.argv[2].strip() #final_fasta_file
name_CHROM = sys.argv[3].strip()

#Setting pandas options
pd.options.display.max_rows = 66666
pd.options.display.max_columns = 66666



dpos = pd.read_csv(positions, sep='\t')
dpos ['pos_alignment'] = dpos .index
dpos.drop('NUM_SAMPLES', 1, inplace=True)

sub_chrom = dpos.loc[dpos['#CHROM'] == name_CHROM] ##CHANGE CHROM

start =  sub_chrom["pos_alignment"].iloc[0]
end = sub_chrom["pos_alignment"].iloc[-1] + 1
print(start, end)


dict_seq = {}

with open(fasta) as f:
    for line1,line2 in itertools.zip_longest(*[f]*2):
        main_seq = list(line2.strip())
        seq = main_seq[start : end]
        genome_name = line1.split(">")[1].strip()
        dict_seq[genome_name] = seq
        
dSeq = pd.DataFrame.from_dict(dict_seq)

sub_chrom.drop('pos_alignment', 1, inplace=True)
sub_chrom.reset_index(drop=True, inplace=True)
dSeq.reset_index(drop=True, inplace=True)


df_c = pd.concat([sub_chrom, dSeq], axis=1)


#df_c.index = df_c.index.str.strip()
#df=df_c.reindex(columns = ["#CHROM","POS", "IowaII-ATCC", "CN24716","CN12405","CN24721","CN11707","CN11760","CN9592","CN45011","China","CN42973","CN42632","EG43713","EG34902","EgyptII","EG44493","EgyptI","FIN3","Swe17","UKP124","UKP95","Nor1","Swe7","ITA8","Fra2","Swe5","Human_SWH2","Human_SWH4","UKP107","UKP106","UKP8","Slo4","Cp-29","C6","Slo7","Slo2","C392","ITA4","Pol1","Pol2","C395","C366","ITA3","CZ44619","CZ44621","Cp-13","Swe6","Swe2","Swe1","Hun2","Hun9","Hun1","Hun3","Hun7","Slo5","Slo1","Slo9","Cp-20","Cp-22","Cp-26","Cp-24","Cp-16","Cp-28","Cp-30","Cp-21","Cp-23","ITA10","Cp-27","Cp-17","ITA1","C393","Cp-15","C391","ITA7","Fra1","Fra4","Fra5","Fra6","Fra7","ITA6","ITA12","ITA13","Prt1","UKP4","UKP6","ITA9","C390","C385","Swe13","Swe4","Swe3","Swe9","Swe11","C320","FIN11","FIN24","FIN4","FIN10","FIN6","FIN16","FIN17","FIN15","FIN23","FIN5","FIN22","FIN12","FIN1","FIN2","FIN13","FIN8","UKP90","UKP121","UKP120","UKP94","UKP134","UKP2","UKP125","UKP7","UKP130","UKP104","UKP1","UKP133","UKP3","UKP102","UKP118","UKP103","US41889","US42556","US42557","US42429","US42554","US41565","US42564","US41560","US42561","US41566","US42435","US41883","US41886","US38783","US44519","US41562"])
#print(df.head())
#print(df.tail())


#sub = (df=='CP044421.1').any()
#df_ch = df.loc[:, sub]
#df_ch = df_ch.T


#color based on condition
color_mapping = {'A': 'red', 'a': 'red', 'T': 'blue', 'C': 'green', 'G': 'yellow', 't': 'blue', 'c': 'green', 'g': 'yellow'}
df.style.applymap(lambda v: f"background-color: {color_mapping.get(v, 'black')}").applymap(lambda _: "color: white").to_excel("convertedfile.xlsx", engine='openpyxl')


