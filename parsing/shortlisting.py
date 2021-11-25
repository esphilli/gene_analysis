import sys
from gene_parsing import GeneParser
from sim_filter import SimParser
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt  
import networkx as nx         

# helper functions for graphs -- G must be a networkx object
def is_subclique(G,nodelist):
    H = G.subgraph(nodelist)
    n = len(nodelist)
    return H.size() == n*(n-1)/2

# change species, timept to compute shortlist for each dataset
species = sys.argv[1]
timept = sys.argv[2]
num_in_list = int(sys.argv[3])
setname = species+'_'+timept+'_top'+str(num_in_list)

if sys.argv[4]=='rep':
    bias1 = float(input('enter bias on Ra: '))
    bias2 = float(input('enter bias on Rn: '))
    g=GeneParser(species+'/'+species+'_'+timept+'.json',num_in_list,bias1,bias2)
    g.parse()
    g.get_stats()
    g.plot_histograms(name=setname)
    g.plot_scatter(name=setname)
    g.filter_rep_score(0.9)
elif sys.argv[4]=='sim':
    s=SimParser(species+'/'+species+'_'+timept+'.json',species+'/'+species+'_'+timept+'_sim.xlsx')
    s.filter_sim_scores(0.9, 0.75)
    s.get_graph(name=setname)
    s.update_stats()
    shortlist = s.no_sim_terms
    for c in s.components:
            if is_subclique(s.graph,c):
                # get term w/highest bit amt
                bit_amts = {}
                for i in c:
                    if i not in bit_amts.keys():
                        bit_amts[i] = s.repeat_stats_df['bits'][i] 
                max_bits = max(bit_amts.values()) 
                max_bits_term = list(key for key, value in bit_amts.items() if value==max_bits)[0]
                #print(max_bits_term)
                shortlist.append(max_bits_term)     
            else:
                # get term w/highest degree
                deg_amts = {}
                degs = nx.degree(s.graph,c)
                for i in degs:
                    deg_amts[i[0]] = i[1]
                max_deg = max(deg_amts.values())
                possible_terms = []
                for key, value in deg_amts.items():
                    if value==max_deg:
                        possible_terms.append(key)
                if len(possible_terms)>1:
                    bit_amts = {}
                    for i in possible_terms:
                        if i not in bit_amts.keys():
                            bit_amts[i] = s.repeat_stats_df['bits'][i]
                    max_bits = max(bit_amts.values())
                    max_deg_term = list(key for key, value in bit_amts.items() if value==max_bits)[0]        
                else:    
                    max_deg_term = possible_terms[0]
                shortlist.append(max_deg_term)
    shortlist_df = s.stats 
    shortlist_df.drop(columns=['Unnamed: 0','level','p-value','fold enrich'],inplace=True)
    for i in s.stats.index:
        if i not in shortlist:
            shortlist_df.drop(index=i,inplace=True)
    with pd.ExcelWriter('~/Desktop/gene_analysis/parsing/'+species+'/'+species+'_'+timept+'_shortlist.xlsx') as writer:    
            shortlist_df.to_excel(writer,sheet_name='shortlist')                  
else:
    print('invalid entry')    




