import math
import sys
import json
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import networkx as nx  
import pygraphviz as viz           

class GeneParser():

    def __init__(self, file, list_len, bias1, bias2):
        #file = sys.argv[1]
        fp = open(file, "r")
        self.jsn = json.loads(fp.read())
        fp.close()
        self.filename = file.replace('.json','')
        self.list_len = list_len
        self.bias1 = bias1
        self.bias2 = bias2

    # calculate representation score
    def repscore(self, rep_amt, norm_rep_amt):
        #rep_score = 100*(rep_amt+norm_rep_amt)
        #rep_score = -rep_amt**2-norm_rep_amt+1
        rep_score = 100*(rep_amt-self.bias1)*(norm_rep_amt-self.bias2)
        return rep_score

    # parse file
    def parse(self):
        # set up data structures
        mcarr_len = 45101
        edges = []
        lvl_data = {}
        lvl_data['go_id'] = {}
        lvl_data['bio_proc'] = {}
        lvl_data['rep'] = {}
        lvl_data['norm_rep'] = {}
        lvl_data['pval'] = {}
        bits_data = {}
        bits_data['go_id'] = {}
        bits_data['bio_proc'] = {}
        bits_data['rep'] = {}
        bits_data['norm_rep'] = {}
        bits_data['rep_score'] = {}
        bits_data['pval'] = {}
        bits_data['enrch'] = {}
        for val in bits_data.values():
            val['0-2']=[]
            val['2-4']=[]
            val['4-6']=[]
            val['6-8']=[]
            val['8-10']=[]
            val['10-12']=[]
            val['12-14']=[]
        stats = {}
        # parse json file
        data = self.jsn['overrepresentation']['group']
        for num, rslt in enumerate(data):
            data = self.jsn['overrepresentation']['group']
            levels = {}
            if type(data) == dict:
                gp = data['result']
            else:
                gp = data[num]['result']
            if type(gp) == dict:
                if gp['term']['label'] == 'UNCLASSIFIED':
                    continue
                else:
                    if gp['input_list']['number_in_list'] == 0:
                        continue
                    elif gp['input_list']['plus_minus'] == '-':
                        continue
                    else:
                        # get bitwise information content for each GO term
                        prob = gp['number_in_reference']/mcarr_len
                        bt = -math.log(prob,2)
                        # partition by level
                        if 0 not in lvl_data['go_id'].keys():
                            lvl_data['go_id'][0] = []
                        lvl_data['go_id'][0].append(gp['term']['id'])
                        if 0 not in lvl_data['bio_proc'].keys():
                            lvl_data['bio_proc'][0] = []
                        lvl_data['bio_proc'][0].append(gp['term']['label']) 
                        if 0 not in lvl_data['pval'].keys():
                            lvl_data['pval'][0] = []
                        lvl_data['pval'][0].append(gp['input_list']['pValue'])       
                        if 0 not in lvl_data['rep'].keys():
                            lvl_data['rep'][0] = []
                        lvl_data['rep'][0].append(gp['input_list']['number_in_list']/gp['number_in_reference'])
                        if 0 not in lvl_data['norm_rep'].keys():
                            lvl_data['norm_rep'][0] = []
                        lvl_data['norm_rep'][0].append(gp['input_list']['number_in_list']/self.list_len)  
                        # partition by bits 
                        rep_amt = gp['input_list']['number_in_list']/gp['number_in_reference']
                        norm_rep_amt = gp['input_list']['number_in_list']/self.list_len
                        rep_score = self.repscore(rep_amt,norm_rep_amt)
                        if bt>=0 and bt<2:
                            bits_data['go_id']['0-2'].append(gp['term']['id'])
                            bits_data['bio_proc']['0-2'].append(gp['term']['label'])
                            bits_data['rep']['0-2'].append(rep_amt)
                            bits_data['norm_rep']['0-2'].append(norm_rep_amt)
                            bits_data['rep_score']['0-2'].append(rep_score)
                            bits_data['pval']['0-2'].append(gp['input_list']['pValue'])
                            bits_data['enrch']['0-2'].append(gp['input_list']['fold_enrichment'])
                        elif bt>=2 and bt<4:
                            bits_data['go_id']['2-4'].append(gp['term']['id'])
                            bits_data['bio_proc']['2-4'].append(gp['term']['label'])
                            bits_data['rep']['2-4'].append(rep_amt)
                            bits_data['norm_rep']['2-4'].append(norm_rep_amt)
                            bits_data['rep_score']['2-4'].append(rep_score)
                            bits_data['pval']['2-4'].append(gp['input_list']['pValue'])
                            bits_data['enrch']['2-4'].append(gp['input_list']['fold_enrichment'])
                        elif bt>=4 and bt<6:
                            bits_data['go_id']['4-6'].append(gp['term']['id'])
                            bits_data['bio_proc']['4-6'].append(gp['term']['label'])
                            bits_data['rep']['4-6'].append(rep_amt)
                            bits_data['norm_rep']['4-6'].append(norm_rep_amt)
                            bits_data['rep_score']['4-6'].append(rep_score)
                            bits_data['pval']['4-6'].append(gp['input_list']['pValue'])
                            bits_data['enrch']['4-6'].append(gp['input_list']['fold_enrichment'])
                        elif bt>=6 and bt<8:
                            bits_data['go_id']['6-8'].append(gp['term']['id'])
                            bits_data['bio_proc']['6-8'].append(gp['term']['label'])
                            bits_data['rep']['6-8'].append(rep_amt)
                            bits_data['norm_rep']['6-8'].append(norm_rep_amt)
                            bits_data['rep_score']['6-8'].append(rep_score)
                            bits_data['pval']['6-8'].append(gp['input_list']['pValue'])
                            bits_data['enrch']['6-8'].append(gp['input_list']['fold_enrichment'])
                        elif bt>=8 and bt<10:
                            bits_data['go_id']['8-10'].append(gp['term']['id'])
                            bits_data['bio_proc']['8-10'].append(gp['term']['label'])
                            bits_data['rep']['8-10'].append(rep_amt)
                            bits_data['norm_rep']['8-10'].append(norm_rep_amt)
                            bits_data['rep_score']['8-10'].append(rep_score)
                            bits_data['pval']['8-10'].append(gp['input_list']['pValue'])
                            bits_data['enrch']['8-10'].append(gp['input_list']['fold_enrichment'])
                        elif bt>=10 and bt<12: 
                            bits_data['go_id']['10-12'].append(gp['term']['id'])
                            bits_data['bio_proc']['10-12'].append(gp['term']['label'])
                            bits_data['rep']['10-12'].append(rep_amt)
                            bits_data['norm_rep']['10-12'].append(norm_rep_amt)
                            bits_data['rep_score']['10-12'].append(rep_score)
                            bits_data['pval']['10-12'].append(gp['input_list']['pValue'])
                            bits_data['enrch']['10-12'].append(gp['input_list']['fold_enrichment'])  
                        elif bt>=12 and bt<14:
                            bits_data['go_id']['12-14'].append(gp['term']['id'])
                            bits_data['bio_proc']['12-14'].append(gp['term']['label'])
                            bits_data['rep']['12-14'].append(rep_amt)
                            bits_data['norm_rep']['12-14'].append(norm_rep_amt)
                            bits_data['rep_score']['12-14'].append(rep_score)
                            bits_data['pval']['12-14'].append(gp['input_list']['pValue'])
                            bits_data['enrch']['12-14'].append(gp['input_list']['fold_enrichment'])       
                        # gather all data
                        stats[gp['term']['id']] = []
                        stats[gp['term']['id']].append(gp['term']['label'])
                        stats[gp['term']['id']].append(gp['term']['level'])
                        stats[gp['term']['id']].append(bt)
                        stats[gp['term']['id']].append(rep_amt)
                        stats[gp['term']['id']].append(norm_rep_amt)
                        stats[gp['term']['id']].append(rep_score)
                        stats[gp['term']['id']].append(gp['input_list']['pValue'])
                        stats[gp['term']['id']].append(gp['input_list']['fold_enrichment']) 
            else:
                for idx, inpt in enumerate(gp):
                    if inpt['input_list']['number_in_list'] == 0:
                        continue
                    elif inpt['input_list']['plus_minus'] == '-':
                        continue
                    else:
                        # get bitwise information content for each GO term  
                        prob = inpt['number_in_reference']/mcarr_len
                        bt = -math.log(prob,2)
                        rep_amt = inpt['input_list']['number_in_list']/inpt['number_in_reference']
                        norm_rep_amt = inpt['input_list']['number_in_list']/self.list_len
                        rep_score = self.repscore(rep_amt,norm_rep_amt)
                        # partition by level
                        if inpt['term']['level'] not in lvl_data['go_id'].keys():
                            lvl_data['go_id'][inpt['term']['level']]=[]
                        lvl_data['go_id'][inpt['term']['level']].append(inpt['term']['id'])
                        if inpt['term']['level'] not in lvl_data['bio_proc'].keys():
                            lvl_data['bio_proc'][inpt['term']['level']]=[]
                        lvl_data['bio_proc'][inpt['term']['level']].append(inpt['term']['label'])
                        if inpt['term']['level'] not in lvl_data['pval'].keys():
                            lvl_data['pval'][inpt['term']['level']]=[]
                        lvl_data['pval'][inpt['term']['level']].append(inpt['input_list']['pValue'])
                        if inpt['term']['level'] not in lvl_data['rep'].keys():
                            lvl_data['rep'][inpt['term']['level']]=[]
                        lvl_data['rep'][inpt['term']['level']].append(inpt['input_list']['number_in_list']/inpt['number_in_reference'])
                        if inpt['term']['level'] not in lvl_data['norm_rep'].keys():
                            lvl_data['norm_rep'][inpt['term']['level']]=[]
                        lvl_data['norm_rep'][inpt['term']['level']].append(inpt['input_list']['number_in_list']/self.list_len)
                        # partition by bits
                        if bt>=0 and bt<2:
                            bits_data['go_id']['0-2'].append(inpt['term']['id'])
                            bits_data['bio_proc']['0-2'].append(inpt['term']['label'])
                            bits_data['rep']['0-2'].append(rep_amt)
                            bits_data['norm_rep']['0-2'].append(norm_rep_amt)
                            bits_data['rep_score']['0-2'].append(rep_score)
                            bits_data['pval']['0-2'].append(inpt['input_list']['pValue'])
                            bits_data['enrch']['0-2'].append(inpt['input_list']['fold_enrichment'])
                        elif bt>=2 and bt<4:
                            bits_data['go_id']['2-4'].append(inpt['term']['id'])
                            bits_data['bio_proc']['2-4'].append(inpt['term']['label'])
                            bits_data['rep']['2-4'].append(rep_amt)
                            bits_data['norm_rep']['2-4'].append(norm_rep_amt)
                            bits_data['rep_score']['2-4'].append(rep_score)
                            bits_data['pval']['2-4'].append(inpt['input_list']['pValue'])
                            bits_data['enrch']['2-4'].append(inpt['input_list']['fold_enrichment'])
                        elif bt>=4 and bt<6:
                            bits_data['go_id']['4-6'].append(inpt['term']['id'])
                            bits_data['bio_proc']['4-6'].append(inpt['term']['label'])
                            bits_data['rep']['4-6'].append(rep_amt)
                            bits_data['norm_rep']['4-6'].append(norm_rep_amt)
                            bits_data['rep_score']['4-6'].append(rep_score)
                            bits_data['pval']['4-6'].append(inpt['input_list']['pValue'])
                            bits_data['enrch']['4-6'].append(inpt['input_list']['fold_enrichment'])
                        elif bt>=6 and bt<8:
                            bits_data['go_id']['6-8'].append(inpt['term']['id'])
                            bits_data['bio_proc']['6-8'].append(inpt['term']['label'])
                            bits_data['rep']['6-8'].append(rep_amt)
                            bits_data['norm_rep']['6-8'].append(norm_rep_amt)
                            bits_data['rep_score']['6-8'].append(rep_score)
                            bits_data['pval']['6-8'].append(inpt['input_list']['pValue'])
                            bits_data['enrch']['6-8'].append(inpt['input_list']['fold_enrichment'])
                        elif bt>=8 and bt<10:
                            bits_data['go_id']['8-10'].append(inpt['term']['id'])
                            bits_data['bio_proc']['8-10'].append(inpt['term']['label'])
                            bits_data['rep']['8-10'].append(rep_amt)
                            bits_data['norm_rep']['8-10'].append(norm_rep_amt)
                            bits_data['rep_score']['8-10'].append(rep_score)
                            bits_data['pval']['8-10'].append(inpt['input_list']['pValue'])
                            bits_data['enrch']['8-10'].append(inpt['input_list']['fold_enrichment'])
                        elif bt>=10 and bt<12: 
                            bits_data['go_id']['10-12'].append(inpt['term']['id'])
                            bits_data['bio_proc']['10-12'].append(inpt['term']['label'])
                            bits_data['rep']['10-12'].append(rep_amt)
                            bits_data['norm_rep']['10-12'].append(norm_rep_amt)
                            bits_data['rep_score']['10-12'].append(rep_score)
                            bits_data['pval']['10-12'].append(inpt['input_list']['pValue'])
                            bits_data['enrch']['10-12'].append(inpt['input_list']['fold_enrichment'])
                        elif bt>=12 and bt<14:
                            bits_data['go_id']['12-14'].append(inpt['term']['id'])
                            bits_data['bio_proc']['12-14'].append(inpt['term']['label'])
                            bits_data['rep']['12-14'].append(rep_amt)
                            bits_data['norm_rep']['12-14'].append(norm_rep_amt)
                            bits_data['rep_score']['12-14'].append(rep_score)
                            bits_data['pval']['12-14'].append(inpt['input_list']['pValue'])
                            bits_data['enrch']['12-14'].append(inpt['input_list']['fold_enrichment'])
                        # gather all data
                        stats[inpt['term']['id']] = []
                        stats[inpt['term']['id']].append(inpt['term']['label'])
                        stats[inpt['term']['id']].append(inpt['term']['level'])
                        stats[inpt['term']['id']].append(bt)
                        stats[inpt['term']['id']].append(rep_amt)
                        stats[inpt['term']['id']].append(norm_rep_amt)
                        stats[inpt['term']['id']].append(rep_score)
                        stats[inpt['term']['id']].append(inpt['input_list']['pValue'])
                        stats[inpt['term']['id']].append(inpt['input_list']['fold_enrichment'])

        return lvl_data, bits_data, stats

    def get_dag_edges(self):
        max_lvl_val = max(levels.values()); {value for key, value in levels.items() if value==max_lvl_val}
        max_lvl_key = max(levels.values()); {key for key, value in levels.items() if value==max_lvl_val}
        for i in range(max_lvl_val):
            for key in list(levels.keys()):
                if key == max_lvl_key:
                    del levels[key]
                    if len(levels) == 0:
                        break
                child_val = max(levels.values()); {value for key, value in levels.items() if value==max_lvl_val}
                child_key = max(levels.values()); {key for key, value in levels.items() if value==max_lvl_val}
                if type(child_key) == int:
                    edges.append([max_lvl_key, child_key])
                else:
                    for child in child_key:    
                        edges.append([max_lvl_key, child])
                max_lvl_key = child_key
                max_lvl_val = child_val

        # create map of edges according to levels
        for idx, inpt in enumerate(gp):
            for j in range(len(edges)):
                parent = edges[j][0]
                child = edges[j][1]
                if parent == result['inpt_'+str(idx)]['lvl']:
                    edges[j][0] = result['inpt_'+str(idx)]
                if child == result['inpt_'+str(idx)]['lvl']:
                    edges[j][1] = result['inpt_'+str(idx)]
   
        return edges

    # create dataframes
    def get_lvl_bit_dfs(self):
        # get data
        lvl_data, bits_data, stats = self.parse()
        # create dataframes for level partitioning
        lvl_go_df = pd.DataFrame.from_dict(lvl_data['go_id'],orient='index')
        lvl_go_df = lvl_go_df.transpose()
        lvl_bioproc_df = pd.DataFrame.from_dict(lvl_data['bio_proc'],orient='index')
        lvl_bioproc_df = lvl_bioproc_df.transpose()
        lvl_pval_df = pd.DataFrame.from_dict(lvl_data['pval'],orient='index')
        lvl_pval_df = lvl_pval_df.transpose()
        lvl_rep_df = pd.DataFrame.from_dict(lvl_data['rep'],orient='index')
        lvl_rep_df = lvl_rep_df.transpose()
        lvl_normrep_df = pd.DataFrame.from_dict(lvl_data['norm_rep'],orient='index')
        lvl_normrep_df = lvl_normrep_df.transpose()
        # create dataframes for bitwise partitioning
        bits_go_df = pd.DataFrame.from_dict(bits_data['go_id'],orient='index')
        bits_go_df = bits_go_df.transpose()
        bits_bioproc_df = pd.DataFrame.from_dict(bits_data['bio_proc'],orient='index')
        bits_bioproc_df = bits_bioproc_df.transpose()
        bits_pval_df = pd.DataFrame.from_dict(bits_data['pval'],orient='index')
        bits_pval_df = bits_pval_df.transpose()
        bits_rep_df = pd.DataFrame.from_dict(bits_data['rep'],orient='index')
        bits_rep_df = bits_rep_df.transpose()
        bits_normrep_df = pd.DataFrame.from_dict(bits_data['norm_rep'],orient='index')
        bits_normrep_df = bits_normrep_df.transpose()
        bits_repscore_df = pd.DataFrame.from_dict(bits_data['rep_score'],orient='index')
        bits_repscore_df = bits_repscore_df.transpose()

        return lvl_go_df,lvl_bioproc_df,lvl_pval_df,lvl_rep_df,lvl_normrep_df,bits_go_df,bits_bioproc_df,bits_pval_df,bits_rep_df,bits_normrep_df,bits_repscore_df

    def get_stats(self):
       # get data
        lvl_data, bits_data, stats = self.parse()
        # create overall stats dataframe
        stats_df = pd.DataFrame.from_dict(stats,orient='index')
        stats_df = stats_df.rename(columns={0:'process',1:'level',2:'bits',3:'rep amt',4:'norm rep amt',5:'rep score',6:'p-value',7:'fold enrich'})
        
        with pd.ExcelWriter('~/Desktop/gene_analysis/parsing/'+self.filename+'_stats.xlsx') as writer:    
            stats_df.to_excel(writer,sheet_name='stats') 

        self.stats = stats_df
        return self.stats 

    # create excel sheets    
    def get_lvl_bit_sheets(self): 
        # get dataframes
        lvl_go_df,lvl_bioproc_df,lvl_pval_df,lvl_rep_df,lvl_normrep_df,bits_go_df,bits_bioproc_df,bits_pval_df,bits_rep_df,bits_normrep_df,bits_repscore_df = self.get_lvl_bit_dfs()

        # write to excel sheets
        with pd.ExcelWriter('~/Desktop/gene_analysis/parsing/'+self.filename+'_lvl.xlsx') as writer:
            lvl_go_df.to_excel(writer, sheet_name='go terms')
            lvl_bioproc_df.to_excel(writer, sheet_name='bio proc')
            lvl_pval_df.to_excel(writer, sheet_name='p value')
            lvl_rep_df.to_excel(writer, sheet_name='rep amt')
            lvl_normrep_df.to_excel(writer, sheet_name='norm rep amt')

        with pd.ExcelWriter('~/Desktop/gene_analysis/parsing/'+self.filename+'_bits.xlsx') as writer2:
            bits_go_df.to_excel(writer2, sheet_name='go terms')
            bits_bioproc_df.to_excel(writer2, sheet_name='bio proc')
            bits_pval_df.to_excel(writer2, sheet_name='p value')
            bits_rep_df.to_excel(writer2, sheet_name='rep amt')
            bits_normrep_df.to_excel(writer2, sheet_name='norm rep amt')
            bits_repscore_df.to_excel(writer2, sheet_name='rep score') 

        return     
            
    # create histograms
    def plot_histograms(self, name=None):
        # get dataframes
        lvl_go_df,lvl_bioproc_df,lvl_pval_df,lvl_rep_df,lvl_normrep_df,bits_go_df,bits_bioproc_df,bits_pval_df,bits_rep_df,bits_normrep_df,bits_repscore_df = self.get_lvl_bit_dfs()    

        # histograms for bitwise partitioning
        colors = ['red','darkorange','limegreen','turquoise','blue','violet','magenta']
        bits_rep_df = bits_rep_df.dropna(axis=1,how='all')
        bits_normrep_df = bits_normrep_df.dropna(axis=1,how='all')
        bits_repscore_df = bits_repscore_df.dropna(axis=1,how='all')
        numcols = len(bits_repscore_df.columns)
        cols = list(bits_repscore_df.columns)

        fig3,ax3=plt.subplots(numcols,1,sharex=True)
        for idx,i in enumerate(cols):
            ax3[idx].hist(bits_rep_df[i],bins=20,color=colors[idx],label=str(i)+' bits')
            ax3[idx].legend(loc='right',bbox_to_anchor=(1.12,0.5),fontsize=9)
        plt.xlim(xmin=0,xmax=1)
        plt.subplots_adjust(hspace=0.6)
        fig3.supxlabel('representation amount')
        fig3.supylabel('frequency')
        fig3.suptitle('representation amount')
        #plt.show()
        if name is not None:
            plt.savefig(name+'repamt.png')
        else:
            plt.savefig('repamt.png')

        fig4,ax4=plt.subplots(numcols,1,sharex=True)
        for idx,i in enumerate(cols):
            ax4[idx].hist(bits_normrep_df[i],bins=20,color=colors[idx],label=str(i)+' bits')
            ax4[idx].legend(loc='right',bbox_to_anchor=(1.12,0.5),fontsize=9)
        plt.xlim(xmin=0,xmax=1)
        plt.subplots_adjust(hspace=0.6)
        fig4.supxlabel('normalized representation amount')
        fig4.supylabel('frequency')
        fig4.suptitle('normalized representation amount')
        #plt.show()
        if name is not None:
            plt.savefig(name+'normrepamt.png')
        else:
            plt.savefig('normrepamt.png') 

        return   

    # create scatter plots
    def plot_scatter(self, name=None):
        bits = self.stats['bits']
        rep = self.stats['rep amt']
        norm_rep = self.stats['norm rep amt']
        rep_score = self.stats['rep score']
        fig,ax=plt.subplots(3,1,sharex=True)
        ax[0].scatter(bits,rep,c='magenta',label='rep amount')
        ax[0].legend(loc='right',bbox_to_anchor=(1.12,0.5),fontsize=10)
        ax[1].scatter(bits,norm_rep,c='purple',label='norm rep amount')
        ax[1].legend(loc='right',bbox_to_anchor=(1.12,0.5),fontsize=10)
        ax[2].scatter(bits,rep_score,c='blue',label='rep score')
        ax[2].legend(loc='right',bbox_to_anchor=(1.12,0.5),fontsize=10)
        plt.xlim(xmin=0,xmax=15)
        plt.subplots_adjust(hspace=0.5)
        fig.supxlabel('bits')
        fig.supylabel('representation')
        fig.suptitle('representation metrics vs bits of GO terms')
        #plt.show()
        if name is not None:
            plt.savefig(name+'rep_scatter.png')
        else:
            plt.savefig('rep_scatter.png')
        return 

    def filter_rep_score(self, qval):
        t = self.stats.quantile(q=qval)
        rep_th = t['rep score']
        filtered_terms = []
        for i in self.stats.index: 
            if self.stats['rep score'][i]>rep_th:
                filtered_terms.append(i)
            else:
                self.stats.drop(index=i,inplace=True) 
        self.filtered_terms = filtered_terms                      

        with pd.ExcelWriter('~/Desktop/gene_analysis/parsing/'+self.filename+'_filtered.xlsx') as writer:    
            self.stats.to_excel(writer,sheet_name='filtered terms')

        return self.stats, self.filtered_terms

    # don't really need this stuff rn!
    def get_max_rep(self):        
        # get GO term, idx, lvl w/max norm rep amt
        colmax = lvl_normrep_df.max()
        idxcolmax = lvl_normrep_df.idxmax()
        #print(colmax)
        #print(idxcolmax)
        totmax = colmax.max()
        lvltotmax = colmax.idxmax()
        idxtotmax = idxcolmax[lvltotmax]
        # get GO term, idx, lvl w/min p-value
        colmin = lvl_pval_df.min()
        idxcolmin = lvl_pval_df.idxmin()
        #print(colmin)
        #print(idxcolmin)
        totmin = colmin.min()
        lvltotmin = colmin.idxmin()
        idxtotmin = idxcolmin[lvltotmin]

        # get data for GO term w/max norm rep amt
        go_max = {}
        go_max['norm_rep_amt'] = totmax
        go_max['p_val'] = lvl_pval_df.at[idxtotmax,lvltotmax]
        go_max['go_term'] = lvl_go_df.at[idxtotmax,lvltotmax]
        go_max['level'] = lvltotmax
        go_max['bio_proc'] = lvl_bioproc_df.at[idxtotmax,lvltotmax]
        print('data for '+self.filename)
        print('GO term w/max norm rep amt')
        print(go_max)
        # get data for GO term w/min p-value
        go_min = {}
        go_min['norm_rep_amt'] = lvl_normrep_df.at[idxtotmin,lvltotmin]
        go_min['p_val'] = totmin
        go_min['go_term'] = lvl_go_df.at[idxtotmin,lvltotmin]
        go_min['level'] = lvltotmin
        go_min['bio_proc'] = lvl_bioproc_df.at[idxtotmin,lvltotmin]
        print('data for '+self.filename)
        print('GO term w/min p-value')
        print(go_min)   

        return

    def get_innerv_data(self):        
        # look for innervation go terms
        innerv_terms = []
        innerv = {}
        for col in lvl_bioproc_df.columns:
            for idx,item in enumerate(lvl_bioproc_df[col]):
                if item is not None:
                    keywords = ['nerve','nervous','neuron','neuro','neural','innervation','synapse','synaptic','axon','dendrite','ganglion']
                    for word in keywords:
                        if word in item:
                            term = lvl_go_df.at[idx,col]
                            innerv_terms.append(term)
                            r = lvl_normrep_df.at[idx,col]
                            p = lvl_pval_df.at[idx,col]
                            #print(item,p)
                            innerv[idx] = []
                            innerv[idx].append(term)
                            innerv[idx].append(item)
                            innerv[idx].append(r)
                            innerv[idx].append(p)
                            break
        innerv_df = pd.DataFrame.from_dict(innerv,orient='index')
        # write to excel sheet
        with pd.ExcelWriter('~/Desktop/gene_analysis/parsing/'+self.filename+'_innerv.xlsx') as writer4:
            innerv_df.to_excel(writer4, sheet_name='go terms for each level')    

        return    

    # create DAG
    def get_dag(self):
        # get data
        lvl_data, bits_data, stats = self.parse()
        edges = self.get_dag_edges()

        return 


if __name__=='__main__':
    g=GeneParser('mouse/mouse_day1.json',100,0.08,0.01)
    g.parse()
    g.get_stats()
    #g.get_lvl_bit_sheets()
    g.plot_histograms()
    g.plot_scatter()
    g.filter_rep_score(0.9)



