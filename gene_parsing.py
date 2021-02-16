from graphviz import Digraph
from abstract_graphs import DAG
from frozendict import frozendict
import sys
import json
#import numpy as np 
import pandas as pd 
#import matplotlib.pyplot as plt 

class GeneParser():
    def __init__(self, file):
        fp = open(file, "r")
        self.jsn = json.loads(fp.read())
        fp.close()
        self.map = {}
        #self.edges = set()  
        self.edges = {} 
        self.filename = file.replace('.json','')
    
    # parse one result tag   
    def parse_gp(self, num, rslt):
        data = self.jsn['overrepresentation']['group']
        result = {}
        edges = []
        if sys.argv[1] != 'all':
            n = int(sys.argv[1])
            gp = data[n]['result']
        else:
            gp = data[num]['result']
        levels = {}
        go = {}
        go_lvl = {}
        if type(gp) == dict:
            if gp['term']['label'] == 'UNCLASSIFIED':
                return result, edges, go, go_lvl
            else:
                result['inpt_0'] = {}
                #result['inpt_0']['list_size'] = gp['input_list']['number_in_list']
                #result['inpt_0']['rep_amt'] = gp['input_list']['number_in_list']/gp['number_in_reference']
                result['inpt_0']['go_id'] = gp['term']['id']
                result['inpt_0']['bio_proc'] = gp['term']['label']
                result['inpt_0']['lvl'] = gp['term']['level']
                result['inpt_0']['gene_list'] = gp['input_list']['mapped_id_list']['mapped_id']
                # get go terms for each input & level
                go['inpt_0'] = gp['term']['id']
                go_lvl[0] = gp['term']['id']
        else:
            for idx, inpt in enumerate(gp):
                result['inpt_'+str(idx)] = {}
                # get labels for each result 
                #result['inpt_'+str(idx)]['list_size'] = inpt['input_list']['number_in_list']
                #result['inpt_'+str(idx)]['rep_amt'] = inpt['input_list']['number_in_list']/inpt['number_in_reference']
                result['inpt_'+str(idx)]['go_id'] = inpt['term']['id']
                result['inpt_'+str(idx)]['bio_proc'] = inpt['term']['label']
                result['inpt_'+str(idx)]['lvl'] = inpt['term']['level']
                result['inpt_'+str(idx)]['gene_list'] = inpt['input_list']['mapped_id_list']['mapped_id']
                # get go terms for each input
                go['inpt_'+str(idx)] = inpt['term']['id']
                # get go terms for each level
                if result['inpt_'+str(idx)]['lvl'] not in go_lvl.keys():
                    go_lvl[result['inpt_'+str(idx)]['lvl']]=[]
                    go_lvl[result['inpt_'+str(idx)]['lvl']].append(result['inpt_'+str(idx)]['go_id'])
                else:
                    if result['inpt_'+str(idx)]['go_id'] not in go_lvl[result['inpt_'+str(idx)]['lvl']]:
                        go_lvl[result['inpt_'+str(idx)]['lvl']].append(result['inpt_'+str(idx)]['go_id'])
                # get info for edges
                if result['inpt_'+str(idx)]['lvl'] not in levels:    
                    levels[idx] = result['inpt_'+str(idx)]['lvl']
            
            # get max and match with next highest and repeat for all
            max_lvl_val = max(levels.values()); {value for key, value in levels.items() if value==max_lvl_val}
            max_lvl_key = max(levels.values()); {key for key, value in levels.items() if value==max_lvl_val}
            for i in range(max_lvl_val):
                for key in list(levels.keys()):
                    if key == max_lvl_key:
                        del levels[key]
                    child_val = max(levels.values()); {value for key, value in levels.items() if value==max_lvl_val}
                    child_key = max(levels.values()); {key for key, value in levels.items() if value==max_lvl_val}
                    if type(child_key) == int:
                        edges.append([max_lvl_key, child_key])
                    else:
                        for child in child_key:    
                            edges.append([max_lvl_key, child])
                    max_lvl_key = child_key
                    max_lvl_val = child_val

            for idx, inpt in enumerate(gp):
                for j in range(len(edges)):
                    parent = edges[j][0]
                    child = edges[j][1]
                    if parent == result['inpt_'+str(idx)]['lvl']:
                        edges[j][0] = result['inpt_'+str(idx)]
                    if child == result['inpt_'+str(idx)]['lvl']:
                        edges[j][1] = result['inpt_'+str(idx)]  

        #print(result)
        #print(edges)
        #print(go)
        return result, edges, go, go_lvl   


    # parse entire file
    def parse_all(self):
        gps = {}
        edge_map = {}
        go_terms = {}
        go_lvl_data = {}
        data = self.jsn['overrepresentation']['group']
        for num, rslt in enumerate(data):
            #print(num)
            result, edges, go, go_lvl = self.parse_gp(num,rslt)
            gps['rslt_'+str(num)] = result
            edge_map['rslt_'+str(num)] = edges
            go_terms['rslt_'+str(num)] = go
            go_lvl_data['rslt_'+str(num)] = go_lvl
        #self.map = gps
        #self.edges = edge_map

        df = pd.DataFrame(gps)
        df = df.T
        df1 = pd.DataFrame(go_terms)
        df1 = df1.T
        df2 = pd.DataFrame(go_lvl_data)
        df2 = df2.T
        with pd.ExcelWriter('~/Desktop/gene_analysis/'+self.filename+'.xlsx') as writer:
            df.to_excel(writer, sheet_name='all data')
            df1.to_excel(writer, sheet_name='go terms for each result')
            df2.to_excel(writer, sheet_name='go terms for each level')
        
        return gps, edge_map, go_terms, go_lvl_data

    # call parsing funcs
    def process(self):   
        if sys.argv[1] == 'all':
            self.parse_all()
        else:
            self.parse_gp(num=None,rslt=None)
        
        # change this stuff
        '''for key in self.map:
            if 'parent_span' in self.map[key].keys():
                parent_key = self.map[key]['parent_span']
                parent = self.map[parent_key]
                pair = (frozendict(parent), frozendict(self.map[key]))
                self.edges.add(pair)
            else:
                self.root = self.map[key]'''

    # eventually add in: 
    #   functions for abstracting down to just the groups to reduce problem
    #   differentiating between functional groups
    #   something about connection between genes up/downregulated & timing of regulation??

    # create DAG
    '''def to_abstract(self):
        dag = DAG()
        for l, r in self.edges:
            dag.add_edge(l, r)
        return dag'''


if __name__=='__main__':
    g=GeneParser('upreg_early_gp1.json')
    g.process()



