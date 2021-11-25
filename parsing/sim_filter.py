from gene_parsing import GeneParser
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt  
import networkx as nx         

# only filtered GO term data can be a SimParser object!!!
class SimParser():

    def __init__(self, filename, simfile):
        self.filename = filename.replace('.json','')
        data = pd.read_excel(self.filename+'_filtered.xlsx')
        self.stats = data.set_index(keys=data.iloc[:,0])
        self.simfile = simfile 
        data2 = pd.read_excel(self.simfile,usecols=['GO1','GO2','RSS','CAS'])
        df = pd.DataFrame(data2)
        df = df.replace(' ',0)
        self.simdata = df

    def filter_sim_scores(self, q_rss, q_cas):
        # get qval quantile of similarity scores
        tr = self.simdata.quantile(q=q_rss)
        tc = self.simdata.quantile(q=q_cas)
        t_rss = tr['RSS']
        t_cas = tc['CAS']
        repeats = [] 
        repeat_count = {}
        no_sim_terms = []
        # filter and get list of pairs & repeats
        for i in self.simdata.index:
            if self.simdata['RSS'][i]>0.9:
                repeats.append((self.simdata['GO1'][i],self.simdata['GO2'][i]))
                if self.simdata['GO1'][i] not in repeat_count.keys():
                    repeat_count[self.simdata['GO1'][i]] = 1
                else:
                    repeat_count[self.simdata['GO1'][i]] += 1
                if self.simdata['GO2'][i] not in repeat_count.keys():
                    repeat_count[self.simdata['GO2'][i]] = 1
                else:
                    repeat_count[self.simdata['GO2'][i]] += 1
            else:
                if self.simdata['RSS'][i]>t_rss and self.simdata['CAS'][i]>t_cas:
                    repeats.append((self.simdata['GO1'][i],self.simdata['GO2'][i]))
                    if self.simdata['GO1'][i] not in repeat_count.keys():
                        repeat_count[self.simdata['GO1'][i]] = 1
                    else:
                        repeat_count[self.simdata['GO1'][i]] += 1
                    if self.simdata['GO2'][i] not in repeat_count.keys():
                        repeat_count[self.simdata['GO2'][i]] = 1
                    else:
                        repeat_count[self.simdata['GO2'][i]] += 1                                    

        for i in self.stats.index:
            if i not in repeat_count.keys():
                no_sim_terms.append(i)

        self.repeats = repeats
        self.repeat_count = repeat_count
        self.no_sim_terms = no_sim_terms
        return self.repeats, self.repeat_count, self.no_sim_terms 

    def get_graph(self, name=None):
        # get info from similar GO terms  
        graph = nx.Graph()
        graph.add_nodes_from(self.repeat_count.keys())
        graph.add_edges_from(self.repeats)
        pos = nx.spring_layout(graph,k=0.3,iterations=40)
        nx.draw_networkx(graph,pos=pos,with_labels=True,node_size=700,node_color='c',width=1.5,font_size=8,font_weight='bold')
        connects = {}
        for node in nx.nodes(graph):
            connects[node] = list(nx.all_neighbors(graph, node)) 
        if name is not None:
            plt.savefig(name+'similarity_graph.png')
        else:    
            plt.savefig('similarity_graph.png')
        #plt.show() 
        self.components = nx.connected_components(graph)  

        self.graph = graph
        self.connects = connects
        return self.graph, self.connects, self.components

    def update_stats(self):
        # get info from similar GO terms 
        repeat_stats = {}
        for term in self.repeat_count.keys():
            repeat_stats[term] = []
            for t in self.stats.index:
                if t == term:
                    repeat_stats[term].append(self.stats['process'][term])
                    repeat_stats[term].append(self.stats['bits'][term])
                    repeat_stats[term].append(self.stats['rep amt'][term])
                    repeat_stats[term].append(self.stats['norm rep amt'][term])
                    repeat_stats[term].append(self.stats['rep score'][term])
                    repeat_stats[term].append(self.repeat_count[term])
                    repeat_stats[term].append(self.connects[term])
        
        #print(repeat_stats)
        fn = self.simfile.replace('.xlsx','')
        repeat_stats_df = pd.DataFrame.from_dict(repeat_stats,orient='index')
        repeat_stats_df = repeat_stats_df.rename(columns={0:'process',1:'bits',2:'rep amt',3:'norm rep amt',4:'rep score',5:'# similar terms',6:'similar terms'})
        #print(repeat_stats_df)
        with pd.ExcelWriter('~/Desktop/gene_analysis/parsing/'+fn+'_info.xlsx') as writer:    
            repeat_stats_df.to_excel(writer,sheet_name='info')

        self.repeat_stats_df = repeat_stats_df
        return self.repeat_stats_df               


'''if __name__=='__main__':
    s=SimParser('top100.json','top100_filtered_sim.xlsx')
    s.filter_sim_scores(0.9)
    s.get_graph()
    s.update_stats()'''



