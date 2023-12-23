import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt

class Network:
    # Attributes
    # p_network = ""


    def __init__(self, string_file):
        self.string_file = string_file
        self.net_diameter = self.net_diameter()
        self.p_network = self.create_net()
        #self.create_net()

    # Create the network

    def create_net(self):

        p_network = nx.Graph()
        file = self.string_file
        with open(file, "r") as string:
            for line in string:
                if '#' not in line:
                    line = line.strip()
                    line = (line.split('\t'))
                    p_network.add_edge(line[0], line[1], weight=float(line[12]))

        self.p_network = p_network
        return p_network


    # Method to draw a network graph
    def network_graph(self):

        # Draw the graph
        nx.draw(self.p_network, with_labels=True)

        plt.show()
        # Create a gml file for cytoscape
        nx.write_gml(self.p_network, 'dreb1a_interactions_Q3_IV.gml')


    # Method to count nodes and edges
    def prot_counts(self):

        p_network=self.p_network
        print('No. of nodes in the graph : ', p_network.number_of_nodes())
        print('No of edges in the graph : ', p_network.number_of_edges())

    def prot_degree(self, prot_name):
        p_network = self.p_network
        print('Degree of', prot_name, 'protein : ', p_network.degree[prot_name])

    def plot_degree_dist(self):
        p_network = self.p_network
        
        degrees = [p_network.degree(n) for n in p_network.nodes()]

        sns.histplot(degrees, stat='count')
        plt.xlabel('Degree')
        plt.title('Degree count vs Degree ')
        plt.show()

    def shortest_path(self, prot_1, prot_2):
        p_network = self.p_network
        print('Shortest path : ', nx.shortest_path(p_network, source=prot_1, target=prot_2))

    def net_diameter(self):
        p_network = self.p_network
        print('Diameter : ', nx.diameter(p_network))
        self.net_diameter = nx.diameter(p_network)
        return nx.diameter(p_network)

    def cluster_co(self, prot):
        p_network = self.p_network
        print('Clustering coefficient : ', round(nx.clustering(p_network, prot), 3))

        # Average clustering
        print('Average clustering coefficient : ', round(sum(nx.clustering(p_network).values())/len(nx.clustering(p_network).values()), 3))




    def betw_cent(self,prot):
        p_network = self.p_network
        print(nx.betweenness_centrality(p_network)[prot]) # not sure

    def betw_cent_hist(self):
        p_network = self.p_network
        print(nx.betweenness_centrality(p_network))  # not sure

        cents = [nx.betweenness_centrality(p_network)[n] for n in nx.betweenness_centrality(p_network).keys()]


        sns.histplot(cents, stat='count', bins=50)
        plt.xlabel('Betweenness Centralities')
        plt.title('The histogram of Betweenness Centralities')
        plt.xlim(0, 1)
        plt.show()

    @staticmethod
    def compare_2(file1, file2):

        pn1 = Network(file1)
        pn2 = Network(file2)

        # Average clustering bar chart
        net1_avg_clust = round(sum(nx.clustering(pn1.p_network).values()) / len(nx.clustering(pn1.p_network).values()), 3)
        net2_avg_clust = round(sum(nx.clustering(pn2.p_network).values()) / len(nx.clustering(pn2.p_network).values()), 3)

        net1_diam = nx.diameter(pn1.p_network)
        net2_diam = nx.diameter(pn2.p_network)

        data = {'Network': ['Net1', 'Net2'],
                'AVG Clust': [net1_avg_clust, net2_avg_clust]}

        network = ['Net1', 'Net2']
        avg_clust = [net1_avg_clust, net2_avg_clust]

        print('Average clustering coefficient p1 : ',
              round(sum(nx.clustering(pn1.p_network).values()) / len(nx.clustering(pn1.p_network).values()), 3),
              '\nAverage clustering coefficient p2 : ',
            round(sum(nx.clustering(pn2.p_network).values()) / len(nx.clustering(pn2.p_network).values()), 3))



s1 = Network("string_interactions_short.tsv")
# s2 = Network("string_interactions_short_2.tsv")
s1.network_graph()
#
# s1.prot_counts()
#
# s1.prot_degree('DREB1G')
# s1.plot_degree_dist()
# s1.shortest_path('IRO2', 'PHR3')
# s1.net_diameter()
# s1.cluster_co('DREB1G')
# s1.betw_cent_hist()
# s1.betw_cent('DREB1G')
# s1.net_diameter

# Network.compare_2("string_interactions_short.tsv", "string_interactions_short_2.tsv")


