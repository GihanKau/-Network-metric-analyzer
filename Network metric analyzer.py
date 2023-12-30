'''
Author          : B.M.G.G.K.Rajapaksha
Date            : 13/01/22
Project topic   : A 'Network metric analyzer' by using python and NetworkX package
'''

import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd


class Network:

    # 1.Constructor
    def __init__(self, string_file):
        self.string_file = string_file
        self.p_network = self.create_net()
        self.no_prot = self.prot_counts()
        self.no_inter = self.prot_counts()
        self.network_diameter = self.net_d()

    # 2.Create network
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

    # 3.Method to draw a network graph
    def network_graph(self):

        # Draw the graph
        nx.draw(self.p_network, with_labels=True,
                node_size=200,
                node_color='green',
                font_size=9,
                width=0.2,
                label='PPI Network')

        plt.show()

        # Create a gml file for cytoscape
        nx.write_gml(self.p_network, 'PPI.gml')

    # 4.Method to count nodes and edges
    def prot_counts(self):

        # no_prot = number of proteins, no_inter = number of interactions

        p_network = self.p_network
        no_prot = p_network.number_of_nodes()
        no_inter = p_network.number_of_edges()

        self.no_prot = no_prot
        self.no_inter = no_inter
        return no_inter

    # 5.Method to count the degree of a given protein
    def prot_degree(self, prot_name):
        p_network = self.p_network
        return p_network.degree[prot_name]

    # 6.Method to plot the degree distribution of a PPI network
    def plot_degree_dist(self):
        p_network = self.p_network

        degrees = [p_network.degree(n) for n in p_network.nodes()]

        sns.histplot(degrees, stat='count', binwidth=10)
        plt.xlabel('Degree value')
        plt.title('Degree count vs Degree value')
        plt.show()

    # 7.Method to find the shortest path bw 2 proteins
    def shortest_path(self, prot_1, prot_2):
        p_network = self.p_network
        return nx.shortest_path(p_network, source=prot_1, target=prot_2)

    # 8.Method to find the network diameter of a PPI network
    def net_d(self):

        p_network = self.p_network
        ecc = nx.eccentricity(p_network)
        self.network_diameter = nx.diameter(p_network, ecc)

        return nx.diameter(p_network)

    # 9.Method to calculate the clustering coefficient of a protein
    def cluster_co(self, prot):
        p_network = self.p_network
        cluster_co = round(nx.clustering(p_network, prot), 4)
        return cluster_co

    # 10.Method to calculate betweenness centrality of a given protein
    def betw_cent(self, prot):
        p_network = self.p_network
        betw_cent = round((nx.betweenness_centrality(p_network)[prot]), 5)
        return betw_cent

    # 11.Method to plot the histogram of betweenness centrality
    def betw_cent_hist(self):

        betw_cents = nx.betweenness_centrality(self.p_network)

        cents = [betw_cents[n] for n in betw_cents.keys()]

        sns.histplot(cents, stat='count', bins=50)
        plt.xlabel('Betweenness Centralities')
        plt.title('Histogram of betweenness centralities')
        plt.xlim(0, 1)
        plt.show()

        # return betw_cents

    # 12.Static method to compare two PPI networks
    @staticmethod
    def compare_2(file1, file2):

        pn1 = Network(file1)
        pn2 = Network(file2)

        # Average clustering bar chart
        net1_avg_clust = round(sum(nx.clustering(pn1.p_network).values()) / len(nx.clustering(pn1.p_network).values()),
                               3)
        net2_avg_clust = round(sum(nx.clustering(pn2.p_network).values()) / len(nx.clustering(pn2.p_network).values()),
                               3)

        net1_diam = pn1.network_diameter
        net2_diam = pn2.network_diameter

        data = {'Network': ['Network 1', 'Network 2'],
                'Avg. clustering coefficient': [net1_avg_clust, net2_avg_clust],
                'Network diameter': [net1_diam, net2_diam]}

        df = pd.DataFrame(data, columns=['Network', 'Avg. clustering coefficient', 'Network diameter'])

        fig, axis = plt.subplots(2, 2, )
        fig.suptitle('Comparison of the two PPI networks', fontsize=20)

        # Left hand plot
        axis[0][0].bar(df['Network'], df['Avg. clustering coefficient'], color=['black', 'blue'], width=0.3)
        axis[0][0].set_ylabel('Average clustering coefficient')
        axis[0][0].set_title("Plot of Average clustering coefficient vs Network ")

        # Right hand plot
        axis[0][1].bar(df['Network'], df['Network diameter'], color=['black', 'blue'], width=0.3)
        axis[0][1].set_ylabel('Network diameter')
        axis[0][1].set_title('Plot of Network diameter vs Network')

        # Histogram
        betw_cents1 = nx.betweenness_centrality(pn1.p_network)
        betw_cents2 = nx.betweenness_centrality(pn2.p_network)

        cents1 = [betw_cents1[n] for n in betw_cents1.keys()]
        cents2 = [betw_cents2[n] for n in betw_cents2.keys()]
        # fig, axes = plt.subplots(1, 2)

        # axes[0] = left, axes[1] = right
        sns.histplot(cents1, ax=axis[1][0], color='red', bins=80, binwidth=0.1)
        axis[1][0].set_xlabel('Betweenness centrality')
        axis[1][0].set_title('Betweenness centrality distribution for Network 1')

        sns.histplot(cents2, ax=axis[1][1], color='green')
        axis[1][1].set_xlabel('Betweenness centrality')
        axis[1][1].set_title('Betweenness centrality distribution for Network 2')

        plt.show()

        return df


# ------------------------------------------- Create objects --------------------------------------------------

net1 = Network("string_interactions_short.tsv")

net1.network_graph()


print('No. of proteins in the graph : ', net1.no_prot)
print('No. of interactions in the graph : ', net1.no_inter)
print('Degree of the protein : ', net1.prot_degree('DREB1G'))

net1.plot_degree_dist()

print('Shortest path : ', net1.shortest_path('IRO2', 'PHR3'))
print('Network diameter : ', net1.network_diameter)
print('Cluster coefficient : ', net1.cluster_co('IRO2'))
print('Betweenness centrality : ', net1.betw_cent("IRO2"))

net1.betw_cent_hist()

# Compare the two networks
print('\nComparison of two networks\n', Network.compare_2("string_interactions_short.tsv", "string_interactions_short_2.tsv"))

# END
