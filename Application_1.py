"""analysis of citation graphs"""
import random
import urllib2
from matplotlib import pyplot as plt

CITATION_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_phys-cite.txt"


def load_graph(graph_url):
    """
    Function that loads a graph given the URL
    for a text representation of the graph

    Returns a dictionary that models a graph
    """
    graph_file = urllib2.urlopen(graph_url)
    graph_text = graph_file.read()
    graph_lines = graph_text.split('\n')
    graph_lines = graph_lines[: -1]

    print "Loaded graph with", len(graph_lines), "nodes"

    answer_graph = {}
    for line in graph_lines:
        neighbors = line.split(' ')
        node = int(neighbors[0])
        answer_graph[node] = set([])
        for neighbor in neighbors[1: -1]:
            answer_graph[node].add(int(neighbor))

    return answer_graph

citation_graph = load_graph(CITATION_URL)


def compute_in_degrees(digraph):
    """function to compute the in-degrees of every node in a digraph given as a dictionary"""
    result = {}
    for node in digraph:
        result[node] = 0
    for node in digraph:
        for connection in digraph[node]:
            result[connection] += 1
    return result


def in_degree_distribution(digraph):
    """function to compute unnormalized in-degree distribution"""
    result = {}
    in_degrees = compute_in_degrees(digraph)
    for key in in_degrees:
        result[in_degrees[key]] = 0
    for key in in_degrees:
        result[in_degrees[key]] += 1
    return result


def normalize(data):
    total = sum(data)
    result = [float(num)/total for num in data]
    return result


def plot_in_degree_distribution(digraph, title, x_axis, y_axis):
    in_degrees = in_degree_distribution(digraph)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.set_xlabel(x_axis)
    ax.set_ylabel(y_axis)
    plt.loglog([key for key in in_degrees], normalize([in_degrees[key] for key in in_degrees]), 'o')

# question 1
# in_degrees = in_degree_distribution(citation_graph)
#
# x_vals = [key for key in in_degrees]
# y_vals = normalize([in_degrees[key] for key in in_degrees])
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title("Log-Log Plot of Distribution of In-Degrees of Graph of Citations")
# ax.set_xlabel("In-Degree")
# ax.set_ylabel("Probability")
# plt.loglog(x_vals, y_vals, 'o', basex=10, basey=10)
# plt.show()

# question 2


def create_double_sets(inp):
    """creates all possible sets of 2 from input"""
    result = set([])
    for item in inp:
        for other in inp:
            if item != other:
                result.add((item, other))
    return result


def di_er(n, p):
    """creates directed graph with n nodes and p chance of having a connection from one node to another. returns a
    dictionary """
    v = set([num for num in range(n)])
    e = set([])
    g = {}
    for node in v:
        g[node] = set([])
    sets = create_double_sets(v)
    for double in sets:
        a = random.random()
        if a < p:
            e.add(double)
    for edge in e:
        g[edge[0]].add(edge[1])
    return g

# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)
# digraph = di_er(3000, .1)
# di_in_degrees = in_degree_distribution(digraph)
# in_degrees = in_degree_distribution(citation_graph)
# # ax1.loglog([key for key in in_degrees],
# #            normalize([in_degrees[key] for key in in_degrees]), 'o')
# # ax2.loglog([key for key in di_in_degrees],
# #            normalize([di_in_degrees[key] for key in di_in_degrees]), 'o')
# plt.loglog([key for key in di_in_degrees],
#            normalize([di_in_degrees[key] for key in di_in_degrees]), 'ro',
#            [key for key in in_degrees], normalize([in_degrees[key] for key in in_degrees]), 'bo')
# plt.show()

#question 3
def avg_out_degree(digraph):
    total = 0
    for node in digraph:
        for connection in digraph[node]:
            total += 1
    result = float(total) / len(digraph)
    return result

# print avg_out_degree(citation_graph)

#question 4
class DPATrial:
    """
    Simple class to encapsulate optimized trials for DPA algorithm

    Maintains a list of node numbers with multiple instances of each number.
    The number of instances of each node number are
    in the same proportion as the desired probabilities

    Uses random.choice() to select a node number from this list for each trial.
    """

    def __init__(self, num_nodes):
        """
        Initialize a DPATrial object corresponding to a
        complete graph with num_nodes nodes

        Note the initial list of node numbers has num_nodes copies of
        each node number
        """
        self._num_nodes = num_nodes
        self._node_numbers = [node for node in range(num_nodes) for dummy_idx in range(num_nodes)]

    def run_trial(self, num_nodes):
        """
        Conduct num_node trials using by applying random.choice()
        to the list of node numbers

        Updates the list of node numbers so that the number of instances of
        each node number is in the same ratio as the desired probabilities

        Returns:
        Set of nodes
        """

        # compute the neighbors for the newly-created node
        new_node_neighbors = set()
        for dummy_idx in range(num_nodes):
            new_node_neighbors.add(random.choice(self._node_numbers))

        # update the list of node numbers so that each node number
        # appears in the correct ratio
        self._node_numbers.append(self._num_nodes)
        self._node_numbers.extend(list(new_node_neighbors))

        # update the number of nodes
        self._num_nodes += 1
        return new_node_neighbors


def dpa(n, m):
    v = set([num for num in range(m)])
    e = set([])
    dpa_trial = DPATrial(m)
    g = {}
    for node in v:
        connections = set(v)
        connections.remove(node)
        for connection in connections:
            e.add((node, connection))
    print "done with complete small graph"
    for i in range(m, n):
        for j in dpa_trial.run_trial(m):
            e.add((i, j))
        v.add(i)
    print "done with dpa stuff"
    for node in v:
        g[node] = set([])
    for edge in e:
        g[edge[0]].add(edge[1])
    return g

# plot_in_degree_distribution(dpa(27770, 13), "Log-Log Plot of In-Degree Distribution of DPA Graph", "In-Degree",
#                             "Probability")

dpa_indegrees = in_degree_distribution(dpa(27770, 13))
indegrees = in_degree_distribution(citation_graph)
plt.loglog([key for key in dpa_indegrees], normalize([dpa_indegrees[key] for key in dpa_indegrees]), 'ro',
           [key for key in indegrees], normalize([indegrees[key] for key in indegrees]), 'bo')
plt.show()
