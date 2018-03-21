"""analysis of computer networks"""
import urllib2
import random
import time
from matplotlib import pyplot as plt
from collections import deque


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


def copy_graph(graph):
    """
    Make a copy of a graph
    """
    new_graph = {}
    for node in graph:
        new_graph[node] = set(graph[node])
    return new_graph


def delete_node(ugraph, node):
    """
    Delete a node from an undirected graph
    """
    neighbors = ugraph[node]
    ugraph.pop(node)
    for neighbor in neighbors:
        ugraph[neighbor].remove(node)


def targeted_order(ugraph):
    """
    Compute a targeted attack order consisting
    of nodes of maximal degree

    Returns:
    A list of nodes
    """
    # copy the graph
    new_graph = copy_graph(ugraph)

    order = []
    while len(new_graph) > 0:
        max_degree = -1
        for node in new_graph:
            if len(new_graph[node]) > max_degree:
                max_degree = len(new_graph[node])
                max_degree_node = node

        neighbors = new_graph[max_degree_node]
        new_graph.pop(max_degree_node)
        for neighbor in neighbors:
            new_graph[neighbor].remove(max_degree_node)

        order.append(max_degree_node)
    return order

NETWORK_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_rf7.txt"
network_graph = load_graph(NETWORK_URL)

def create_double_sets(inp):
    """creates all possible sets of 2 from input"""
    result = set([])
    for item in inp:
        for other in inp:
            if item != other and (other, item) not in result:
                result.add((item, other))
    return result


def u_er(n, p):
    v = set([num for num in range(n)])
    e = set([])
    g = {}
    for node in v:
        g[node] = set([])
    doubles = create_double_sets(v)
    for double in doubles:
        a = random.random()
        if a < p:
            e.add(double)
    for edge in e:
        g[edge[0]].add(edge[1])
        g[edge[1]].add(edge[0])
    return g


class UPATrial:
    """
    Simple class to encapsulate optimizated trials for the UPA algorithm

    Maintains a list of node numbers with multiple instance of each number.
    The number of instances of each node number are
    in the same proportion as the desired probabilities

    Uses random.choice() to select a node number from this list for each trial.
    """

    def __init__(self, num_nodes):
        """
        Initialize a UPATrial object corresponding to a
        complete graph with num_nodes nodes

        Note the initial list of node numbers has num_nodes copies of
        each node number
        """
        self._num_nodes = num_nodes
        self._node_numbers = [node for node in range(num_nodes) for dummy_idx in range(num_nodes)]

    def run_trial(self, num_nodes):
        """
        Conduct num_nodes trials using by applying random.choice()
        to the list of node numbers

        Updates the list of node numbers so that each node number
        appears in correct ratio

        Returns:
        Set of nodes
        """

        # compute the neighbors for the newly-created node
        new_node_neighbors = set()
        for _ in range(num_nodes):
            new_node_neighbors.add(random.choice(self._node_numbers))

        # update the list of node numbers so that each node number
        # appears in the correct ratio
        self._node_numbers.append(self._num_nodes)
        for dummy_idx in range(len(new_node_neighbors)):
            self._node_numbers.append(self._num_nodes)
        self._node_numbers.extend(list(new_node_neighbors))

        # update the number of nodes
        self._num_nodes += 1
        return new_node_neighbors


def upa(n, m):
    v = set([num for num in range(m)])
    e = set([])
    upa_trial = UPATrial(m)
    g = {}
    for node in v:
        connections = set(v)
        connections.remove(node)
        for connection in connections:
            if (connection, node) not in e:
                e.add((node, connection))
    for i in range(m, n):
        for j in upa_trial.run_trial(m):
            e.add((i, j))
        v.add(i)
    for node in v:
        g[node] = set([])
    for edge in e:
        g[edge[0]].add(edge[1])
        g[edge[1]].add(edge[0])
    return g


def random_order(graph):
    nodes = [key for key in graph]
    result = []
    while len(nodes) > 0:
        result.append(random.choice(nodes))
        nodes.remove(result[-1])
    return result


def bfs_visited(ugraph, start_node):
    """breadth-first search keeping track of what gets visited"""
    que = deque([])
    visited = {start_node}
    que.append(start_node)
    while len(que) > 0:
        node = que.popleft()
        for neighbor in ugraph[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                que.append(neighbor)

    return visited


def cc_visited(ugraph):
    """uses breadth-first search to find sets of connected nodes"""
    remaining = set([node for node in ugraph])
    connected = []
    while len(remaining) > 0:
        node = random.choice(tuple(remaining))
        other = bfs_visited(ugraph, node)
        connected.append(other)
        remaining = remaining.difference(other)
    return connected


def largest_cc_size(ugraph):
    """gives size of largest set of connected nodes in graph"""
    connected = cc_visited(ugraph)
    largest = 0
    for item in connected:
        if len(item) > largest:
            largest = len(item)

    return largest


def compute_resilience(ugraph, attack_order):
    """returns list of largest set of connected nodes in a graph as nodes in attack_order are removed"""
    updated_g = copy_graph(ugraph)
    result = [largest_cc_size(updated_g)]
    for node in attack_order:
        delete_node(updated_g, node)
        result.append(largest_cc_size(updated_g))

    return result


# er_graph = u_er(1239, 0.00397)
# upa_graph = upa(1239, 2)
#
# network_resilience = compute_resilience(network_graph, random_order(network_graph))
# er_resilience = compute_resilience(er_graph, random_order(er_graph))
# upa_resilience = compute_resilience(upa_graph, random_order(upa_graph))
# x_axis = [i for i in range(1240)]
#
# plt.plot(x_axis, network_resilience, 'r', label='Computer Network Graph')
# plt.plot(x_axis, er_resilience, 'b', label='ER Graph (p=0.00397)')
# plt.plot(x_axis, upa_resilience, 'y', label='UPA Graph (m=2)')
# plt.legend(loc='upper right')
# plt.xlabel("Number of Nodes Removed")
# plt.ylabel("Size of Largest Connected Component")
# plt.title("Effect of Node Removal on Largest Connected Component in Various Graphs")
# plt.show()


def fast_targeted_order(g1):
    g = copy_graph(g1)
    degree_sets = []
    for k in range(len(g)):
        degree_sets.insert(k, set([]))
    for node in g:
        d = len(g[node])
        degree_sets[d].add(node)
    lst = []
    i = 0
    for k in reversed(range(len(g))):
        while len(degree_sets[k]) > 0:
            u = random.choice(list(degree_sets[k]))
            degree_sets[k].remove(u)
            for neighbor in g[u]:
                d = len(g[neighbor])
                degree_sets[d].remove(neighbor)
                degree_sets[d - 1].add(neighbor)
            lst.insert(i, u)
            i += 1
            delete_node(g, u)
    return lst

# x_axis = [num for num in range(10, 1000, 10)]
# norm_time = []
# fast_time = []
#
# for n in x_axis:
#     g = upa(n, 5)
#
#     start_time = time.time()
#     targeted_order(g)
#     total_time = time.time() - start_time
#     norm_time.append(total_time)
#
#     start_time = time.time()
#     fast_targeted_order(g)
#     total_time = time.time() - start_time
#     fast_time.append(total_time)
#
# plt.plot(x_axis, norm_time, 'r', label='targeted_order')
# plt.plot(x_axis, fast_time, 'b', label='fast_targeted_order')
# plt.legend(loc='upper left')
# plt.xlabel("Number of Nodes")
# plt.ylabel("Running Time (seconds)")
# plt.title("Running Times of fast_targeted order vs targeted_ order on Desktop Python")
#
# plt.show()


er_graph = u_er(1239, 0.00397)
upa_graph = upa(1239, 2)

network_attack = fast_targeted_order(network_graph)
er_attack = fast_targeted_order(er_graph)
upa_attack = fast_targeted_order(upa_graph)
print "Done computing attack sequences"

network_resilience = compute_resilience(network_graph, network_attack)
print "Done with Network Graph"
er_resilience = compute_resilience(er_graph, er_attack)
print "Done with ER Graph"
upa_resilience = compute_resilience(upa_graph, upa_attack)
print "Done with UPA Graph"

x_axis = [i for i in range(1240)]

plt.plot(x_axis, network_resilience, 'r', label='Computer Network Graph')
plt.plot(x_axis, upa_resilience, 'y', label='UPA Graph (m=2)')
plt.plot(x_axis, er_resilience, 'b', label='ER Graph (p=0.00397)')
plt.legend(loc='upper right')
plt.xlabel("Number of Nodes Removed")
plt.ylabel("Number of Nodes in Largest Connected Group")
plt.title("Effect of Targeted Removal of Nodes on Various Graphs")

plt.show()
