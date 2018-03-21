"""first project for algorithmic thinking"""

EX_GRAPH0 = {0: {1, 2}, 1: set([]), 2: set([])}
EX_GRAPH1 = {0: {1, 4, 5}, 1: {2, 6}, 2: {3}, 3: {0}, 4: {1}, 5: {2}, 6: set([])}
EX_GRAPH2 = {0: {1, 4, 5}, 1: {2, 6}, 2: {3, 7}, 3: {7}, 4: {1}, 5: {2}, 6: set([]), 7: {3}, 8: {1, 2},
             9: {0, 3, 4, 5, 6, 7}}


def make_complete_graph(num_nodes):
    """function to create a full directional graph with num_nodes nodes"""
    result = {}
    for num in range(num_nodes):
        connections = set(range(num_nodes))
        connections.remove((num))
        result[num] = connections
    return result


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
