"""implementation of project 2 for algorithmic thinking"""
import random
from collections import deque


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
    updated_g = ugraph
    result = [largest_cc_size(updated_g)]
    for node in attack_order:
        updated_g.pop(node)
        for thing in updated_g:
            if node in updated_g[thing]:
                updated_g[thing].remove(node)
        result.append(largest_cc_size(updated_g))

    return result
