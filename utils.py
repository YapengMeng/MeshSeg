from collections import deque
import heapq
import networkx as nx
import numpy as np

def dijkstra(graph, start, weight='dist'):
    # 字典 每个节点的最短路径长度 初始时，所有节点的最短路径长度被设为无穷大。
    shortest_paths = {vertex: float('infinity') for vertex in graph.nodes}
    # 起始节点的最短路径长度为0
    shortest_paths[start] = 0
    # 优先队列 （当前最短路径长度，节点）
    priority_queue = [(0, start)]
    
    while priority_queue:
        # 每次弹出当前最短路径的节点
        current_distance, current_vertex = heapq.heappop(priority_queue)
        # 如果当前节点的路径长度大于已知的最短路径长度，则跳过
        if current_distance > shortest_paths[current_vertex]:
            continue
        # 遍历当前节点的相邻节点
        for neighbor in graph.neighbors(current_vertex):
            distance = current_distance + graph[current_vertex][neighbor][weight]
            # 如果计算出的距离小于已知的最短路径，则更新最短路径并将其加入优先队列
            if distance < shortest_paths[neighbor]:
                shortest_paths[neighbor] = distance
                heapq.heappush(priority_queue, (distance, neighbor))
    
    return shortest_paths


# 创建图节点的索引映射
def create_index_mapping(graph):
    mapping = {}  # 节点到索引的映射
    reverse_mapping = {}  # 索引到节点的映射
    for i, node in enumerate(graph.nodes()):
        mapping[node] = i
        reverse_mapping[i] = node
    return mapping, reverse_mapping

# 广度优先搜索（BFS），用于Edmonds-Karp算法中寻找增广路径
def bfs(rGraph, s, t, parent, n):
    visited = [False] * n  # 标记节点是否被访问过
    queue = []  # BFS的队列

    queue.append(s)
    visited[s] = True

    while queue:
        u = queue.pop(0)

        for ind in range(n):
            # 寻找从u到ind的可行边
            if visited[ind] == False and rGraph[u][ind] > 0:
                queue.append(ind)
                visited[ind] = True
                parent[ind] = u

    return True if visited[t] else False

# Edmonds-Karp算法实现最大流
def edmonds_karp(graph, source, sink, mapping):
    n = len(graph)
    rGraph = np.zeros((n, n))  # 剩余网络

    # 初始化剩余网络
    for u, v, data in graph.edges(data=True):
        rGraph[mapping[u]][mapping[v]] = data['fuzzy_c']

    parent = [-1] * n
    max_flow = 0

    # 当存在增广路径时，计算最大流
    while bfs(rGraph, mapping[source], mapping[sink], parent, n):
        path_flow = float("Inf")
        s = mapping[sink]

        # 计算增广路径上的最小残留容量
        while s != mapping[source]:
            path_flow = min(path_flow, rGraph[parent[s]][s])
            s = parent[s]

        max_flow += path_flow

        # 更新剩余网络
        v = mapping[sink]
        while v != mapping[source]:
            u = parent[v]
            rGraph[u][v] -= path_flow
            rGraph[v][u] += path_flow
            v = parent[v]

    return max_flow

# 计算最小割
def minimum_cut(graph, source, sink):
    mapping, reverse_mapping = create_index_mapping(graph)
    max_flow = edmonds_karp(graph, source, sink, mapping)

    # 找到残余网络中可达源点的节点
    residual_graph = nx.adjacency_matrix(graph, weight='fuzzy_c').toarray().astype(int)
    visited = [False] * len(graph)
    bfs(residual_graph, mapping[source], mapping[sink], visited, len(graph))

    reachable = [reverse_mapping[i] for i in range(len(visited)) if visited[i]]
    non_reachable = [reverse_mapping[i] for i in range(len(visited)) if not visited[i]]

    return max_flow, (reachable, non_reachable)

