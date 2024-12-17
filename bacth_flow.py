import heapq
from collections import defaultdict, deque
import networkx as nx
import matplotlib.pyplot as plt
import copy
def ma_search(graph):
    """
    Implements the MA_Search algorithm.
    Args:
        graph: A dictionary where graph[u][v] represents the weight of the edge (u, v).
    Returns:
        An ordering of vertices {v1, ..., vn}.
    """
    nodes = list(graph.keys())
    v1 = nodes[0]
    A = {v1}
    ordering = [v1]

    for _ in range(1, len(nodes)):
        max_node = None
        max_weight = -1
        for x in graph:
            if x not in A:
                weight = sum(graph[z].get(x, 0) for z in A)
                if weight > max_weight:
                    max_node = x
                    max_weight = weight
        ordering.append(max_node)
        A.add(max_node)

    return ordering

def batch_flow(graph):
    """
    Implements the BatchFlow algorithm to compute the maximum s-t flow.
    Args:
        graph: A dictionary where graph[u][v] represents the weight of the edge (u, v).
    Returns:
        A flow function F: a dictionary where F[u][v] is the flow from u to v.
    """
    V = ma_search(graph)
    # print(f'ordering is {V}')
    n = len(V)
    s, t = V[-2], V[-1]
    
    # Initialize data structures
    f = defaultdict(lambda: defaultdict(float))
    # 初始化每一个 f[u][v] 为 0
    for u in graph:
        for v in graph[u]:
            f[u][v] = 0

    C = sum(graph[v].get(t, 0) for v in V[:-1])  # w(A_{n-1}, v_n)
    print(f'min s-t cut is {C}, s is {s}, t is {t}')
    
    ori_graph = copy.deepcopy(graph)
    f[s][t] = graph[s][t] # s,t 有边, 先充满
    graph[s][t] -= f[s][t]
    graph[t][s] -= f[s][t]
    # Push flow to s and t
    for j in range(n - 2-1, -1, -1):
        
        f[s][V[j]] = min(graph[s].get(V[j], 0), C - sum(f[s][vl] for vl in V[j + 1:]))
        if f[s][V[j]] >0 and f[s][V[j]] < ori_graph[s][V[j]]:
            print(f'{V[j]} is the last vertex with s-flow pushed, flow is {f[s][V[j]]}/{ori_graph[s][V[j]]}')
        
        graph[s][V[j]] -= f[s][V[j]]
        graph[V[j]][s] -= f[s][V[j]]
        
        f[V[j]][t] = graph[V[j]].get(t, 0)
        
        graph[V[j]][t] -= f[V[j]][t]
        graph[t][V[j]] -= f[V[j]][t]
    
    # 初始化之后输出 s-flow 和 t-flow 的大小
    print(f's-flow is {sum(f[s][v] for v in f[s])}, t-flow is {sum(f[v][t] for v in f[t])}')

    tmp_C_0 = []
    # Balance flow at intermediate vertices
    for i in range(n - 2-1, 0, -1):
        C_s = sum(f[vl][V[i]] for vl in V[i + 1:])
        C_t = sum(f[V[i]][vl] for vl in V[i + 1:])
        C_0 = abs(C_s - C_t)
        # print(f'i is {i}, V[i] is {V[i]}, C_s is {C_s}, C_t is {C_t}, C_0 is {C_0}')
        tmp_C_0.append([C_0,sum(graph[V[i]][v] for v in V[:i])])
        if C_s >= C_t:
            vis1 = [0]*len(V)
            while C_0 > 0:
                k = -1
                for kp in range(i-1, -1, -1):
                    if graph[V[i]][V[kp]]>0 and graph[V[kp]][V[i]]>0 and vis1[kp]==0:
                        k = kp
                        vis1[k] = 1
                        break
                if k==-1:
                    print(f'No k is found! No vertices can be pushed from {V[i]}! remains flow is {C_0}!')
                    break
                push_flow = min(graph[V[i]][V[k]], C_0)
                if push_flow == 0:
                    print(f'No push flow! C_0 is {C_0}')
                    break
                f[V[i]][V[k]] += push_flow
                C_0 -= push_flow
                
                graph[V[i]][V[k]] -= push_flow
                graph[V[k]][V[i]] -= push_flow
            if C_0!=0:
                print(f'abnormal termination in C_s >= C_t!, flow is not balanced! C_0 is {C_0}')
        else:
            while C_0 > 0:
                k = -1
                for kp in range(i-1, -1, -1):
                    if graph[V[i]][V[kp]]>0 and graph[V[kp]][V[i]]>0:
                        k = kp
                        break
                if k==-1:
                    print(f'No k is found! No vertices can be pushed from {V[i]}!remains flow is {C_0}!')
                    break
                push_flow = min(graph[V[k]][V[i]], C_0)
                if push_flow == 0:
                    print(f'No push flow! what\'s wrong?, C_0 is {C_0}')
                    break
                f[V[k]][V[i]] += push_flow
                C_0 -= push_flow
                graph[V[k]][V[i]] -= push_flow
                graph[V[i]][V[k]] -= push_flow
            if C_0!=0:
                print(f'abnormal termination in C_s < C_t!, flow is not balanced! C_0 is {C_0}')

    print(f'all |C_s-C_t| is {tmp_C_0}')
    N = len(tmp_C_0)
    for i in range(N):
        for j in range(i):
            if tmp_C_0[j][0] >= tmp_C_0[i][0] + sum(graph[V[l]][V[j]] for l in range(i+1)):
                print(f'true for {i}')
                break
    # Construct the flow function
    F = defaultdict(lambda: defaultdict(float))
    for u in f:
        for v in f[u]:
            F[u][v] = f[u][v]
        
    return F

def is_legal_flow(graph, flow, s, t, ordering):
    """
    Checks if the computed flow is legal.
    Args:
        graph: Original graph.
        flow: Computed flow.
        s: Source vertex.
        t: Sink vertex.
    Returns:
        True if the flow is legal, False otherwise.
    """
    for u in graph:
        for v in graph[u]:
            if flow[u][v] > graph[u][v]:
                return False
            if flow[u][v] < 0:
                return False
    illegal_vertices = []

    for u in ordering:
        if u == s:
            continue
        if u == t:
            continue
        inflow = sum(flow[v][u] for v in flow if u in flow[v])
        outflow = sum(flow[u][v] for v in flow[u])
        if inflow != outflow:
            illegal_vertices.append((u,inflow,outflow))

    return illegal_vertices

def visualize_flow(graph, flow):
    """
    Visualizes the flow graph.
    Args:
        graph: Original graph.
        flow: Computed flow.
    """
    G = nx.DiGraph()
    for u in graph:
        for v, capacity in graph[u].items():
            if graph[u][v] > 0:
                G.add_edge(u, v, capacity=capacity, flow=flow[u].get(v, 0))

    pos = nx.spring_layout(G)
    edge_labels = {(u, v): f"{d['flow']}/{d['capacity']}" for u, v, d in G.edges(data=True)}
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=2000, font_size=10, font_weight='bold')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)
    plt.show()

import random

def generate_random_graph(num_nodes, max_weight=10):
    graph = {f'u{i}': {} for i in range(num_nodes)}  # Create nodes 'v0' to 'v9'
    
    # Add random edges with random weights
    for u in graph:
        for v in graph:
            r = random.random()
            if u != v and r < 0.7:  # 30% chance to create an edge
                weight = random.randint(1, max_weight)
                graph[u][v] = weight
                graph[v][u] = weight
            else:
                graph[u][v] = 0
                graph[v][u] = 0
    
    return graph

# Generate a random graph with 10 nodes
graph = generate_random_graph(50)
# 将生成的图写入 data.txt 中, 不覆盖原有的数据, 而是追加
# with open('data.txt', 'a') as f:
#     f.write(str(graph)+'\n')

# 从文件 inbalanced_example.txt 中读取图数据
# with open('inbalanced_example_100.txt', 'r') as f:
#     graph = eval(f.read())

ori_g = copy.deepcopy(graph)

ordering = ma_search(graph)
s, t = ordering[-2], ordering[-1]

# Compute flow using networkx for comparison
G = nx.DiGraph()
for u in graph:
    for v, capacity in graph[u].items():
        G.add_edge(u, v, capacity=capacity)
flow_value, flow_dict = nx.maximum_flow(G, s, t)
print(f"Maximum flow value using networkx: {flow_value}")
# for u in flow_dict:
#     for v in flow_dict[u]:
#         print(f"Flow from {u} to {v}: {flow_dict[u][v]}")

# Compute flow using BatchFlow
import time
start = time.time()
flow = batch_flow(graph)
end = time.time()
print(f"Time taken by BatchFlow: {end - start} seconds")
# for u in flow:
#     for v in flow[u]:
#         print(f"Flow from {u} to {v}: {flow[u][v]}/{ori_g[u][v]}")
print(f'is_legal_flow is {is_legal_flow(ori_g, flow, s, t, ordering)}')
# Compute total flow from s to t
max_flow_value_from_s = sum(flow[s][v] for v in flow[s])
max_flow_value_to_t = sum(flow[v][t] for v in flow[t])
print(f"Maximum flow value using BatchFlow: {max_flow_value_from_s},max_flow_value_to_t is {max_flow_value_to_t,flow[s][t]}")
print(ordering[:10])


# Visualize flows
print("Visualizing BatchFlow result:")
# visualize_flow(ori_g, flow)
# visualize_flow(flow, flow)
# print("Visualizing networkx flow result:")
# visualize_flow(graph, flow_dict)