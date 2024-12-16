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
    print(f'ordering is {V}')
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
    
    # original_graph = defaultdict(lambda: defaultdict(float))
    ori_graph = copy.deepcopy(graph)
    f[s][t] = graph[s][t] # s,t 有边, 先充满
    # residual_graph[s][t] -= f[s][t]
    graph[s][t] -= f[s][t]
    graph[t][s] -= f[s][t]
    # Push flow to s and t
    for j in range(n - 2-1, -1, -1):
        # vj = ordering[j]
        f[s][V[j]] = min(graph[s].get(V[j], 0), C - sum(f[s][vl] for vl in V[j + 1:]))
        # residual_graph[s][V[j]] -= f[s][V[j]]
        # print(f'V[{j}] is {V[j]}, s is {s}, t is {t}')
        # print(f'graph[s] is {graph[s]}')
        # if j==3:
        #     print('uuuuu888888',f[s][V[j]],graph[s][V[j]],graph[V[j]][s])
        graph[s][V[j]] -= f[s][V[j]]
        graph[V[j]][s] -= f[s][V[j]]
        # if j==3:
        #     print('uuuuu888888',f[s][V[j]],graph[s][V[j]],graph[V[j]][s])
        f[V[j]][t] = graph[V[j]].get(t, 0)
        # residual_graph[V[j]][t] -= f[V[j]][t]
        graph[V[j]][t] -= f[V[j]][t]
        graph[t][V[j]] -= f[V[j]][t]
    # xx='u8'
    # yy='u4'
    # print(f'graph[u8][u4]:{graph[xx][yy],graph[yy][xx]}')
    # Balance flow at intermediate vertices
    for i in range(n - 2-1, 0, -1):
        # vi = ordering[i]
        C_s = sum(f[vl][V[i]] for vl in V[i + 1:])
        C_t = sum(f[V[i]][vl] for vl in V[i + 1:])
        C_0 = abs(C_s - C_t)

        if C_s >= C_t:
            while C_0 > 0:
                # print(f'i is {i},graph[{V[i]}] is {graph[V[i]]}')
                # print(f'i is {i},{[kp for kp in range(i) if graph[V[i]][V[kp]]>0 or graph[V[kp]][V[i]]>0]}')
                k = -1
                for kp in range(i-1, -1, -1):
                    print(f'i, kp, graph[{V[i]}][{V[kp]}] is {i, kp, graph[V[i]][V[kp]]}')
                    if graph[V[i]][V[kp]]>0 and graph[V[kp]][V[i]]>0:
                        k = kp
                        print(f'get k is {k}!')
                        break
                # k = max(kp for kp in range(i) if graph[V[i]][V[kp]]>0 or graph[V[kp]][V[i]]>0)
                # print(f'k is {k}')
                # vk = ordering[k]
                # print(f'graph[{V[k]}][{V[i]}] is {graph[V[k]][V[i]]},push_flow is {push_flow}')
                if k==-1:
                    print(f'No k is found! No vertices can be pushed from {V[i]}!')
                    break
                push_flow = min(graph[V[i]][V[k]], C_0)
                if push_flow == 0:
                    print(f'No push flow! what\'s wrong?, C_0 is {C_0}')
                    print(f'graph[{V[i]}] is {graph[V[i]]}, i is {i, V[i]}, k is {k, V[k]}')
                    print(f'graph[{V[i]}][{V[k]}] is {graph[V[i]][V[k]]}')
                    break
                # print(f'push_flow is {push_flow}')
                f[V[i]][V[k]] += push_flow
                C_0 -= push_flow
                
                # print(f'C_0 is {C_0}')
                # residual_graph[V[i]][V[k]] -= push_flow
                graph[V[i]][V[k]] -= push_flow
                graph[V[k]][V[i]] -= push_flow
            if C_0!=0:
                print('abnormal termination!')
        else:
            while C_0 > 0:
                # print(f'i is {i},graph[{V[i]}] is {graph[V[i]]}')
                # print(f'i is {i},{[kp for kp in range(i) if graph[V[i]][V[kp]]>0 or graph[V[kp]][V[i]]>0]}')
                k = -1
                for kp in range(i-1, -1, -1):
                    if graph[V[i]][V[kp]]>0 and graph[V[kp]][V[i]]>0:
                        # print(f'graph[{V[i]}][{V[kp]}] is {graph[V[i]][V[kp]]}')
                        k = kp
                        # break
                # k = max(kp for kp in range(i) if graph[V[i]][V[kp]]>0 or graph[V[kp]][V[i]]>0)
                # vk = ordering[k]
                if k==-1:
                    print(f'No k is found! No vertices can be pushed from {V[i]}!')
                    break
                push_flow = min(graph[V[k]][V[i]], C_0)
                if push_flow == 0:
                    print(f'No push flow! what\'s wrong?, C_0 is {C_0}')
                    break
                # print(f'graph[{V[k]}][{V[i]}] is {graph[V[k]][V[i]]},push_flow is {push_flow}')
                # print(f'push_flow is {push_flow}')
                f[V[k]][V[i]] += push_flow
                C_0 -= push_flow
                # print(f'C_0 is {C_0}')
                graph[V[k]][V[i]] -= push_flow
                graph[V[i]][V[k]] -= push_flow
            if C_0!=0:
                print('abnormal termination!')

    # Construct the flow function
    F = defaultdict(lambda: defaultdict(float))
    for u in f:
        for v in f[u]:
            F[u][v] = f[u][v]
    # Compute total flow from s to t
    max_flow_value_from_s = sum(F[s][v] for v in F[s])
    # print([F[s][v] for v in F[s]],[F[v][t] for v in F[t]])
    # print([v for v in F[t]],ori_graph[t])
    # for u in ori_graph[t]:
    #     print(f'{F[u][t]}/{ori_graph[u][t]}')
    print(len(F[t]),len(ori_graph[t]))
    max_flow_value_to_t = sum(F[v][t] for v in F[t])

    print(f"Maximum flow value using BatchFlow: {max_flow_value_from_s},max_flow_value_to_t is {max_flow_value_to_t,F[s][t]}")
    
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
        print(u, inflow,outflow)
        if inflow != outflow:
            illegal_vertices.append((u,inflow,outflow))
            # return False

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

def generate_random_graph(num_nodes, max_weight=1000):
    graph = {f'u{i}': {} for i in range(num_nodes)}  # Create nodes 'v0' to 'v9'
    
    # Add random edges with random weights
    for u in graph:
        for v in graph:
            r = random.random()
            # print(r)
            if u != v and r < 0.7:  # 30% chance to create an edge
                weight = random.randint(1, max_weight)
                graph[u][v] = weight
                graph[v][u] = weight
            else:
                graph[u][v] = 0
                graph[v][u] = 0
    
    return graph

# Generate a random graph with 10 nodes
graph = generate_random_graph(100)
# print("Generated Random Graph:")
# for u in graph:
#     for v, weight in graph[u].items():
#         print(f"Edge from {u} to {v} with weight {weight}")
# original_graph = defaultdict(lambda: defaultdict(float))
# original_graph.update(graph)
# graph = {'u0': {'u0': 0, 'u1': 0, 'u2': 24, 'u3': 36, 'u4': 0, 'u5': 30, 'u6': 22, 'u7': 99, 'u8': 28, 'u9': 98, 'u10': 0, 'u11': 66, 'u12': 36, 'u13': 6, 'u14': 0, 'u15': 54, 'u16': 0, 'u17': 0, 'u18': 62, 'u19': 31}, 'u1': {'u0': 0, 'u1': 0, 'u2': 0, 'u3': 16, 'u4': 46, 'u5': 59, 'u6': 0, 'u7': 6, 'u8': 57, 'u9': 37, 'u10': 47, 'u11': 96, 'u12': 5, 'u13': 71, 'u14': 12, 'u15': 50, 'u16': 63, 'u17': 20, 'u18': 31, 'u19': 36}, 'u2': {'u0': 24, 'u1': 0, 'u2': 0, 'u3': 0, 'u4': 8, 'u5': 67, 'u6': 35, 'u7': 0, 'u8': 5, 'u9': 24, 'u10': 100, 'u11': 0, 'u12': 62, 'u13': 0, 'u14': 1, 'u15': 0, 'u16': 0, 'u17': 46, 'u18': 36, 'u19': 0}, 'u3': {'u0': 36, 'u1': 16, 'u2': 0, 'u3': 
# 0, 'u4': 9, 'u5': 17, 'u6': 15, 'u7': 81, 'u8': 20, 'u9': 89, 'u10': 95, 'u11': 34, 'u12': 0, 'u13': 22, 'u14': 0, 'u15': 0, 'u16': 94, 'u17': 17, 'u18': 0, 'u19': 87}, 'u4': {'u0': 0, 'u1': 46, 'u2': 8, 'u3': 9, 'u4': 0, 'u5': 38, 'u6': 0, 'u7': 0, 'u8': 0, 'u9': 0, 'u10': 100, 'u11': 22, 'u12': 21, 'u13': 28, 'u14': 0, 'u15': 43, 'u16': 35, 'u17': 0, 'u18': 89, 'u19': 51}, 'u5': 
# {'u0': 30, 'u1': 59, 'u2': 67, 'u3': 17, 'u4': 38, 'u5': 0, 'u6': 0, 'u7': 0, 'u8': 17, 'u9': 97, 'u10': 55, 'u11': 65, 'u12': 0, 'u13': 0, 'u14': 28, 'u15': 0, 'u16': 5, 'u17': 19, 'u18': 84, 'u19': 0}, 'u6': {'u0': 22, 'u1': 0, 'u2': 35, 'u3': 15, 'u4': 0, 'u5': 0, 'u6': 0, 'u7': 9, 'u8': 77, 'u9': 90, 'u10': 40, 'u11': 0, 'u12': 60, 'u13': 60, 'u14': 72, 'u15': 33, 'u16': 42, 'u17': 63, 'u18': 63, 'u19': 0}, 'u7': {'u0': 99, 'u1': 6, 'u2': 0, 'u3': 81, 'u4': 0, 'u5': 0, 'u6': 9, 'u7': 0, 'u8': 99, 'u9': 10, 'u10': 0, 'u11': 0, 'u12': 49, 'u13': 0, 'u14': 47, 'u15': 0, 'u16': 30, 'u17': 43, 'u18': 72, 'u19': 40}, 'u8': {'u0': 28, 'u1': 57, 'u2': 5, 'u3': 20, 'u4': 0, 'u5': 17, 'u6': 77, 'u7': 99, 'u8': 0, 'u9': 85, 'u10': 78, 'u11': 0, 'u12': 0, 'u13': 92, 'u14': 24, 'u15': 78, 'u16': 37, 'u17': 0, 'u18': 54, 'u19': 0}, 'u9': {'u0': 98, 'u1': 37, 'u2': 24, 'u3': 89, 'u4': 0, 'u5': 97, 'u6': 90, 'u7': 10, 'u8': 85, 'u9': 0, 'u10': 46, 'u11': 7, 
# 'u12': 0, 'u13': 82, 'u14': 59, 'u15': 24, 'u16': 5, 'u17': 50, 'u18': 0, 'u19': 0}, 'u10': {'u0': 0, 'u1': 47, 'u2': 100, 'u3': 95, 'u4': 100, 'u5': 55, 'u6': 40, 'u7': 0, 'u8': 78, 'u9': 46, 'u10': 0, 'u11': 44, 'u12': 0, 'u13': 96, 'u14': 67, 'u15': 68, 'u16': 33, 'u17': 62, 'u18': 65, 'u19': 34}, 'u11': {'u0': 66, 'u1': 96, 'u2': 0, 'u3': 34, 'u4': 22, 'u5': 65, 'u6': 0, 'u7': 
# 0, 'u8': 0, 'u9': 7, 'u10': 44, 'u11': 0, 'u12': 29, 'u13': 17, 'u14': 60, 'u15': 50, 'u16': 9, 
# 'u17': 0, 'u18': 40, 'u19': 26}, 'u12': {'u0': 36, 'u1': 5, 'u2': 62, 'u3': 0, 'u4': 21, 'u5': 0, 'u6': 60, 'u7': 49, 'u8': 0, 'u9': 0, 'u10': 0, 'u11': 29, 'u12': 0, 'u13': 31, 'u14': 56, 'u15': 5, 'u16': 0, 'u17': 9, 'u18': 52, 'u19': 51}, 'u13': {'u0': 6, 'u1': 71, 'u2': 0, 'u3': 22, 
# 'u4': 28, 'u5': 0, 'u6': 60, 'u7': 0, 'u8': 92, 'u9': 82, 'u10': 96, 'u11': 17, 'u12': 31, 'u13': 0, 'u14': 0, 'u15': 43, 'u16': 57, 'u17': 97, 'u18': 26, 'u19': 4}, 'u14': {'u0': 0, 'u1': 12, 'u2': 1, 'u3': 0, 'u4': 0, 'u5': 28, 'u6': 72, 'u7': 47, 'u8': 24, 'u9': 59, 'u10': 67, 'u11': 
# 60, 'u12': 56, 'u13': 0, 'u14': 0, 'u15': 78, 'u16': 22, 'u17': 0, 'u18': 53, 'u19': 64}, 'u15': {'u0': 54, 'u1': 50, 'u2': 0, 'u3': 0, 'u4': 43, 'u5': 0, 'u6': 33, 'u7': 0, 'u8': 78, 'u9': 24, 'u10': 68, 'u11': 50, 'u12': 5, 'u13': 43, 'u14': 78, 'u15': 0, 'u16': 1, 'u17': 0, 'u18': 94, 'u19': 0}, 'u16': {'u0': 0, 'u1': 63, 'u2': 0, 'u3': 94, 'u4': 35, 'u5': 5, 'u6': 42, 'u7': 30, 'u8': 37, 'u9': 5, 'u10': 33, 'u11': 9, 'u12': 0, 'u13': 57, 'u14': 22, 'u15': 1, 'u16': 0, 'u17': 10, 'u18': 86, 'u19': 69}, 'u17': {'u0': 0, 'u1': 20, 'u2': 46, 'u3': 17, 'u4': 0, 'u5': 19, 'u6': 63, 'u7': 43, 'u8': 0, 'u9': 50, 'u10': 62, 'u11': 0, 'u12': 9, 'u13': 97, 'u14': 0, 'u15': 0, 'u16': 10, 'u17': 0, 'u18': 12, 'u19': 17}, 'u18': {'u0': 62, 'u1': 31, 'u2': 36, 'u3': 0, 'u4': 89, 'u5': 84, 'u6': 63, 'u7': 72, 'u8': 54, 'u9': 0, 'u10': 65, 'u11': 40, 'u12': 52, 'u13': 26, 'u14': 53, 'u15': 94, 'u16': 86, 'u17': 12, 'u18': 0, 'u19': 13}, 'u19': {'u0': 31, 'u1': 36, 'u2': 0, 'u3': 87, 'u4': 51, 'u5': 0, 'u6': 0, 'u7': 40, 'u8': 0, 'u9': 0, 'u10': 34, 'u11': 26, 'u12': 51, 'u13': 4, 'u14': 64, 'u15': 0, 'u16': 69, 'u17': 17, 'u18': 13, 'u19': 0}} 


ori_g = copy.deepcopy(graph)
print(graph)
# graph['u0']['u1']=999
# graph={'u0': {'u0': 0, 'u1': 9, 'u2': 8, 'u3': 5, 'u4': 0, 'u5': 7, 'u6': 2, 'u7': 10, 'u8': 10, 'u9': 3}, 'u1': {'u0': 2, 'u1': 0, 'u2': 0, 'u3': 0, 'u4': 3, 'u5': 10, 'u6': 4, 'u7': 0, 'u8': 7, 'u9': 0}, 'u2': {'u0': 6, 'u1': 0, 'u2': 0, 'u3': 0, 'u4': 2, 'u5': 8, 'u6': 3, 'u7': 3, 'u8': 5, 
# 'u9': 0}, 'u3': {'u0': 9, 'u1': 2, 'u2': 7, 'u3': 0, 'u4': 0, 'u5': 3, 'u6': 0, 'u7': 6, 'u8': 7, 'u9': 2}, 'u4': {'u0': 0, 'u1': 3, 'u2': 0, 'u3': 0, 'u4': 0, 'u5': 0, 'u6': 10, 'u7': 6, 'u8': 10, 'u9': 3}, 'u5': {'u0': 2, 'u1': 0, 'u2': 8, 'u3': 5, 'u4': 1, 'u5': 0, 'u6': 0, 'u7': 6, 'u8': 0, 'u9': 4}, 'u6': {'u0': 8, 'u1': 0, 'u2': 0, 'u3': 2, 'u4': 5, 'u5': 9, 'u6': 0, 'u7': 0, 'u8': 0, 'u9': 0}, 'u7': {'u0': 7, 'u1': 10, 'u2': 6, 'u3': 9, 'u4': 9, 'u5': 1, 'u6': 9, 'u7': 0, 'u8': 7, 'u9': 10}, 'u8': {'u0': 0, 'u1': 2, 'u2': 0, 'u3': 5, 'u4': 0, 'u5': 6, 'u6': 0, 'u7': 5, 'u8': 0, 'u9': 9}, 'u9': {'u0': 1, 'u1': 7, 'u2': 5, 'u3': 5, 'u4': 7, 'u5': 9, 'u6': 1, 
# 'u7': 9, 'u8': 3, 'u9': 0}}

# Example usage
# graph = {
#     'a': {'b': 3, 'c': 2},
#     'b': {'a': 3, 'c': 1, 'd': 2},
#     'c': {'a': 2, 'b': 1, 'd': 4},
#     'd': {'b': 2, 'c': 4}
# }
# 不存在的边边权赋值为 0
# graph = {
#     'a': {'b': 3, 'c': 2, 'd': 0},
#     'b': {'a': 3, 'c': 1, 'd': 2},
#     'c': {'a': 2, 'b': 1, 'd': 4},
#     'd': {'a': 0, 'b': 2, 'c': 4}
# }

ordering = ma_search(graph)
s, t = ordering[-2], ordering[-1]


# Compute flow using networkx for comparison
G = nx.DiGraph()
for u in graph:
    for v, capacity in graph[u].items():
        G.add_edge(u, v, capacity=capacity)
print(G)
flow_value, flow_dict = nx.maximum_flow(G, s, t)
print(f"Maximum flow value using networkx: {flow_value}")
max_flow_value_to_t = sum(flow_dict[v][t] for v in flow_dict[t])
print(max_flow_value_to_t)
# print("Flow details using networkx:")
# for u in flow_dict:
#     for v in flow_dict[u]:
#         print(f"Flow from {u} to {v}: {flow_dict[u][v]}")

# Compute flow using BatchFlow
flow = batch_flow(graph)
# for u in flow:
#     for v in flow[u]:
#         print(f"Flow from {u} to {v}: {flow[u][v]}/{ori_g[u][v]}")
print(f'is_legal_flow is {is_legal_flow(ori_g, flow, s, t, ordering)}')



# s, t = 'a', 'd'  # Example source and sink


# Visualize flows
# print("Visualizing BatchFlow result:")
# visualize_flow(flow, flow)
# print("Visualizing networkx flow result:")
# visualize_flow(graph, flow_dict)
