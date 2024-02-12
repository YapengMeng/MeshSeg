import argparse
import heapq
import itertools
import os
import numpy as np
import trimesh
import networkx as nx
from collections import deque
from utils import dijkstra, minimum_cut

parser = argparse.ArgumentParser(description="网格分割程序")

# 添加参数
parser.add_argument('--input', '-i', default='input\dino.ply', help='输入待分割的mesh')
parser.add_argument('--output', '-o', default=None, help='输出保存网格')
parser.add_argument('--mode', '-m', choices=['k', '01'], default='01', help='分割模式, K路分割或01分割')
parser.add_argument('--hierarchical', '-H', default=0, type=int, help='是否进行层次化分解,0代表不层次化, 大于0代表层次化分解层数')
parser.add_argument('--epsilon', '-e', default=0.1, type=float, help='超参数epsilon 用于界定模糊区域')
parser.add_argument('--delta', '-d', default=0.5, type=float, help='超参数delta 用于平衡测地线距离和角距离')
args = parser.parse_args()



class MeshSegmenter():
    def __init__(self, args) -> None:
        # 加载超参数
        self.epsilon = args.epsilon
        self.delta = args.delta
        self.mode = args.mode
        # 加载网格
        self.mesh = trimesh.load(args.input)
        
    
    def reloadMesh(self, vertices, faces):
        """
            由vertices和faces直接加载mesh
        """
        self.mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

    def doSegment(self):
        """
            外部调用此函数开始计算分割
        """
        # 计算基本信息
        self.calculateBase()
        # 创建图 带权对偶图
        self.mesh_face_graph = nx.Graph()
        # 为图添加节点
        for index in range(self.num_face):
            self.mesh_face_graph.add_node(index)
        # 计算距离作为图中边的权重
        self.calculateDistance()
        # 计算种子和各个面片到种子的距离
        self.calculate_seed_and_P()
        # 初步判定各个面片的归属
        self.seg_preliminary()
        # 处理模糊区域
        self.seg_fuzzy()
        # 上色
        self.set_color()
        # 显示结果
        self.showMesh()
        return self.seg_result

    def calculateBase(self):
        """
            获得面片列表，计算面片数量、法向量和质心
        """
        # 获得面片列表
        self.faces = self.mesh.faces  # mesh.faces，F个面，每个面为[a,b,c]格式，代表顶点index集合
        self.num_face = len(self.faces)
        print("num mesh ",self.num_face)
        # 获取网格的面片法向量
        self.face_normals = self.mesh.face_normals  # F个面，F个法向量
        # 计算所有面片的质心
        self.face_centroids = self.calculate_face_centroids(self.mesh)

    def calculateDistance(self):
        """计算角距离和测地线距离,保存到带权对偶图的dist属性中,使用Dijkstra 算法计算所有节点对之间的最短路径长度,保存到self.path_lengths中"""
        angle_dist_sum = 0
        geo_dist_sum = 0
        num_edge = 0  # 邻接面片数量

        # 第一次遍历邻接面 计算绝对距离
        for idx, face_adj_idx in enumerate(self.mesh.face_adjacency):  # [a,b] 相邻面的index集合, 数组长度和mesh.face_adjacency_edges相同
            # 获取相邻面片的索引
            face_1, face_2 = face_adj_idx
            # 相邻面片的法线
            n_1 = self.face_normals[face_1]
            n_2 = self.face_normals[face_2]
            # 相邻面片的质心
            c_1 = self.face_centroids[face_1]
            c_2 = self.face_centroids[face_2]
            # 相邻面片公共边的组成顶点
            p_0 = self.mesh.vertices[self.mesh.face_adjacency_edges[idx][0]]
            p_1 = self.mesh.vertices[self.mesh.face_adjacency_edges[idx][1]]

            # 计算面片间的角距离
            angle_dist = self.angle_dist_between_faces(n_1, n_2, c_1, c_2)
            # 计算面片间的测地线距离
            # face_adjacency_edges return:edges – Vertex indices which correspond to face_adjacency 范围：0~num_v-1(18882)
            geo_dist = self.geo_dist_between_faces(c_1, c_2, p_0, p_1)

            # print(angle_dist, geo_dist)

            # 累加
            angle_dist_sum += angle_dist
            geo_dist_sum += geo_dist
            num_edge += 1

            # 将边添加到图中，边的权重分别是angle_dist,geo_dist
            self.mesh_face_graph.add_edge(face_adj_idx[0], face_adj_idx[1], angle_dist=angle_dist, geo_dist=geo_dist)


        angle_dist_ave = angle_dist_sum / num_edge
        geo_dist_ave = geo_dist_sum / num_edge

        # 第二次遍历邻接面 计算归一化距离和最终距离在dist属性中
        for idx, face_adj_idx in enumerate(self.mesh.face_adjacency):
            angle_dist_norm = self.mesh_face_graph[face_adj_idx[0]][face_adj_idx[1]]["angle_dist"] / angle_dist_ave
            geo_dist_norm = self.mesh_face_graph[face_adj_idx[0]][face_adj_idx[1]]["geo_dist"] / geo_dist_ave
            self.mesh_face_graph[face_adj_idx[0]][face_adj_idx[1]]["dist"] = self.delta * geo_dist_norm + (1 - self.delta) * angle_dist_norm
        
        # 使用 Dijkstra 算法计算所有节点对之间的最短路径长度
        # self.path_lengths = dict(nx.all_pairs_dijkstra_path_length(self.mesh_face_graph, weight='dist'))
        self.path_lengths = self.calculate_all_pairs_shortest_path_length(self.mesh_face_graph, weight='dist')
        # 补上自己到自己距离为0
        for node, node_dict in self.path_lengths.items():
            node_dict[node]=0

    
    
    def calculate_all_pairs_shortest_path_length(self, graph, weight='dist'):
        all_pairs_shortest_path_length = {}
        for node in graph.nodes:
            all_pairs_shortest_path_length[node] = dijkstra(graph, node, weight)
        return all_pairs_shortest_path_length

    def calculate_seed_and_P(self):
        if self.mode == "k":  # K路分解种子生成
            K_max = 10 # 先指定Kmax
            seed_list = []
            # 选距离所有面距离之和最小的面为seed_start
            distance_sum_list = []
            for s in range(self.num_face):
                sum_t = 0
                for t in range(self.num_face):
                    sum_t+=self.path_lengths[s][t]
                distance_sum_list.append(sum_t)
            seed_start = distance_sum_list.index(min(distance_sum_list))
            print(f"初始种子为{seed_start}, 距离所有点距离之和为{distance_sum_list[seed_start]}")
            seed_list.append(seed_start)
            # 找后续种子 最大化到之前种子的最小距离
            G = []
            for seed_idx in range(1, K_max):  # 开始找第seed_idx个seed
                distance_min_list = []
                for face_idx in range(self.num_face):  # 遍历每个面片
                    distance_face_to_seeds = []
                    for exist_seed_idx_in_seed_list in range(0, seed_idx): # 遍历之前的种子
                        distance_face_to_seeds.append(self.path_lengths[seed_list[exist_seed_idx_in_seed_list]][face_idx])
                    distance_min_list.append(min(distance_face_to_seeds))
                new_seed_idx = distance_min_list.index(max(distance_min_list))
                seed_list.append(new_seed_idx)
                G.append(distance_min_list[new_seed_idx])
                print(f"新增种子为{new_seed_idx}, 距离之前的种子距离最小值为{distance_min_list[new_seed_idx]}")
            # 选最大化G(k)-G(k+1)的K
            G_t = [G[i] - G[i+1] for i in range(K_max-2)]
            print(G_t)
            K = G_t.index(max(G_t)) + 2
            seed_list = seed_list[:K]
            print(f"最终确定{K}个种子，分别是：", seed_list)

            # 计算P
            P = []
            for face_idx in range(self.num_face):
                distance_reciprocal_list = []  # 到各个种子点的距离倒数列表
                for seed_idx in range(K):
                    distance_reciprocal_list.append( float(1 / (self.path_lengths[seed_list[seed_idx]][face_idx] + 1e-6)))
                distance_reciprocal_sum = sum(distance_reciprocal_list)
                p_list = []
                for seed_idx in range(K):
                    p_list.append(distance_reciprocal_list[seed_idx] / distance_reciprocal_sum)
                P.append(p_list)
            
            self.K = K
            self.seed_list = seed_list
            self.P = P
        
        elif self.mode == "01":  # 01分解种子生成
            # 初始化最远节点对和最大距离, 得到种子
            max_distance = 0
            farthest_nodes = None

            # 遍历路径长度来找到最远的节点对
            for node, lengths in self.path_lengths.items():
                for target_node, distance in lengths.items():
                    if distance > max_distance:
                        max_distance = distance
                        farthest_nodes = (node, target_node)

            # 输出最远的节点对和它们之间的距离
            print(f"最远的两个节点是 {farthest_nodes}，它们之间的距离是 {max_distance}")

            # 初始两个种子
            seed1, seed2 = farthest_nodes
            # 迭代更新种子
            while True:
                # 更新P
                P = []
                sum1_obj = []
                sum2_obj = []
                for face_idx in range(self.num_face):
                    distance1 = self.path_lengths[seed1][face_idx] + 1e-6
                    distance2 = self.path_lengths[seed2][face_idx] + 1e-6
                    p1 = (1/distance1)/(1/distance1 + 1/distance2)
                    p2 = (1/distance2)/(1/distance1 + 1/distance2)
                    P.append((p1, p2))
                # 固定P计算新种子
                for seed_idx in range(self.num_face):
                    sum1 = 0
                    sum2 = 0
                    for face_idx in range(self.num_face):
                        sum1 += P[face_idx][0]*self.path_lengths[seed_idx][face_idx]
                        sum2 += P[face_idx][1]*self.path_lengths[seed_idx][face_idx]
                    sum1_obj.append(sum1)
                    sum2_obj.append(sum2)
                sum1_min = 1e10
                sum2_min = 1e10
                seed1_new = seed1
                seed2_new = seed2
                for idx in range(self.num_face):
                    if sum1_obj[idx]<sum1_min:
                        sum1_min = sum1_obj[idx]
                        seed1_new = idx
                    if sum2_obj[idx]<sum2_min:
                        sum2_min = sum2_obj[idx]
                        seed2_new = idx
                print(f"种子点更新为 {seed1_new}和{seed2_new}")

                if seed1_new == seed2_new:
                    break
                
                if seed1_new == seed1 and seed2_new == seed2:
                    break
                else:
                    seed1 = seed1_new
                    seed2 = seed2_new
            self.seed_list = [seed1, seed2]
            self.P = P
            self.K = 2
    
    def seg_preliminary(self):
        """
            初步判定各个面片的归属
            能清晰分类的面片写入self.seg_result列表中 该列表包含num_seed个列表 存放了每个seed所属面片的索引
            模糊面片写入self.fuzzy_face_idx_dict字典中 该字典以seed_idx作为索引 例如"01"."02"."12" 以face_idx列表作为值
        """
        self.fuzzy_face_idx_dict = {}  # 模糊面片索引,使用字符串"01"."02"."12"来索引
        # self.fuzzy_face_idx = []
        self.seg_result = [[] for _ in range(self.K)]  # 分割结果 包含num_seed个列表 代表每个seed所属面片的索引

        if self.mode == "01" or (self.mode == "k" and self.K == 2):  # 实际分两路
            fuzzy_idx_temp = []
            for idx in range(self.num_face):
                # 清晰面片判定规则，从属于两个区域的概率差足够大
                if abs(self.P[idx][0] - self.P[idx][1]) > 2 * self.epsilon:
                    if self.P[idx][0] > self.P[idx][1]:
                        self.seg_result[0].append(idx)
                    else:
                        self.seg_result[1].append(idx)
                else:
                    fuzzy_idx_temp.append(idx)
                    # self.fuzzy_face_idx.append(idx)
            self.fuzzy_face_idx_dict["01"] = fuzzy_idx_temp
        elif self.K > 2:  # 实际分多路
            # 遍历所有排列
            for seed1, seed2 in itertools.combinations(list(range(self.K)), 2):
                self.fuzzy_face_idx_dict[f"{seed1}{seed2}"] = []
            for idx in range(self.num_face):
                p = self.P[idx]
                # 生成带有原始索引的元组列表
                indexed_list = [(num, index) for index, num in enumerate(p)]
                # 获取最大的三个数及其索引
                top_three = heapq.nlargest(3, indexed_list)  # 按[(数值，索引)]排列

                # 打印结果
                # for rank, (num, index) in enumerate(top_three, start=1):
                #     print(f"{rank}st largest number is {num} at index {index}")

                # 清晰面片判定规则 top1 和 top2之间的间距够大
                if abs(top_three[0][0] - top_three[1][0]) > 2 * self.epsilon:
                    self.seg_result[top_three[0][1]].append(idx)
                else:
                    seed1, seed2 = (top_three[0][1], top_three[1][1]) if top_three[0][1] < top_three[1][1] else (top_three[1][1], top_three[0][1])
                    # seed2 = top_three[0][1] if top_three[0][1] > top_three[1][1] else top_three[1][1]
                    self.fuzzy_face_idx_dict[f"{seed1}{seed2}"].append(idx)

    def seg_fuzzy(self):
        """
            判定模糊面片的归属
        """
        # 两两模糊区域解决模糊
        print(list(self.fuzzy_face_idx_dict.keys()))
        for seed_idx1, seed_idx2 in itertools.combinations(list(range(self.K)), 2):  # 类似0 1
            # print("seg_fuzzy", seed1, seed2)
            seed1 = self.seed_list[seed_idx1]
            seed2 = self.seed_list[seed_idx2]
            # print("real seed face idx:", seed1, seed2)
            self.fuzzy_face_idx = self.fuzzy_face_idx_dict[f"{seed_idx1}{seed_idx2}"]
            # print(self.fuzzy_face_idx)
            if len(self.fuzzy_face_idx)==0:  # 如果两个种子间没有fuzzy区域则跳过
                continue
            # print(self.fuzzy_face_idx)
            fuzzy_face_graph = nx.Graph()
            
            for face_adj_idx in self.mesh.face_adjacency:
                # 获取相邻面片的索引
                face_1, face_2 = face_adj_idx
                if face_1 in self.fuzzy_face_idx and face_2 in self.fuzzy_face_idx:  # 若两个面都在fuzzy区域，则在代价图中添加代价
                    fuzzy_face_graph.add_edge(face_adj_idx[0], face_adj_idx[1], fuzzy_c=1 / (1 + self.mesh_face_graph[face_adj_idx[0]][face_adj_idx[1]]["angle_dist"]))
                elif (face_1 in self.fuzzy_face_idx and face_2 not in self.fuzzy_face_idx) or (face_1 not in self.fuzzy_face_idx and face_2 in self.fuzzy_face_idx):
                    # 若两个面中有一个面在fuzzy区域，先添加这两个面之间的代价，再添加源/汇到清晰归属的面之间的代价为无穷大
                    fuzzy_face_graph.add_edge(face_adj_idx[0], face_adj_idx[1], fuzzy_c=1 / (1 + self.mesh_face_graph[face_adj_idx[0]][face_adj_idx[1]]["angle_dist"]))
                    if face_1 in self.fuzzy_face_idx:
                        fuzzy_idx = face_1
                        clear_idx = face_2
                    else:
                        fuzzy_idx = face_2
                        clear_idx = face_1
                    if self.path_lengths[seed1][clear_idx] < self.path_lengths[seed2][clear_idx]:
                        fuzzy_face_graph.add_edge(seed1, clear_idx, fuzzy_c=float("Inf"))
                    else:
                        fuzzy_face_graph.add_edge(seed2, clear_idx, fuzzy_c=float("Inf"))

            # 计算最小割
            cut_value, partition = nx.minimum_cut(fuzzy_face_graph, seed1, seed2, capacity='fuzzy_c')
            # cut_value, partition = minimum_cut(fuzzy_face_graph, seed1, seed2)
            # 解包分割的两个节点集
            reachable, non_reachable = partition

            for idx in self.fuzzy_face_idx:
                if idx in reachable:
                    self.seg_result[seed_idx1].append(idx)
                else:
                    self.seg_result[seed_idx2].append(idx)

    
    def calculate_face_centroids(self, mesh):
        """
        计算网格中每个面片的质心。
        :param mesh: Trimesh 网格对象。
        :return: 每个面片的质心坐标数组。
        """
        centroids = np.zeros((len(mesh.faces), 3))
        for i, face in enumerate(mesh.faces):
            # 计算每个面片的质心
            centroids[i] = mesh.vertices[face].mean(axis=0)
        return centroids

    def angle_dist_between_faces(self, v1, v2, c1, c2):
        """计算相邻面的角距离,输入两个面的法向量v1, v2, 两个面的质心c1, c2"""
        unit_vector_1 = v1 / np.linalg.norm(v1)
        unit_vector_2 = v2 / np.linalg.norm(v2)
        dot_product = np.dot(unit_vector_1, unit_vector_2)
        isConvex = np.dot(v1, c2-c1)<=0
        # print(isConvex)
        if isConvex:  # 凸的
            eta = self.epsilon
        else:
            eta = 1
        angle_dist = eta*(1 - dot_product)
        return angle_dist

    def geo_dist_between_faces(self, c1, c2, p0, p1):
        """计算相邻面的测地线距离, 输入两个面的质心c1,c2 , 两个面公共边的组成顶点p0,p1"""
        axis = p1 - p0
        va = c1 - p0
        vb = c2 - p0
        len_axis = np.sqrt(np.dot(axis,axis))
        len_a = np.sqrt(np.dot(va,va))
        len_b = np.sqrt(np.dot(vb,vb))
        ang_a = np.arccos(np.dot(va,axis)/(len_a*len_axis))
        ang_b = np.arccos(np.dot(vb,axis)/(len_b*len_axis))
        geo_dist = len_a*len_a + len_b*len_b -2*len_a*len_b*np.cos(ang_a+ang_b)
        
        return geo_dist

    def set_color(self):
        """
            设置颜色 RGB为随机值 不透明度都设为255
        """
        colors = np.zeros((self.num_face, 4))
        color_choice = np.random.randint(0, 255, (self.K, 4))  # 随机设置RGBA
        color_choice[:,3] = 255  # 不透明度都设为255
        for seed_idx in range(self.K):
            for face_index in self.seg_result[seed_idx]:
                colors[face_index] = color_choice[seed_idx]

        # 设置颜色
        self.mesh.visual.face_colors = colors
    
    def showMesh(self):
        self.mesh.show()
    
    def saveMesh(self, path):
        self.mesh.export(path)

if __name__=="__main__":
    seg = MeshSegmenter(args=args)
    seg_result_ori_idx = seg.doSegment()  # 第一次分割
    # single_k = seg.K  # 每次分解的K

    if args.hierarchical > 0:
        # 备份原始顶点和面
        ori_v = seg.mesh.vertices
        ori_f = seg.mesh.faces
        seg_result = seg_result_ori_idx  # 使用原始index标注的分割结果
        

        for hierarchical_idx in range(args.hierarchical):
            new_seg_result = []
            for face_idx in seg_result:  # 分割结果中的面片索引
                # 从原始Mesh中获取对应的面片      
                part_faces = ori_f[face_idx]
                # 将部分面片做一次分割
                seg.reloadMesh(vertices=ori_v, faces=part_faces)
                part_sef_result = seg.doSegment()
                # 将分割结果的idx转化成原始idx放进new_seg_result
                for part_list in part_sef_result:
                    new_seg_result.append([face_idx[i] for i in part_list])
            seg_result = new_seg_result
        # 最终上色展示
        k = len(seg_result)
        print(f"最终分成{k}个部分")
        seg.reloadMesh(vertices=ori_v, faces=ori_f)
        seg.calculateBase()
        seg.K = k
        seg.seg_result = seg_result
        seg.set_color()
        seg.showMesh()
    
    # 保存文件
    if args.output is None:
        input_dir, output_name = os.path.split(args.input)
        output_dir = "./output"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_path = os.path.join(output_dir, output_name.split(".")[-2]+f"_mode_{args.mode}_hierarchical_{args.hierarchical}."+output_name.split(".")[-1])
    else:
        output_path = args.output
    seg.saveMesh(path=output_path)
    print(f"保存文件到{output_path}")
    