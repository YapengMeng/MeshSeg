## 项目说明
复现论文Hierarchical Mesh Decomposition using Fuzzy Clustering and Cuts (Katz and Tal 2003)中的网格分割算法。  
实现两种分解模式（01分解/K路分解）和可选的层次化分解。

## 安装步骤
Windows平台可以直接导入`.yaml`环境文件
  ```
  conda create -f env.yaml
  ```
其他平台推荐从conda命令构建安装，可以参考以下命令
  ```
  conda create -n mesh_seg python=3.9
  pip install trimesh
  pip install networkx[default]
  pip install pyglet==1.5
  ```

## 使用说明
2次层次化01分解
  ```
  python main.py -i input_path -m 01 -H 2
  ```
1次层次化K路分解
  ```
  python main.py -i input_path -m k -H 1
  ```
单次K路分解
  ```
  python main.py -i input_path -m k -H 0
  ```