## 项目说明 / Project Description
复现论文Hierarchical Mesh Decomposition using Fuzzy Clustering and Cuts (Katz and Tal 2003)中的网格分割算法。  
Reimplement the mesh segmentation algorithm described in the paper "Hierarchical Mesh Decomposition using Fuzzy Clustering and Cuts" by Katz and Tal (2003).

实现两种分解模式（01分解/K路分解）和可选的层次化分解。
Implement two decomposition modes (binary decomposition / K-way decomposition) and an optional hierarchical decomposition.

## 安装步骤 / Installation Steps
Windows平台可以直接导入`.yaml`环境文件  
Windows users can directly import the `.yaml` environment file:

  ```
  conda create -f env.yaml
  ```
其他平台推荐从conda命令构建安装，可以参考以下命令  
For other platforms, it is recommended to build and install from the conda command, refer to the following commands:

  ```
  conda create -n mesh_seg python=3.9
  pip install trimesh
  pip install networkx[default]
  pip install pyglet==1.5
  ```


## 使用说明 / Usage Instructions
2次层次化01分解  
Binary decomposition with 2 levels of hierarchy:

  ```
  python main.py -i input_path -m 01 -H 2
  ```
1次层次化K路分解  
K-way decomposition with 1 level of hierarchy:

  ```
  python main.py -i input_path -m k -H 1
  ```
单次K路分解  
Single K-way decomposition:

  ```
  python main.py -i input_path -m k -H 0
  ```

## 结果概览 Results Overview

我们分析了几种方法：
Our analysis covers several approaches:

- **单次01分解 Single 01 Decomposition**: 展示恐龙模型的初始分段。
  Demonstrates the initial segmentation of the dinosaur model.
- **层次化01分解 Hierarchical 01 Decomposition**: 进一步分解模型以获得更详细的分段。
  Further decomposes the model for more detailed segmentation.
- **K-路分解 K-way Decomposition**: 自动选择K=3，将模型沿其长度分成三个主要部分。
  Automatically selects K=3, dividing the model into three main parts along its length.
- **层次化K-路分解 Hierarchical K-way Decomposition**: 通过计算不同部分的特定K值，推进分解过程，总共得到15个分割部分。
  Advances the decomposition process by calculating specific K values for different parts, resulting in a total of 15 segmented parts.

每种技术的效果都通过显示分割恐龙模型的PNG图片来说明。分割结果清晰地划分了恐龙的解剖学部分，突出了层次方法的精确性。
Each technique's effectiveness is illustrated with PNG images showing the segmented dinosaur model. The segmentation results offer a clear division of the dinosaur's anatomy, highlighting the precision of hierarchical methods.

分析成功地展示了层次和K-方式分解方法在网格分割中的潜力。包括恐龙的前肢和大腿在内的分割部分被清晰地分开，展示了这些技术在详细模型分析中的实际应用。
The analysis successfully demonstrates the potential of hierarchical and K-way decomposition methods in mesh segmentation. The segmented parts, including the dinosaur's forelegs and thighs, are distinctly separated, showcasing the practical application of these techniques in detailed model analysis.

![单次01分解 Single 01 Decomposition]("/assets/Single 01 Decomposition.png")
![层次化01分解 Hierarchical 01 Decomposition]("/assets/Hierarchical 01 Decomposition.png")
![K-路分解 K-way Decomposition]("/assets/K-way Decomposition.png")
![层次化K-路分解 Hierarchical K-way Decomposition]("/assets/Hierarchical K-way Decomposition.png")
