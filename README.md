# 细胞间通讯分析NicheNet

细胞间相互通讯的过程包括**配体---受体---信号蛋白---转录调节子（TF）---靶基因**。现有的软件大部分只能根据ligands-receptors pair计算相互作用的强弱，对于下游的target genes的分析不足。NicheNet通过整合多个数据库，并且训练网络权重，能够预测配体的活性，并且计算上游的ligands的调控潜力。

![preview](https://tva1.sinaimg.cn/large/008i3skNly1gwdqd283d9j310j0u0q87.jpg)

NicheNet主要的功能包括：

1. 细胞互作的分析（ligands-targets&ligands-receptor）
2. 配体调控潜力的计算（ligands）


