# 细胞间通讯分析NicheNet

细胞间相互通讯的过程包括**配体---受体---信号蛋白---转录调节子（TF）---靶基因**。现有的软件大部分只能根据ligands-receptors pair计算相互作用的强弱，对于下游的target genes的分析不足。NicheNet通过整合多个数据库，并且训练网络权重，能够预测配体的活性，并且计算上游的ligands的调控潜力。

![image-20211113192822218](https://tva1.sinaimg.cn/large/008i3skNly1gwdqyj9gfoj30ve0peaei.jpg)

NicheNet主要的功能包括：

1. 细胞互作的分析（ligands-targets&ligands-receptor）
2. 配体调控潜力的计算（ligands）

Reference：Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods 17, 159–162 (2020). https://doi.org/10.1038/s41592-019-0667-5
