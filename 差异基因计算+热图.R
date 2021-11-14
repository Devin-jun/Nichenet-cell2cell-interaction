#单细胞数据根据细胞类型计算差异基因的热图

new_m_label <- as.character(stad.m@active.ident)
for(i in 1:length(stad.m@active.ident))
{if (new_m_label[i] %in% c('T0','T2','T4','T6'))
    new_m_label[i] <- 'T'}

stad.m@meta.data['cluster3'] = new_m_label #macrophage分组信息  
DE_table_sender = FindMarkers(object = stad.m, ident.1 = 'T', ident.2 = 'N', min.pct = 0.10, group.by = 'cluster3') %>% rownames_to_column("gene")

#绘制一下这些差异基因的表达情况
mavg_genes <- AverageExpression(stad.m，feature = m_geneset_up)
pheatmap::pheatmap(mavg_genes$RNA,scale='row', cluster_rows=FALSE,
	cluster_col=FALSE,color = colorRampPalette(colors = c("blue","white","red"))(100),cellwidth=20,cellheight=20)
#注意scale的方向