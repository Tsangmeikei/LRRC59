

# 读取雷达图数据文件
df = read.csv("TMB.csv",header = TRUE,row.names="TMB.score",sep = ",")


# 绘图
ggradar(df,
        grid.max = max(df[,-1]),                 
        grid.mid = max(df[,-1])/2,               
        grid.min = 0,                           
        grid.label.size = 4,                     
        axis.label.size = 5,                     
        group.colours = rainbow(length(df[,1])), 
        background.circle.colour = "white",      
        group.point.size = 2,                    
        group.line.width = 2,                    
        plot.legend = T,                        
        legend.position = "right",              
        legend.title = "",                      
        legend.text.size = 10,                   
        plot.title   = "Title",                 
        plot.extent.x.sf = 1.2,                  
        plot.extent.y.sf = 1.2)