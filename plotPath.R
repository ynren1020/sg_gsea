#' plot enriched pathways 
#' A function to plot enriched pathways, y axis should be normalized enrichment 
#' score, x axis should rank by NES.
#' @param 
#' 
#' @details
#' 
#' @import  
#' 
#' @return 
#' 
#' @example 
#' 
#' @export
#' 
#' 

library(ggplot2)
library(RColorBrewer)
library(ggrepel)

pre_plot <- function(dat){
    dat <- read.delim(dat, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    dat <- dat[, c("pathway","NES","padj")]
    dat$RANK <- order(dat$NES)
    dat$logP <- -log10(dat$padj + 0.00001)
    dat$logP.sign <- sign(dat$NES)*dat$logP
    return(dat)
    
}

dat <- pre_plot("T025160_OV_cor.txt")



# NES vs log10FDR
#"T025160_OV_cor.txt"

plotPath <- function(gsea.df){
    
    dat <- pre_plot(gsea.df)

ptest<-ggplot(dat,aes(x=logP.sign,y=NES,label = ifelse(logP >1.3,as.character(pathway),'')))+
    geom_point(aes(color = logP))+
    #geom_text(aes(label = ifelse(color=="dark",as.character(NAME),'')), vjust= -1)+
    #ggplot2::geom_point(aes(color=ifelse(FDR.q.val < 0.001, "FDR < = 0.001", ifelse(FDR.q.val < 0.01, "FDR < 0.01", ifelse(FDR.q.val ))))+ #mapping aes(color),here define different color based on isSuper!!!
    #ggplot2::scale_colour_manual(name='',values=c('FDR > 0.05'='grey','FDR < = 0.05'='red'))+ #add color manually!!!!
    scale_color_gradient(low = "blue", high = "red") +
    #geom_vline(xintercept = 25982,linetype="dashed",size=0.5,alpha=0.2)+
    #geom_hline(yintercept = 1.3,linetype="dashed",size=0.5,alpha=0.2)+
    #geom_hline(yintercept = 2,linetype="dashed",size=0.5,alpha=0.2)+
    labs(x="sign(NES)*(-log10(FDR))",y="Normal Enrichement Score",color='-log10(q)')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = c(0.85, 0.2),
          text = element_text(size=15))

set.seed(00)
ppd<-ptest+geom_text_repel(data=dat,nudge_x=-1,direction="y",force=2,max.iter=4000,hjust=1,segment.size=0.2,size = 2)

ppd
}
#ggsave("T025160_OV_cor_NES.pdf", width = 6, height = 6)

test <- plotPath("T025160_OV_cor.txt")
