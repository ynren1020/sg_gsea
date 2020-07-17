#' plot enriched pathways 
#' A function to plot enriched pathways, y axis should be normalized enrichment 
#' score, x axis should rank by NES.
#' @param gsea.df A data frame, output from sg_gsea function
#' @param pathway.list Pathways user provided to label, a character variable
#' @param direction Specify pos or neg or both enriched pathway to be labelled
#' @param n a number defining how many pos or neg or both pathways to label
#' @param ... other parameters to adjust plot 
#' 
#' @details
#' 
#' @import ggplot2
#' @import RColorBrewer 
#' @import ggrepel
#'     
#' @return A pdf plot with NES as y axis and sign(NES)*(-log10(FDR)) as x axis
#' 
#' @example 
#' 
#' test <- plotPath("T025160_OV_cor.txt")
#' test <- plotPath("T025160_OV_cor.txt", direction = "user", 
#' pathway.list = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
#' test <- plotPath("ENST00000417354_OV_cor_hallmark.txt")
#' 
#' @export
#' 
#' 


library(ggplot2)
library(RColorBrewer)
library(ggrepel)


plotPath <- function(gsea.df,pathway.list=NULL,direction="both",n=3,...){
    
    # process data for plot ---
        dat <- read.delim(gsea.df, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        dat <- dat[, c("pathway","NES","padj")]
        dat$RANK <- order(dat$NES)
        dat$logP <- -log10(dat$padj + 0.00001)
        dat$logP.sign <- sign(dat$NES)*dat$logP
      
    # pathway to label ---
    if (is.null(pathway.list)) {
        if (direction == "both"){
          label.max <- slice_max(dat, order_by = NES, n=n)
          label.min <- slice_min(dat, order_by = NES, n=n)
          label <- rbind(label.max, label.min)
        
        } else if (direction == "pos"){
               label <- slice_max(dat, order_by = NES, n=n)
        } else {  
               label <- slice_min(dat, order_by = NES, n=n)
         } 
    } else {
        print("Label pathways user provided.")
        label <- data.frame(pathway=pathway.list)
    }
   # label top or bottom pathways --- 
    for (i in 1:nrow(dat)){
        dat$name[i] <- ifelse(dat$pathway[i]%in%label$pathway, "Yes","No")
    }
    
    
ptest<-ggplot(dat,aes(x=logP.sign,y=NES,label = ifelse(name == "Yes",as.character(pathway),"")))+
    geom_point(aes(color = logP))+
    #geom_text(aes(label = ifelse(color=="dark",as.character(NAME),'')), vjust= -1)+
    #ggplot2::geom_point(aes(color=ifelse(FDR.q.val < 0.001, "FDR < = 0.001", ifelse(FDR.q.val < 0.01, "FDR < 0.01", ifelse(FDR.q.val ))))+ #mapping aes(color),here define different color based on isSuper!!!
    #ggplot2::scale_colour_manual(name='',values=c('FDR > 0.05'='grey','FDR < = 0.05'='red'))+ #add color manually!!!!
    scale_color_gradient(low = "blue", high = "red") +
    geom_vline(xintercept = 1.30103,linetype="dashed",size=0.5,alpha=0.2)+
    geom_vline(xintercept = -1.30103,linetype="dashed",size=0.5,alpha=0.2)+
    geom_hline(yintercept = 0,linetype="dashed",size=0.5,alpha=0.2)+
    labs(x="sign(NES)*(-log10(FDR))",y="Normal Enrichement Score",color='-log10(q)')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = c(0.89, 0.25),
          text = element_text(size=12))

set.seed(00)
ppd<-ptest+geom_text_repel(data=dat,nudge_x=-1,direction="y",force=2,max.iter=4000,hjust=1,segment.size=0.2,size = 2)

ppd
}
#ggsave("T025160_OV_cor_NES.pdf", width = 6, height = 6)


