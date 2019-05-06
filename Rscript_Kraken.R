args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(scales)

df<-read.csv(file=args[2], header=TRUE, sep=",")

give.n <- function(x){
  return(data.frame(y = sum(x)+150000, label = paste0("n = ",length(x)))) 
  # experiment with the multiplier to find the perfect position
}

for (f in unique(df$Sample)){
  subdf <- df[df$Sample == f,]
  p<-ggplot(subdf,aes(x=factor(Species),y=Length, fill=Coverage)) + geom_bar(colour="black", stat="identity", position="stack", linetype=0) + coord_flip()
  p<-p+scale_fill_gradient(name=expression(atop("Coverage", paste(""["Grey > 60"]))), limits=c(0,60), low="white",high="red")
  p<-p+scale_y_continuous(labels = comma, expand=c(0,250000))
  p<-p+stat_summary(fun.data = give.n, geom = "text", fun.y = min,  position = position_dodge(width = 1.0), color="blue")	
  p<-p+ggtitle(as.character(f))
  p<-p+xlab("Species identified with Kraken")
  p<-p+ylab(expression(atop("Length", paste(""["n = number of contigs"]))))
  filename = paste(as.character(f), ".pdf",sep="")
  ggsave(file=filename, plot=p)
}


p3<-ggplot(df,aes(x=factor(Species),y=Length, fill=Coverage)) + geom_bar(colour="black", stat="identity", position="stack", linetype=0) + coord_flip()+
  scale_fill_gradient(name=expression(atop("Coverage", paste(""["Grey > 60"]))), limits=c(0,60), low="white",high="red")+facet_grid(~Sample)+
  theme(panel.spacing = unit(0.2, "line"))+	
  theme(strip.background=element_blank(),strip.text.x = element_text(angle = 90, size = 8, hjust=0.4), strip.text.y = element_text(size = 6), plot.title=element_text())+
  theme(axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), panel.grid=element_blank(), panel.background=element_blank())+
  labs(title="All samples")+
  xlab("All species identified with Kraken for the entire MiSeq run")
ggsave(file="All_samples.pdf", plot=p3, width=14, height=14)