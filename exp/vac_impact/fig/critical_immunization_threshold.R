library(ggplot2)
library(gridExtra)

%Critical_Immunization_Threshold_0Nat
VE100<-function(x) {(1-(1/x))/1.00}
VE96<-function(x) {(1-(1/x))/0.96}
VE80<-function(x) {(1-(1/x))/0.80}
VE64<-function(x) {(1-(1/x))/0.64}
VE60<-function(x) {(1-(1/x))/0.60}
VE40<-function(x) {(1-(1/x))/0.40}
VE20<-function(x) {(1-(1/x))/0.20}
DAT<-data.frame(x = c(1.0, 3))
XY<-seq(min(DAT$x),max(DAT$x), length.out=11)
D2<-data.frame(XY, VE100(XY), VE80(XY))
D3<-data.frame(XY, VE80(XY), VE60(XY))
D4<-data.frame(XY, VE60(XY), VE40(XY))
D5<-data.frame(XY, VE40(XY), VE20(XY))
D6<-data.frame(XY, VE20(XY), Inf)
Critical_Immunization_Threshold_0Nat<-ggplot(data.frame(x = c(1.0, 3)), aes(x = x)) + stat_function(fun = VE100, aes(color="tomato3"), lwd=1) + stat_function(fun = VE80, aes(color="yellow1"), lwd=1) + stat_function(fun = VE60, aes(color="cyan3"), lwd=1) + stat_function(fun = VE40, aes(color="blue"), lwd=1) + stat_function(fun = VE20, aes(color="purple"), lwd=1) + labs(y="Critical Immunization Percentage", x= "Effective Reproduction Number, R", title= "Critical Immunization Threshold No Natural Immunity") + theme(plot.title=element_text(hjust=0.5)) + scale_y_continuous(expand=c(0,0), name="Critical Percentage of Population", labels = function(x) paste0(x*100, "%")) + scale_x_continuous(expand= c(0,0)) + coord_cartesian(ylim=c(0,1.0))+theme(axis.line=element_line()) + geom_ribbon(aes(x=XY, ymin=VE100(XY), ymax=VE80(XY)), data=D2, fill="#E7861B", alpha=0.5)+ geom_ribbon(aes(x=XY, ymin=VE80(XY), ymax=VE60(XY)), data=D2, fill="#AFA100", alpha=0.5)+ geom_ribbon(aes(x=XY, ymin=VE60(XY), ymax=VE40(XY)), data=D2, fill="#00BF7D",alpha=0.5) + geom_ribbon(aes(x=XY, ymin=VE40(XY), ymax=VE20(XY)), data=D2, fill="#00B8E5", alpha=0.5) + geom_ribbon(aes(x=XY, ymin=VE20(XY), ymax=Inf), data=D2, fill="#AC88FF", alpha=0.5) + stat_function(fun = VE96, aes(color="grey0"), lwd=1) + stat_function(fun = VE64, aes(color="grey1"), lwd=1,linetype=2) + scale_color_identity(name=bquote(VE[SI]),breaks=c("tomato3","grey0", "yellow1","grey1","cyan3","blue","purple"), labels=c("100%","Dose 2 (96%)","80%","Dose 1 (64%)","60%","40%","20%"), guide = guide_legend(override.aes = list(linetype = c(1,1,1,2,1,1,1))))+theme(legend.key.size = unit(0.95, "cm"))

%Plot with partial population immunity of 30%
%Critical_Immunization_Threshold_30Nat
VE100q<-function(x) {(0.70-(1/x))/1.00}
VE96q<-function(x) {(0.70-(1/x))/0.96}
VE80q<-function(x) {(0.70-(1/x))/0.80}
VE64q<-function(x) {(0.70-(1/x))/0.64}
VE60q<-function(x) {(0.70-(1/x))/0.60}
VE40q<-function(x) {(0.70-(1/x))/0.40}
VE20q<-function(x) {(0.70-(1/x))/0.20}
DAT<-data.frame(x = c(1.0, 3))
XY<-seq(min(DAT$x),max(DAT$x), length.out=11)
D2<-data.frame(XY, VE100q(XY), VE80q(XY))
D3<-data.frame(XY, VE80q(XY), VE60q(XY))
D4<-data.frame(XY, VE60q(XY), VE40q(XY))
D5<-data.frame(XY, VE40q(XY), VE20q(XY))
D6<-data.frame(XY, VE20q(XY), Inf)

Critical_Immunization_Threshold_30Nat<-ggplot(data.frame(x = c(1.0, 3)), aes(x = x)) + stat_function(fun = VE100q, aes(color="tomato3"), lwd=1) + stat_function(fun = VE80q, aes(color="yellow1"), lwd=1) + stat_function(fun = VE60q, aes(color="cyan3"), lwd=1) + stat_function(fun = VE40q, aes(color="blue"), lwd=1) + stat_function(fun = VE20q, aes(color="purple"), lwd=1) + labs(y="Critical Immunization Percentage", x= "Effective Reproduction Number, R", title= "Critical Immunization Threshold 30% Natural Immunity") + theme(plot.title=element_text(hjust=0.5)) + scale_y_continuous(expand=c(0,0), name="Critical Percentage of Population", labels = function(x) paste0(x*100, "%")) + scale_x_continuous(expand= c(0,0)) + coord_cartesian(ylim=c(0,1.0))+theme(axis.line=element_line()) + geom_ribbon(aes(x=XY, ymin=VE100q(XY), ymax=VE80q(XY)), data=D2, fill="#E7861B", alpha=0.5)+ geom_ribbon(aes(x=XY, ymin=VE80q(XY), ymax=VE60q(XY)), data=D2, fill="#AFA100", alpha=0.5)+ geom_ribbon(aes(x=XY, ymin=VE60q(XY), ymax=VE40q(XY)), data=D2, fill="#00BF7D",alpha=0.5) + geom_ribbon(aes(x=XY, ymin=VE40q(XY), ymax=VE20q(XY)), data=D2, fill="#00B8E5", alpha=0.5) + geom_ribbon(aes(x=XY, ymin=VE20q(XY), ymax=Inf), data=D2, fill="#AC88FF", alpha=0.5) + stat_function(fun = VE96q, aes(color="grey0"), lwd=1) + stat_function(fun = VE64q, aes(color="grey1"), lwd=1,linetype=2) + scale_color_identity(name=bquote(VE[SI]),breaks=c("tomato3","grey0", "yellow1","grey1","cyan3","blue","purple"), labels=c("100%","Dose 2 (96%)","80%","Dose 1 (64%)","60%","40%","20%"), guide = guide_legend(override.aes = list(linetype = c(1,1,1,2,1,1,1))))+theme(legend.key.size = unit(0.95, "cm"))

%Both
grid.arrange(Critical_Immunization_Threshold_0Nat,Critical_Immunization_Threshold_30Nat,ncol=2)
