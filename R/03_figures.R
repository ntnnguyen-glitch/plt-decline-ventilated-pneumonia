# 03 figures: Figure 1 flow chart, Figure 2 Kaplan-Meier, Figure 3 restricted cubic spline, Figure 4 subgroup forest plot (TIFF 300 dpi)
suppressMessages({library(data.table);library(ggplot2);library(splines);library(survival);library(patchwork)})
R<-readRDS("R.rds"); D<-as.data.table(readRDS("Da.rds")); fl<-readRDS("fl.rds"); FT<-"Times New Roman"; dir.create("fig",showWarnings=FALSE)
sv<-function(p,f,w,h){ggsave(file.path("fig",paste0(f,".tiff")),p,device=ragg::agg_tiff,width=w,height=h,units="in",res=300,compression="none",bg="white"); ggsave(file.path("fig",paste0(f,".png")),p,width=w,height=h,dpi=150,bg="white")}
fm<-function(n) format(n,big.mark=","); pv<-function(p) if(p<0.001) "P < 0.001" else sprintf("P = %.3f",p)
## ---- Figure 1: flow chart (arrow heads stop short of the boxes)
bx<-function(x,y,w,h,lab,id) data.table(id=id,x=x,y=y,w=w,h=h,lab=lab)
mk<-function(x0,db,main,exc){B1<-list(); E1<-list(); n<-length(main); ys<-seq(10,10-1.9*(n-1),by=-1.9)
 for(i in 1:n) B1[[i]]<-bx(x0,ys[i],3.3,1.0,main[i],paste0(db,"m",i))
 for(i in seq_along(exc)) if(!is.na(exc[i])) E1[[length(E1)+1]]<-bx(x0+3.6,ys[i]-0.95,2.95,0.85,exc[i],paste0(db,"e",i))
 list(b=rbindlist(B1),e=rbindlist(E1))}
cn<-D[,.N,by=.(db,x)]
M<-mk(0,"M",c(sprintf("MIMIC-IV v3.1\nFirst ICU stays of adults\nn = %s",fm(fl$m0)),sprintf("Pneumonia\nn = %s",fm(fl$m_pna)),sprintf("Invasive ventilation started\n−6 to 24 h after ICU admission\nn = %s",fm(fl$m_imv)),sprintf("Alive at the 96-h landmark\nn = %s",fm(fl$m_lm)),sprintf("Baseline (−6 to 24 h) and 48–96 h\nplatelet counts available\nn = %s",fm(fl$m_pn)),sprintf("Analysed cohort\nn = %s",fm(fl$m_ther))),
 c(NA,NA,sprintf("Death before 96 h:\nn = %s",fm(fl$m_imv-fl$m_lm)),sprintf("Missing platelet counts:\nn = %s",fm(fl$m_lm-fl$m_pn)),sprintf("Haematological malignancy n = %s;\nECMO n = %s; therapeutic\nanticoagulation n = %s",fm(fl$m_pn-fl$m_hem),fm(fl$m_hem-fl$m_ecmo),fm(fl$m_ecmo-fl$m_ther))))
E<-mk(8.4,"E",c(sprintf("eICU-CRD v2.0\nFirst ICU stays of adults\nn = %s",fm(fl$e0)),sprintf("Pneumonia\nn = %s",fm(fl$e_pna)),sprintf("Invasive ventilation started\n−6 to 24 h after ICU admission\nn = %s",fm(fl$e_imv)),sprintf("Alive at the 96-h landmark\nn = %s",fm(fl$e_lm)),sprintf("Baseline (−6 to 24 h) and 48–96 h\nplatelet counts available\nn = %s",fm(fl$e_pn)),sprintf("Analysed cohort\nn = %s",fm(fl$e_ther))),
 c(NA,NA,sprintf("Death before 96 h:\nn = %s",fm(fl$e_imv-fl$e_lm)),sprintf("Missing platelet counts:\nn = %s",fm(fl$e_lm-fl$e_pn)),sprintf("Haematological malignancy n = %s;\nECMO n = %s; therapeutic\nanticoagulation n = %s",fm(fl$e_pn-fl$e_hem),fm(fl$e_hem-fl$e_ecmo),fm(fl$e_ecmo-fl$e_ther))))
yg<-10-1.9*5-1.75
gb<-rbind(bx(-0.95,yg,1.75,0.95,sprintf("Decline <30%%\nn = %s",fm(cn[db=="MIMIC-IV"&x==0]$N)),"Mg0"),bx(0.95,yg,1.75,0.95,sprintf("Decline ≥30%%\nn = %s",fm(cn[db=="MIMIC-IV"&x==1]$N)),"Mg1"),
 bx(7.45,yg,1.75,0.95,sprintf("Decline <30%%\nn = %s",fm(cn[db=="eICU-CRD"&x==0]$N)),"Eg0"),bx(9.35,yg,1.75,0.95,sprintf("Decline ≥30%%\nn = %s",fm(cn[db=="eICU-CRD"&x==1]$N)),"Eg1"))
BX<-rbind(M$b,M$e,E$b,E$e,gb); g<-0.13; ar<-list()
vs<-function(b) for(i in 1:(nrow(b)-1)) ar[[length(ar)+1]]<<-data.table(x=b$x[i],y=b$y[i]-b$h[i]/2,xe=b$x[i+1],ye=b$y[i+1]+b$h[i+1]/2+g)
hz<-function(b,e) for(i in 1:nrow(e)) ar[[length(ar)+1]]<<-data.table(x=b$x[1],y=e$y[i],xe=e$x[i]-e$w[i]/2-g,ye=e$y[i])
sp<-function(b,gx){lb<-b[.N]; for(i in 1:nrow(gx)) ar[[length(ar)+1]]<<-data.table(x=lb$x,y=lb$y-lb$h/2,xe=gx$x[i],ye=gx$y[i]+gx$h[i]/2+g)}
vs(M$b); vs(E$b); hz(M$b,M$e); hz(E$b,E$e); sp(M$b,gb[1:2]); sp(E$b,gb[3:4]); AR<-rbindlist(ar)
p1<-ggplot()+geom_rect(data=BX,aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),fill="white",colour="black",linewidth=.4)+geom_text(data=BX,aes(x,y,label=lab),family=FT,size=2.55,lineheight=.92)+
 geom_segment(data=AR,aes(x=x,y=y,xend=xe,yend=ye),arrow=arrow(length=unit(1.5,"mm"),type="closed"),linewidth=.35)+coord_cartesian(xlim=c(-2,13.6),ylim=c(yg-0.7,10.7),expand=FALSE)+theme_void(base_family=FT)+theme(plot.background=element_rect(fill="white",colour=NA))
sv(p1,"Figure1",8.4,7.4)
## ---- Figure 2: Kaplan-Meier from the 96-h landmark with numbers at risk
D[,grp:=factor(ifelse(x==1,"Decline ≥30%","Decline <30%"),levels=c("Decline <30%","Decline ≥30%"))]
kmp<-function(k,tag){d<-D[db==k]; f<-survfit(Surv(t28/24+4,d28)~grp,d); s<-summary(f,times=seq(4,28,0.25),extend=TRUE)
 z<-data.table(t=s$time,cum=1-s$surv,lo=1-s$upper,hi=1-s$lower,g=sub("grp=","",s$strata)); r<-summary(f,times=c(4,7,14,21,28),extend=TRUE); rk<-data.table(t=r$time,n=r$n.risk,g=sub("grp=","",r$strata))
 pl<-R$t1c[db==k]$p[1]; hr<-R$t1a[db==k&grp==1]
 a<-ggplot(z,aes(t,cum,colour=g))+geom_step(linewidth=.7)+scale_colour_manual(values=c("#4D4D4D","#B2182B"),name=NULL)+scale_x_continuous(breaks=c(4,7,14,21,28),limits=c(3.5,28.5))+scale_y_continuous(labels=function(v) paste0(round(100*v),"%"),limits=c(0,NA))+
  labs(title=paste0(tag," ",k),x=NULL,y="Cumulative mortality")+annotate("text",x=4.3,y=max(z$cum)*0.97,hjust=0,family=FT,size=3.1,label=sprintf("Log-rank %s\nAdjusted HR %s",pv(pl),hr$hr))+
  theme_classic(base_family=FT,base_size=10.5)+theme(legend.position=c(.75,.18),plot.title=element_text(face="bold",size=11))
 b<-ggplot(rk,aes(t,g,label=fm(n)))+geom_text(family=FT,size=3)+scale_x_continuous(breaks=c(4,7,14,21,28),limits=c(3.5,28.5))+labs(x="Days since ICU admission",y=NULL,subtitle="Number at risk")+theme_classic(base_family=FT,base_size=10)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.subtitle=element_text(size=9))
 a/b+plot_layout(heights=c(4,1))}
p2<-(kmp("MIMIC-IV","A")|kmp("eICU-CRD","B")); sv(p2,"Figure2",10,4.8)
## ---- Figure 3: restricted cubic spline, model 3 (MIMIC-IV, eICU-CRD, pooled)
rc<-rbindlist(lapply(c("MIMIC-IV","eICU-CRD","Pooled"),function(k){z<-copy(R$rcs[[k]]); z[,lab:=c("MIMIC-IV"="A MIMIC-IV","eICU-CRD"="B eICU-CRD","Pooled"="C Pooled (stratified by database)")[k]]; z}))
pt<-R$t1c[grepl("RCS",a)]; lab<-sapply(c("MIMIC-IV","eICU-CRD","Pooled"),function(k) sprintf("P overall %s\nP non-linear %s",sub("P ","",pv(pt[db==k]$p[1])),sub("P ","",pv(pt[db==k]$p[2])))); names(lab)<-unique(rc$lab)
an<-data.table(lab=names(lab),txt=lab)
p3<-ggplot(rc,aes(g,or))+geom_ribbon(aes(ymin=lo,ymax=hi),fill="grey85")+geom_line(colour="#B2182B",linewidth=.8)+geom_hline(yintercept=1,linetype=2)+geom_vline(xintercept=30,linetype=3,colour="grey40")+facet_wrap(~lab,nrow=1)+
 geom_text(data=an,aes(x=-Inf,y=Inf,label=txt),hjust=-0.08,vjust=1.3,family=FT,size=3)+scale_y_log10()+labs(x="Relative platelet decline from baseline to 48–96 h (%; negative values indicate a rise)",y="Adjusted OR for 28-day mortality (log scale)")+
 theme_classic(base_family=FT,base_size=10.5)+theme(strip.background=element_blank(),strip.text=element_text(face="bold",hjust=0))
sv(p3,"Figure3",10,3.8)
## ---- Figure 4: subgroup forest plot (pooled)
sg<-copy(R$t4b); sg[,lev2:=ifelse(lev=="Yes","Yes","No")]
nm<-c("Age ≥65 years"="Age ≥65 years","Male sex"="Male sex","Sepsis"="Sepsis","Vasopressor use (0–24 h)"="Vasopressor use, 0–24 h","Baseline platelet <150 ×10⁹/L"="Baseline platelet count <150 × 10⁹/L","Pharmacological thromboprophylaxis"="Pharmacological thromboprophylaxis","COPD"="COPD","ARDS"="ARDS","Severity above database mean"="Severity score above database mean","Database: eICU-CRD"="eICU-CRD (vs MIMIC-IV)")
sg[,sg2:=nm[sg]]; rows<-list(); for(s in unique(sg$sg2)){z<-sg[sg2==s]; rows[[length(rows)+1]]<-data.table(lab=s,est=NA,lo=NA,hi=NA,n="",ci="",pint=sprintf("%.3f",z$pint[1]),b=TRUE); for(l in c("Yes","No")){w<-z[lev2==l]; rows[[length(rows)+1]]<-data.table(lab=paste0("    ",l),est=w$est,lo=w$lo,hi=w$hi,n=paste0(fm(w$ev),"/",fm(w$n)),ci=sprintf("%.2f (%.2f–%.2f)",w$est,w$lo,w$hi),pint="",b=FALSE)}}
FR<-rbindlist(rows); FR[,y:=rev(seq_len(.N))]
th<-theme_void(base_family=FT,base_size=10)
ylim<-c(0.5,nrow(FR)+1.2); hd<-function(t,x) annotate("text",x=x,y=nrow(FR)+1,label=t,family=FT,fontface="bold",size=3.1,hjust=0)
pa<-ggplot(FR)+geom_text(aes(0,y,label=lab,fontface=ifelse(b,"bold","plain")),family=FT,size=3,hjust=0)+geom_text(aes(1.45,y,label=n),family=FT,size=3,hjust=0)+hd("Subgroup",0)+hd("Deaths/N",1.45)+coord_cartesian(xlim=c(0,2.1),ylim=ylim,expand=FALSE)+th
pb<-ggplot(FR[!is.na(est)])+geom_vline(xintercept=1,linetype=2)+geom_errorbar(aes(y=y,xmin=lo,xmax=hi),width=0.25,orientation="y")+geom_point(aes(est,y),shape=15,size=2)+scale_x_log10(breaks=c(0.5,1,2,3))+coord_cartesian(xlim=c(0.5,3.2),ylim=ylim,expand=FALSE)+labs(x="Adjusted OR (log scale)")+theme_classic(base_family=FT,base_size=10)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
pc<-ggplot(FR)+geom_text(aes(0,y,label=ci),family=FT,size=3,hjust=0)+geom_text(aes(1.25,y,label=pint),family=FT,size=3,hjust=0)+hd("OR (95% CI)",0)+hd("P for interaction",1.25)+coord_cartesian(xlim=c(0,2.1),ylim=ylim,expand=FALSE)+th
p4<-pa+pb+pc+plot_layout(widths=c(1.25,1,1.1)); sv(p4,"Figure4",10,6.4)
cat("figures done\n")
