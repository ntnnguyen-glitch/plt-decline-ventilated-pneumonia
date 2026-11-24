# 02 analysis: multiple imputation (20 datasets per database), primary and secondary outcomes, confounding control,
# sensitivity and subgroup analyses; estimates pooled with Rubin's rules
suppressMessages({library(data.table);library(survival);library(splines);library(WeightIt);library(MatchIt);library(mice);library(pROC);library(EValue);library(sandwich);library(lmtest);library(quantreg)})
set.seed(2026); D<-as.data.table(readRDS("D.rds")); D[db=="eICU-CRD",vfd:=NA]; R<-list(); M<-20
D[,sevz:=(sev-mean(sev,na.rm=TRUE))/sd(sev,na.rm=TRUE),by=db]
D[,`:=`(lplt0=log(plt0),lpf=log(pmin(pmax(pf,20),700)),old=as.integer(age>=65),lowplt=as.integer(plt0<150),DB=as.integer(db=="eICU-CRD"),hisev=as.integer(sevz>0))]
iv<-c("d28","hd","x","dec","age","male","chf","ckd","copd","dm","liver","cancer","sepsis","vaso","proph","ards","antipl","vaso96","rrt96","sevz","lplt0","inr","lact","creat","bun","hb","wbc","map","hr","lpf","alb","bili","cr96","bili96","hb96")
## ---- missing data description
cc<-c("inr","lact","lpf","sevz","creat","hb","wbc","map","hr")
D[,complete:=as.integer(complete.cases(D[,cc,with=FALSE]))]
R$miss<-D[,.(n=.N,age=median(age),male=mean(male),sepsis=mean(sepsis),vaso=mean(vaso),x=mean(x),d28=mean(d28)),by=.(db,complete)]
R$misspct<-D[,lapply(.SD,function(z) round(100*mean(is.na(z)),1)),by=db,.SDcols=c("inr","lact","creat","hb","wbc","pf","sev","map","hr","alb","bili","cr96","bili96","hb96")]
## ---- multiple imputation within each database
IMP<-list(); for(k in c("MIMIC-IV","eICU-CRD")){d<-D[db==k]; im<-mice(d[,..iv],m=M,maxit=10,seed=2026,printFlag=FALSE)
 IMP[[k]]<-lapply(1:M,function(i){z<-as.data.table(complete(im,i)); oth<-setdiff(names(d),names(z)); z<-cbind(z,d[,..oth]); z[,`:=`(lcr96=log(cr96),lbili96=log(bili96),hbd=hb96-hb,lbili=log(bili))]; z})}
IMP[["Pooled"]]<-lapply(1:M,function(i) rbind(IMP[["MIMIC-IV"]][[i]],IMP[["eICU-CRD"]][[i]],fill=TRUE))
saveRDS(IMP,"IMP.rds"); DB3<-c("MIMIC-IV","eICU-CRD","Pooled")
## ---- helpers
M2<-"age+male+chf+ckd+copd+dm+liver+cancer+sepsis"; M3<-paste(M2,"+sevz+lplt0+inr+lact+creat+bun+hb+wbc+map+hr+vaso+lpf+proph+ards")
cvp<-function(cv,k) if(k=="Pooled") paste(cv,"+DB") else cv
rub<-function(est,se,expo=TRUE){m<-length(est); q<-mean(est); W<-mean(se^2); B<-if(m>1) var(est) else 0; Tv<-W+(1+1/m)*B; df<-if(B<=1e-12) 1e6 else (m-1)*(1+W/((1+1/m)*B))^2
 h<-qt(.975,df)*sqrt(Tv); p<-2*pt(-abs(q/sqrt(Tv)),df); if(expo) c(est=exp(q),lo=exp(q-h),hi=exp(q+h),p=p) else c(est=q,lo=q-h,hi=q+h,p=p)}
co<-function(fit,tr,V=NULL){b<-coef(fit); if(is.null(V)) V<-vcov(fit); nm<-grep(tr,names(b),value=TRUE)[1]; c(b[[nm]],sqrt(V[nm,nm]))}
poolfit<-function(L,f,tr="^x$",expo=TRUE){r<-t(sapply(L,function(d) f(d))); rub(r[,1],r[,2],expo)}
lgf<-function(form,tr="^x$") function(d) co(glm(as.formula(form),binomial,d),tr)
cxf<-function(cv,k,tr="^x$",y="Surv(t28,d28)") function(d) co(coxph(as.formula(paste(y,"~x+",cv,if(k=="Pooled") "+strata(DB)" else "")),d),tr)
fm<-function(v) sprintf("%.2f (%.2f–%.2f)",v[1],v[2],v[3]); fb<-function(v) sprintf("%.2f (%.2f to %.2f)",v[1],v[2],v[3])
L<-function(k,sub=NULL) {z<-IMP[[k]]; if(!is.null(sub)) z<-lapply(z,function(d) d[eval(parse(text=sub))]); z}
## ---- primary
R$t1a<-rbindlist(lapply(DB3,function(k){z<-L(k); d<-z[[1]]
 cr<-poolfit(z,lgf("d28~x")); m2<-poolfit(z,lgf(paste("d28~x+",cvp(M2,k)))); m3<-poolfit(z,lgf(paste("d28~x+",cvp(M3,k)))); hr<-poolfit(z,cxf(cvp(M3,k),k))
 rbindlist(lapply(0:1,function(g) data.table(db=k,grp=g,n=sum(d$x==g),ev=sum(d$d28[d$x==g]))))[,`:=`(crude=c("1",fm(cr)),m2=c("1",fm(m2)),m3=c("1",fm(m3)),p3=c(NA,m3[4]),hr=c("1",fm(hr)),ph=c(NA,hr[4]),est3=c(NA,m3[1]),lo3=c(NA,m3[2]),hi3=c(NA,m3[3]))]}))
cat("PRIMARY (MICE m=20)\n"); print(R$t1a[grp==1,.(db,crude,m2,m3,p3,hr,ph)])
## ---- dose-response
R$t1b<-rbindlist(lapply(DB3,function(k){z<-lapply(L(k),function(d){d<-copy(d); d[,x10:=dec/10]; d[,xc:=cut(dec,c(-Inf,10,30,50,Inf),right=FALSE,labels=c("a","b","c","d"))]; d[,xt:=as.integer(xc)]; d})
 a<-poolfit(z,lgf(paste("d28~x10+",cvp(M3,k)),"^x10$"),"^x10$"); h<-poolfit(z,function(d) co(coxph(as.formula(paste("Surv(t28,d28)~x10+",cvp(M3,k),if(k=="Pooled") "+strata(DB)" else "")),d),"^x10$"))
 cs<-lapply(c("xcb","xcc","xcd"),function(tr) poolfit(z,lgf(paste("d28~xc+",cvp(M3,k)),paste0("^",tr,"$")))); pt<-poolfit(z,lgf(paste("d28~xt+",cvp(M3,k)),"^xt$"))
 rbind(data.table(db=k,a="Per 10% decline, OR (Model 3)",v=fm(a),p=a[4]),data.table(db=k,a="Per 10% decline, HR (Cox Model 3)",v=fm(h),p=h[4]),
  data.table(db=k,a=c("Decline 10-29.9% vs <10%, OR","Decline 30-49.9% vs <10%, OR","Decline >=50% vs <10%, OR"),v=sapply(cs,fm),p=sapply(cs,function(v) v[4])),data.table(db=k,a="P for trend (4 categories)",v="—",p=pt[4]))}))
## ---- KM, RCS (pooled Wald), added value
asm<-function(l) as.mira(l)
# restricted cubic spline basis (Harrell parameterisation): non-linear terms for knots k
rcsb<-function(x,k){n<-length(k); nr<-(k[n]-k[1])^2; sapply(1:(n-2),function(j) (pmax(x-k[j],0)^3-pmax(x-k[n-1],0)^3*(k[n]-k[j])/(k[n]-k[n-1])+pmax(x-k[n],0)^3*(k[n-1]-k[j])/(k[n]-k[n-1]))/nr)}
pw<-function(a,b) tryCatch(D1(a,b)$result[1,4],error=function(e) D2(a,b)$result[1,4])
R$t1c<-rbindlist(lapply(DB3,function(k){z<-L(k); d<-z[[1]]; st<-if(k=="Pooled") "+strata(DB)" else ""
 lr<-survdiff(as.formula(paste("Surv(t28,d28)~x",st)),d); plr<-pchisq(lr$chisq,1,lower.tail=FALSE)
 dd0<-if(k=="Pooled") D else D[db==k]; lo<-quantile(dd0$dec,.01); hi<-quantile(dd0$dec,.99); kn<-quantile(pmin(pmax(dd0$dec,lo),hi),c(.05,.35,.65,.95))
 zz<-lapply(z,function(d){d<-copy(d); d[,dc:=pmin(pmax(dec,lo),hi)]; nl<-rcsb(d$dc,kn); d[,`:=`(nl1=nl[,1],nl2=nl[,2])]; d})
 f1<-asm(lapply(zz,function(d) glm(as.formula(paste("d28~dc+nl1+nl2+",cvp(M3,k))),binomial,d))); f0<-asm(lapply(zz,function(d) glm(as.formula(paste("d28~",cvp(M3,k))),binomial,d))); fl<-asm(lapply(zz,function(d) glm(as.formula(paste("d28~dc+",cvp(M3,k))),binomial,d)))
 f3<-asm(lapply(zz,function(d) glm(as.formula(paste("d28~x+",cvp(M3,k))),binomial,d)))
 po<-pw(f1,f0); pnl<-pw(f1,fl); pl<-pw(f3,f0)
 a0<-median(sapply(seq_along(zz),function(i) auc(zz[[i]]$d28,fitted(f0$analyses[[i]]),quiet=TRUE))); a1<-median(sapply(seq_along(zz),function(i) auc(zz[[i]]$d28,fitted(f3$analyses[[i]]),quiet=TRUE)))
 g<-seq(quantile(zz[[1]]$dc,.02),quantile(zz[[1]]$dc,.98),length=200); Xg<-cbind(g,rcsb(g,kn)); X0<-c(0,rcsb(0,kn)); LL<-sweep(Xg,2,X0)
 cu<-lapply(f1$analyses,function(f){i<-c("dc","nl1","nl2"); list(lp=as.vector(LL%*%coef(f)[i]),se=sqrt(rowSums((LL%*%vcov(f)[i,i])*LL)))})
 cv<-t(sapply(seq_along(g),function(j) rub(sapply(cu,function(u) u$lp[j]),sapply(cu,function(u) u$se[j]))))
 R$rcs[[k]]<<-data.table(db=k,g=g,or=cv[,1],lo=cv[,2],hi=cv[,3])
 data.table(db=k,a=c("Kaplan–Meier log-rank (from 96-h landmark)","RCS (4 knots), overall P","RCS, P for non-linearity","Adding decline ≥30% to Model 3 (pooled Wald test)"),v=c("—","—","—",sprintf("AUC %.3f → %.3f",a0,a1)),p=c(plr,po,pnl,pl))}))
## ---- secondary outcomes
qr1<-function(form) function(d){f<-rq(as.formula(form),tau=.5,data=d); s<-summary(f,se="boot",R=200)$coefficients; c(s["x",1],s["x",2])}
so<-function(k,lab,y,type){z<-L(k); d<-z[[1]]; if(all(is.na(d[[y]]))) return(NULL); z<-lapply(z,function(w) w[!is.na(get(y))]); d<-z[[1]]; fr<-sprintf("%.1f%%",100*mean(d[[y]]))
 if(type=="lg"){r<-poolfit(z,lgf(paste(y,"~x+",cvp(M3,k)))); data.table(db=k,o=lab,f=fr,m="Logistic Model 3, OR",v=fm(r),p=r[4])}
 else if(type=="cox"){hz<-if(y=="d90") 2160 else 8760; r<-poolfit(z,function(w){w<-copy(w); w[,tt:=fifelse(get(y)==1,pmin(dday*24,hz),hz)-96]; w[tt<=0,tt:=0.5]; co(coxph(as.formula(paste("Surv(tt,",y,")~x+",cvp(M3,k))),w),"^x$")}); data.table(db=k,o=lab,f=fr,m="Cox Model 3, HR",v=fm(r),p=r[4])}
 else {r<-poolfit(z,qr1(paste(y,"~x+",cvp(M3,k))),expo=FALSE); data.table(db=k,o=lab,f=sprintf("%.1f (%.1f–%.1f)",median(d[[y]]),quantile(d[[y]],.25),quantile(d[[y]],.75)),m="Median regression Model 3, difference in medians (days)",v=fb(r),p=r[4])}}
R$t2<-rbindlist(unlist(lapply(DB3,function(k) list(so(k,"In-hospital mortality","hd","lg"),so(k,"90-day mortality","d90","cox"),so(k,"1-year mortality","d365","cox"),so(k,"Ventilator-free days to day 28","vfd","q"),so(k,"ICU length of stay","los","q"),so(k,"Platelet transfusion (96 h–day 28)","pltx","lg"),so(k,"Red-cell transfusion (96 h–day 28)","rbcx","lg"))),recursive=FALSE))
R$t2<-R$t2[!(db=="Pooled"&o%in%c("90-day mortality","1-year mortality","Ventilator-free days to day 28"))]
cat("SECONDARY\n"); print(R$t2[,.(db,o,v,p=signif(p,3))])
## ---- confounding control (estimated within each imputed dataset, then pooled)
psf<-function(k) as.formula(paste("x~",cvp(M3,k)))
smd<-function(d,w,k){X<-model.matrix(psf(k),d)[,-1]; max(abs(apply(X,2,function(z){m1<-weighted.mean(z[d$x==1],w[d$x==1]); m0<-weighted.mean(z[d$x==0],w[d$x==0]); (m1-m0)/sqrt((var(z[d$x==1])+var(z[d$x==0]))/2)})))}
R$t3<-rbindlist(lapply(DB3,function(k){z<-L(k); out<-list(); sm<-matrix(NA,M,3); np<-numeric(M)
 for(i in 1:M){d<-copy(z[[i]]); w1<-weightit(psf(k),d,method="glm",estimand="ATE",stabilize=TRUE); d[,iw:=pmin(w1$weights,quantile(w1$weights,.99))]; w2<-weightit(psf(k),d,method="glm",estimand="ATO"); d[,ow:=w2$weights]
  g1<-glm(d28~x,quasibinomial,d,weights=iw); out$iptw<-rbind(out$iptw,co(g1,"^x$",vcovHC(g1,type="HC0")))
  c1<-coxph(as.formula(paste("Surv(t28,d28)~x",if(k=="Pooled") "+strata(DB)" else "")),d,weights=iw,robust=TRUE); out$iptwh<-rbind(out$iptwh,co(c1,"^x$"))
  g2<-glm(as.formula(paste("d28~x+",cvp(M3,k))),quasibinomial,d,weights=iw); out$dr<-rbind(out$dr,co(g2,"^x$",vcovHC(g2,type="HC0")))
  g3<-glm(d28~x,quasibinomial,d,weights=ow); out$ow<-rbind(out$ow,co(g3,"^x$",vcovHC(g3,type="HC0")))
  mt<-matchit(psf(k),d,method="nearest",distance="glm",link="linear.logit",caliper=0.2,ratio=1); md<-as.data.table(match.data(mt)); np[i]<-sum(md$x==1)
  g4<-glm(d28~x,binomial,md); out$psm<-rbind(out$psm,co(g4,"^x$",vcovCL(g4,cluster=~subclass))); g5<-glm(as.formula(paste("d28~x+",cvp(M3,k))),binomial,md); out$psmdr<-rbind(out$psmdr,co(g5,"^x$",vcovCL(g5,cluster=~subclass)))
  if(k=="eICU-CRD"){g6<-glm(as.formula(paste("d28~x+",M3)),binomial,d); out$cl<-rbind(out$cl,co(g6,"^x$",vcovCL(g6,cluster=~hosp)))}
  sm[i,]<-c(smd(d,rep(1,nrow(d)),k),smd(d,d$iw,k),smd(d,d$ow,k))}
 P<-lapply(out,function(r) rub(r[,1],r[,2])); m2<-poolfit(z,lgf(paste("d28~x+",cvp(M2,k)))); m3<-poolfit(z,lgf(paste("d28~x+",cvp(M3,k)))); ev<-evalues.OR(m3[1],m3[2],m3[3],rare=mean(z[[1]]$d28)<0.15)
 rows<-data.table(db=k,a=c("Model 2, OR","Model 3, OR","IPTW (ATE, stabilized, truncated 99th pct), OR","IPTW, HR (Cox, robust SE)","Doubly robust (IPTW + Model 3), OR","Overlap weights (ATO), OR",sprintf("PSM 1:1, caliper 0.2 SD of logit PS (median %d pairs), OR",round(median(np))),"PSM + Model 3 covariate adjustment in matched sample, OR","E-value (Model 3 OR; CI bound)"),
  v=c(fm(m2),fm(m3),fm(P$iptw),fm(P$iptwh),fm(P$dr),fm(P$ow),fm(P$psm),fm(P$psmdr),sprintf("%.2f (%.2f)",ev[2,1],ev[2,2])),p=c(m2[4],m3[4],P$iptw[4],P$iptwh[4],P$dr[4],P$ow[4],P$psm[4],P$psmdr[4],NA))
 rows<-rbind(rows,data.table(db=k,a=sprintf("Max absolute SMD (mean across imputations): unweighted %.3f; IPTW %.3f; overlap %.3f",mean(sm[,1]),mean(sm[,2]),mean(sm[,3])),v="—",p=NA))
 if(k=="eICU-CRD") rows<-rbind(rows,data.table(db=k,a="Model 3 OR, hospital-clustered robust SE",v=fm(P$cl),p=P$cl[4])); rows}))
cat("CONFOUNDING\n"); print(R$t3[,.(db,a,v,p=signif(p,3))])
## ---- sensitivity analyses
TV<-"+vaso96+rrt96+lcr96+lbili96+hbd"
S<-list(); add<-function(k,lab,sub=NULL,pre=NULL,f="x",cv=NULL){z<-L(k,sub); if(!is.null(pre)) z<-lapply(z,function(d){d<-copy(d); eval(parse(text=pre)); d}); if(is.null(cv)) cv<-cvp(M3,k)
 r<-poolfit(z,lgf(paste("d28~",f,"+",cv),paste0("^",f,"$"))); S[[length(S)+1]]<<-data.table(db=k,a=lab,n=nrow(z[[1]]),v=fm(r),p=r[4])}
for(k in DB3){
 add(k,"Decline threshold ≥20% (vs <20%)",pre="d[,x2:=as.integer(dec>=20)]",f="x2"); add(k,"Decline threshold ≥50% (vs <50%)",pre="d[,x2:=as.integer(dec>=50)]",f="x2")
 add(k,"New-onset thrombocytopenia (<150 ×10⁹/L) among baseline ≥150",sub="plt0>=150",pre="d[,x2:=as.integer(nadir<150)]",f="x2")
 add(k,"Nadir window 48–72 h",sub="!is.na(nadir72)",pre="d[,x2:=as.integer((plt0-nadir72)/plt0>=0.3)]",f="x2"); add(k,"Baseline = lowest platelet count 0–24 h",pre="d[,x2:=as.integer((pltmin24-nadir)/pltmin24>=0.3)]",f="x2")
 add(k,"Model 3 + organ support and laboratory change at 48–96 h",cv=paste(cvp(M3,k),TV)); add(k,"Model 3 + antiplatelet therapy",cv=paste(cvp(M3,k),"+antipl"))
 add(k,"Model 3 + albumin and bilirubin",cv=paste(cvp(M3,k),"+alb+lbili"))
 add(k,"Pneumonia in first 3 discharge diagnoses (MIMIC-IV) or recorded within 24 h (eICU-CRD)",sub="pna_top3==1"); add(k,"Excluding aspiration pneumonia",sub="asp==0")
 add(k,"Excluding platelet transfusion 0–96 h",sub="pltx96==0"); add(k,"Excluding therapeutic anticoagulation started 24–96 h",sub="ther96==0")
 add(k,"Excluding sepsis (non-septic pneumonia)",sub="sepsis==0"); add(k,"Excluding cirrhosis",sub="cirr==0"); add(k,"Excluding DIC or HIT diagnosis",sub="dic==0&hit==0")
 add(k,"Baseline platelet ≥150 ×10⁹/L only",sub="plt0>=150"); add(k,"Excluding death before day 7",sub="is.na(dday)|dday>7")}
S<-rbindlist(S)
# single median imputation and complete-case analysis
num<-c("sevz","inr","lact","creat","bun","hb","wbc","map","hr","lpf")
for(k in DB3){d<-if(k=="Pooled") copy(D) else copy(D[db==k]); cc1<-d[complete.cases(d[,..num])]
 for(v in num) d[,(v):=fifelse(is.na(get(v)),median(get(v),na.rm=TRUE),get(v)),by=db]
 r1<-rub(lgf(paste("d28~x+",cvp(M3,k)))(d)[1],lgf(paste("d28~x+",cvp(M3,k)))(d)[2]); r2<-rub(lgf(paste("d28~x+",cvp(M3,k)))(cc1)[1],lgf(paste("d28~x+",cvp(M3,k)))(cc1)[2])
 S<-rbind(S,data.table(db=k,a=c("Single median imputation","Complete-case analysis"),n=c(nrow(d),nrow(cc1)),v=c(fm(r1),fm(r2)),p=c(r1[4],r2[4])))}
R$t4a<-S; cat("SENSITIVITY\n"); print(S[,.(db,a,n,v,p=signif(p,3))])
## ---- subgroups (pooled cohort)
sg<-list(old="Age ≥65 years",male="Male sex",sepsis="Sepsis",vaso="Vasopressor use (0–24 h)",lowplt="Baseline platelet <150 ×10⁹/L",proph="Pharmacological thromboprophylaxis",copd="COPD",ards="ARDS",hisev="Severity above database mean",DB="Database: eICU-CRD")
cvr<-"age+male+chf+ckd+copd+dm+liver+cancer+sepsis+sevz+lplt0+inr+lact+creat+hb+wbc+map+vaso+lpf+proph+ards+DB"
R$t4b<-rbindlist(lapply(names(sg),function(v){cv<-gsub(paste0("\\+",v,"(?![a-z0-9])|^",v,"\\+"),"",cvr,perl=TRUE); if(v=="old") cv<-sub("^age\\+","",cv); if(v=="lowplt") cv<-sub("\\+lplt0","",cv); if(v=="hisev") cv<-sub("\\+sevz","",cv)
 z<-L("Pooled"); if(v=="hisev") z<-lapply(z,function(d) d[!is.na(hisev)]); it<-paste0("^x:",v,"$"); pi<-poolfit(z,lgf(paste("d28~x*",v,"+",cv),it),it)[4]
 rbindlist(lapply(1:0,function(l){zl<-lapply(z,function(d) d[get(v)==l]); r<-poolfit(zl,lgf(paste("d28~x+",cv))); data.table(sg=sg[[v]],lev=ifelse(l==1,"Yes","No"),n=nrow(zl[[1]]),ev=sum(zl[[1]]$d28),v=fm(r),est=r[1],lo=r[2],hi=r[3],p=r[4],pint=pi)}))}))
## ---- baseline table data
R$t0p<-NULL; saveRDS(R,"R.rds"); saveRDS(D,"Da.rds"); cat("ANALYSIS DONE\n"); print(R$t4b[,.(sg,lev,n,v,p=signif(p,3),pint=signif(pint,3))]); print(R$t1b[,.(db,a,v,p=signif(p,3))]); print(R$t1c[,.(db,a,v,p=signif(p,3))]); print(R$miss)
