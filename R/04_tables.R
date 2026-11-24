# 04 tables: Tables 1-3 and supplementary tables (formatted data in T.rds; Word copies in Tables.docx)
suppressMessages({library(data.table);library(officer);library(flextable)})
R<-readRDS("R.rds"); D<-as.data.table(readRDS("D.rds")); T<-list()
fm<-function(n) format(n,big.mark=","); p3<-function(p) ifelse(is.na(p),"",ifelse(p<0.001,"<0.001",sprintf("%.3f",p)))
cp<-function(v,p) ifelse(is.na(p),v,paste0(v,"; P ",ifelse(p<0.001,"< 0.001",paste0("= ",sprintf("%.3f",p)))))
## Table 1
vb<-list(c("age","Age, years","c"),c("male","Male sex","b"),c("chf","Heart failure","b"),c("copd","COPD","b"),c("ckd","Chronic kidney disease","b"),c("dm","Diabetes mellitus","b"),c("liver","Liver disease","b"),c("cancer","Solid malignancy","b"),c("sepsis","Sepsis","b"),c("ards","ARDS","b"),c("asp","Aspiration pneumonia","b"),
 c("plt0","Baseline platelet count, × 10⁹/L","c"),c("nadir","Nadir platelet count 48–96 h, × 10⁹/L","c"),c("dec","Relative platelet decline, %","c"),c("inr","INR","c"),c("lact","Lactate, mmol/L","c"),c("creat","Creatinine, mg/dL","c"),c("hb","Haemoglobin, g/dL","c"),c("wbc","White blood cell count, × 10⁹/L","c"),c("pf","PaO₂/FiO₂, mmHg","c"),c("sev","Severity score*","c"),
 c("vaso","Vasopressor use, 0–24 h","b"),c("proph","Pharmacological thromboprophylaxis, 0–96 h","b"),c("antipl","Antiplatelet therapy, −24 to 96 h","b"),c("pltx96","Platelet transfusion, 0–96 h","b"),c("d28","28-day mortality","b"),c("hd","In-hospital mortality","b"))
cell<-function(z,t) if(t=="b") sprintf("%s (%.1f)",fm(sum(z,na.rm=TRUE)),100*mean(z,na.rm=TRUE)) else sprintf("%.1f (%.1f–%.1f)",median(z,na.rm=TRUE),quantile(z,.25,na.rm=TRUE),quantile(z,.75,na.rm=TRUE))
pp<-function(z,g,t) if(t=="b") suppressWarnings(chisq.test(table(z,g))$p.value) else wilcox.test(z~g)$p.value
t1<-rbindlist(lapply(vb,function(v){r<-data.table(Characteristic=v[2]); for(k in c("MIMIC-IV","eICU-CRD")){d<-D[db==k]; if(v[1]=="sev"&&k=="eICU-CRD") lab<-"" ; r[[paste(k,"<30%")]]<-cell(d[x==0][[v[1]]],v[3]); r[[paste(k,"≥30%")]]<-cell(d[x==1][[v[1]]],v[3]); r[[paste(k,"P")]]<-p3(pp(d[[v[1]]],d$x,v[3]))}; r}))
t1<-rbind(data.table(Characteristic="n",`MIMIC-IV <30%`=fm(sum(D$db=="MIMIC-IV"&D$x==0)),`MIMIC-IV ≥30%`=fm(sum(D$db=="MIMIC-IV"&D$x==1)),`MIMIC-IV P`="",`eICU-CRD <30%`=fm(sum(D$db=="eICU-CRD"&D$x==0)),`eICU-CRD ≥30%`=fm(sum(D$db=="eICU-CRD"&D$x==1)),`eICU-CRD P`=""),t1)
T$t1<-t1; T$miss<-R$misspct; T$missc<-R$miss
## Table 2 primary + secondary
g<-function(k,a) R$t1a[db==k&grp==1][[a]]; gp<-function(k,a) R$t1a[db==k&grp==1][[a]]
ev<-function(k,gr){z<-R$t1a[db==k&grp==gr]; sprintf("%s/%s (%.1f)",fm(z$ev),fm(z$n),100*z$ev/z$n)}
b1<-function(k,pat){z<-R$t1b[db==k&grepl(pat,a,fixed=TRUE)]; cp(z$v,z$p)}
s2<-function(k,oc){z<-R$t2[db==k&grepl(oc,o,fixed=TRUE)]; if(!nrow(z)) "Not available" else cp(z$v,z$p)}
row<-function(lab,f) data.table(Outcome=lab,`MIMIC-IV`=f("MIMIC-IV"),`eICU-CRD`=f("eICU-CRD"))
T$t2<-rbind(row("Primary outcome: 28-day mortality",function(k) ""),row("   Deaths, decline <30%, n/N (%)",function(k) ev(k,0)),row("   Deaths, decline ≥30%, n/N (%)",function(k) ev(k,1)),
 row("   Model 1 (unadjusted), OR (95% CI)",function(k) g(k,"crude")),row("   Model 2, OR (95% CI)",function(k) g(k,"m2")),row("   Model 3, OR (95% CI)",function(k) cp(g(k,"m3"),g(k,"p3"))),row("   Cox model 3 from the 96-h landmark, HR (95% CI)",function(k) cp(g(k,"hr"),g(k,"ph"))),
 row("   Per 10% decline, model 3 OR (95% CI)",function(k) b1(k,"Per 10% decline, OR")),row("   Decline 10–29.9% vs <10%, OR (95% CI)",function(k) b1(k,"10-29.9%")),row("   Decline 30–49.9% vs <10%, OR (95% CI)",function(k) b1(k,"30-49.9%")),row("   Decline ≥50% vs <10%, OR (95% CI)",function(k) b1(k,">=50%")),row("   P for trend across categories",function(k) p3(R$t1b[db==k&grepl("trend",a)]$p)),
 row("Secondary outcomes (decline ≥30% vs <30%, model 3)",function(k) ""),row("   In-hospital mortality, OR (95% CI)",function(k) s2(k,"In-hospital")),row("   90-day mortality, HR (95% CI)",function(k) s2(k,"90-day")),row("   1-year mortality, HR (95% CI)",function(k) s2(k,"1-year")),
 row("   Ventilator-free days to day 28, difference in medians (95% CI), days",function(k) s2(k,"Ventilator-free")),row("   ICU length of stay, difference in medians (95% CI), days",function(k) s2(k,"ICU length")),row("   Platelet transfusion, 96 h to day 28, OR (95% CI)",function(k) s2(k,"Platelet transfusion")),row("   Red-cell transfusion, 96 h to day 28, OR (95% CI)",function(k) s2(k,"Red-cell")))
## Table 3 confounding control + sensitivity
c3<-function(k,pat){z<-R$t3[db==k&grepl(pat,a,fixed=TRUE)]; if(!nrow(z)) "Not applicable" else cp(z$v,z$p)}
s3<-function(k,pat){z<-R$t4a[db==k&a==pat]; sprintf("%s; n = %s",cp(z$v,z$p),fm(z$n))}
ps<-function(k) sub(".*median (\\d+) pairs.*","\\1",R$t3[db==k&grepl("PSM 1:1",a)]$a)
smd<-function(k){z<-R$t3[db==k&grepl("SMD",a)]$a; m<-regmatches(z,gregexpr("[0-9.]+",z))[[1]]; sprintf("%s → %s (IPTW)",m[1],m[2])}
evl<-function(k) R$t3[db==k&grepl("E-value",a)]$v
T$t3<-rbind(row("Confounding control",function(k) ""),row("   Model 3 (main analysis, multiple imputation)",function(k) c3(k,"Model 3, OR")),row("   IPTW, OR",function(k) c3(k,"IPTW (ATE")),row("   IPTW, Cox HR",function(k) c3(k,"IPTW, HR")),row("   Doubly robust (IPTW + model 3), OR",function(k) c3(k,"Doubly robust")),row("   Overlap weighting, OR",function(k) c3(k,"Overlap")),
 row("   Propensity score matching 1:1, OR",function(k) paste0(c3(k,"PSM 1:1"),"; median ",ps(k)," pairs")),row("   Matching plus model 3 adjustment, OR",function(k) c3(k,"PSM + Model 3")),row("   Hospital-clustered standard errors, OR",function(k) c3(k,"hospital-clustered")),row("   Maximum absolute SMD, unweighted → weighted",smd),row("   E-value, point estimate (lower CI limit)",evl),
 row("Exposure definition",function(k) ""),row("   Decline ≥20% vs <20%",function(k) s3(k,"Decline threshold ≥20% (vs <20%)")),row("   Decline ≥50% vs <50%",function(k) s3(k,"Decline threshold ≥50% (vs <50%)")),row("   Nadir window 48–72 h",function(k) s3(k,"Nadir window 48–72 h")),row("   Baseline = lowest count within 0–24 h",function(k) s3(k,"Baseline = lowest platelet count 0–24 h")),row("   New-onset thrombocytopenia (<150 × 10⁹/L), baseline ≥150",function(k) s3(k,"New-onset thrombocytopenia (<150 ×10⁹/L) among baseline ≥150")),
 row("Additional adjustment",function(k) ""),row("   Organ support and laboratory change at 48–96 h",function(k) s3(k,"Model 3 + organ support and laboratory change at 48–96 h")),row("   Antiplatelet therapy",function(k) s3(k,"Model 3 + antiplatelet therapy")),row("   Albumin and bilirubin",function(k) s3(k,"Model 3 + albumin and bilirubin")),
 row("Restricted populations",function(k) ""),row("   Pneumonia as an early or leading diagnosis†",function(k) s3(k,"Pneumonia in first 3 discharge diagnoses (MIMIC-IV) or recorded within 24 h (eICU-CRD)")),row("   No aspiration pneumonia",function(k) s3(k,"Excluding aspiration pneumonia")),row("   No platelet transfusion 0–96 h",function(k) s3(k,"Excluding platelet transfusion 0–96 h")),row("   No therapeutic anticoagulation started 24–96 h",function(k) s3(k,"Excluding therapeutic anticoagulation started 24–96 h")),row("   No DIC or HIT diagnosis",function(k) s3(k,"Excluding DIC or HIT diagnosis")),row("   No cirrhosis",function(k) s3(k,"Excluding cirrhosis")),row("   Baseline platelet count ≥150 × 10⁹/L",function(k) s3(k,"Baseline platelet ≥150 ×10⁹/L only")),row("   No sepsis",function(k) s3(k,"Excluding sepsis (non-septic pneumonia)")),row("   Alive on day 7",function(k) s3(k,"Excluding death before day 7")),
 row("Missing data",function(k) ""),row("   Single median imputation",function(k) s3(k,"Single median imputation")),row("   Complete-case analysis",function(k) s3(k,"Complete-case analysis")))
## Supplementary: pooled analyses and subgroup data
T$sp<-data.table(Analysis=c("Model 3, OR","Cox model 3, HR","IPTW, OR","Overlap weighting, OR","Propensity score matching plus adjustment, OR","Per 10% decline, OR","Decline ≥50% vs <10%, OR","P for trend","RCS P overall / P non-linear","AUC of model 3 without → with decline ≥30% (LR test)"),
 Pooled=c(cp(g("Pooled","m3"),g("Pooled","p3")),cp(g("Pooled","hr"),g("Pooled","ph")),c3("Pooled","IPTW (ATE"),c3("Pooled","Overlap"),c3("Pooled","PSM + Model 3"),b1("Pooled","Per 10% decline, OR"),b1("Pooled",">=50%"),p3(R$t1b[db=="Pooled"&grepl("trend",a)]$p),
 paste(p3(R$t1c[db=="Pooled"]$p[2]),"/",p3(R$t1c[db=="Pooled"]$p[3])),paste0(R$t1c[db=="Pooled"]$v[4],"; P ",ifelse(R$t1c[db=="Pooled"]$p[4]<0.001,"< 0.001",sprintf("= %.3f",R$t1c[db=="Pooled"]$p[4])))))
T$inc<-R$t1c[,.(db,a,v,p)]
T$sg<-R$t4b[,.(Subgroup=sg,Level=lev,`Deaths/N`=paste0(fm(ev),"/",fm(n)),`OR (95% CI)`=v,P=p3(p),`P for interaction`=p3(pint))]
saveRDS(T,"T.rds")
set_flextable_defaults(font.family="Times New Roman",font.size=9,padding=2); ft<-function(d) fit_to_width(autofit(bold(theme_booktabs(flextable(as.data.frame(d))),part="header")),6.7)
doc<-read_docx(); for(n in c("t1","t2","t3")){doc<-body_add_par(doc,n); doc<-body_add_flextable(doc,ft(T[[n]])); doc<-body_add_break(doc)}; print(doc,target="Tables.docx"); cat("tables done\n")
