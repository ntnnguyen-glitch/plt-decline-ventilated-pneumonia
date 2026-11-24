# 01 cohort: adults with pneumonia, invasive ventilation started within 24 h of ICU admission, alive at the 96-h landmark
# Exposure: early platelet decline = (baseline platelet [first value -6 to 24 h] - nadir 48-96 h) / baseline >= 30%
suppressMessages({library(DBI);library(data.table)})
con<-dbConnect(duckdb::duckdb(shared_home=FALSE),":memory:"); invisible(dbExecute(con,sprintf("ATTACH '%s' AS w (READ_ONLY)",Sys.getenv("DUCKDB","w.duckdb"))))
Q<-function(s) as.data.table(dbGetQuery(con,s)); H<-function(a,b) as.numeric(difftime(a,b,units="hours")); fl<-list(); LM<-96
# ================= MIMIC-IV =================
m<-Q("select * from w.st"); fl$m0<-nrow(m)
dx<-Q("select hadm_id,icd_code from w.dx"); fx<-function(p) as.integer(m$hadm_id%in%dx[grepl(p,icd_code)]$hadm_id)
m[,pna:=fx("^(48[0-6]|J1[2-8]|J69|5070)")]; m<-m[pna==1]; fl$m_pna<-nrow(m)
pe<-Q("select stay_id,starttime,endtime from w.pe where itemid=225792"); pe<-merge(pe[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); pe[,`:=`(s=H(starttime,intime),e=H(endtime,intime))]
m<-m[stay_id%in%pe[s>=-6&s<=24]$stay_id]; fl$m_imv<-nrow(m)
dth<-as.POSIXct(m$deathtime,tz="UTC"); dd<-as.POSIXct(ifelse(is.na(m$dod),NA,paste(as.character(m$dod),"12:00:00")),tz="UTC"); dth[is.na(dth)]<-dd[is.na(dth)]; m[,death_h:=pmax(H(dth,intime),0)]
# out-of-hospital deaths have a date only: day counts use calendar dates (date of death - date of ICU admission)
m[,dday:=fifelse(!is.na(deathtime),as.numeric(H(as.POSIXct(deathtime,tz="UTC"),intime))/24,fifelse(!is.na(dod),as.numeric(as.Date(dod)-as.Date(intime)),NA_real_))]
m<-m[is.na(death_h)|death_h>LM]; fl$m_lm<-nrow(m)
ml<-as.data.table(readRDS("h_mlab.rds")); ml<-merge(ml[hadm_id%in%m$hadm_id&!is.na(v)],m[,.(hadm_id,stay_id,intime)],by="hadm_id"); ml[,h:=H(charttime,intime)]
pb<-ml[itemid==51265&v>=1&v<=1500&h>=-6&h<=24][order(h)][,.(plt0=v[1],pltmin24=min(v)),by=stay_id]; pn<-ml[itemid==51265&v>=1&v<=1500&h>48&h<=LM,.(nadir=min(v)),by=stay_id]
pn72<-ml[itemid==51265&v>=1&v<=1500&h>48&h<=72,.(nadir72=min(v)),by=stay_id]; pd7<-ml[itemid==51265&v>=1&v<=1500&h>LM&h<=168,.(plt7=min(v)),by=stay_id]
m<-merge(m,pb,by="stay_id"); fl$m_pb<-nrow(m); m<-merge(m,pn,by="stay_id"); fl$m_pn<-nrow(m)
m[,`:=`(nadir72=pn72$nadir72[match(stay_id,pn72$stay_id)],plt7=pd7$plt7[match(stay_id,pd7$stay_id)])]
fi<-function(it,lo,hi,f){x<-ml[itemid==it&v>=lo&v<=hi&h>=-6&h<=24,.(v=f(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(inr=fi(51237,.5,20,max),ptt=fi(51275,10,200,max))]
m[,`:=`(hem=fx("^(20[0-8]|C8[1-9]|C9[0-6])"),cirr=fx("^(5712|5715|5716|K703|K717|K74[3-6])"))]
pr<-as.data.table(readRDS("h_mproc.rds")); m[,ecmo:=as.integer(hadm_id%in%pr[grepl("^(3965|5A15)",icd_code)]$hadm_id)]
m<-m[hem==0]; fl$m_hem<-nrow(m); m<-m[ecmo==0]; fl$m_ecmo<-nrow(m)
rx<-as.data.table(readRDS("h_mrx.rds")); rx<-merge(rx[hadm_id%in%m$hadm_id],m[,.(hadm_id,stay_id,intime)],by="hadm_id"); rx[,h:=H(starttime,intime)]; rx[,d:=tolower(drug)]; rx[,r:=toupper(route)]; rx<-rx[!grepl("flush|lock|dwell",d)]; rx[,dv:=suppressWarnings(as.numeric(sub("-.*","",dose)))]
rx[,lmwh_ther:=(grepl("enoxaparin",d)&!is.na(dv)&dv>40)|(grepl("dalteparin",d)&!is.na(dv)&dv>5000)|(grepl("fondaparinux",d)&!is.na(dv)&dv>2.5)]
ther<-rx[(grepl("heparin",d)&!grepl("enoxaparin|dalteparin",d)&r%in%c("IV","IV DRIP"))|grepl("argatroban|bivalirudin|apixaban|rivaroxaban|warfarin|dabigatran|edoxaban",d)|lmwh_ther]
proph<-rx[!lmwh_ther&((grepl("heparin",d)&!grepl("enoxaparin|dalteparin",d)&r%in%c("SC","SUBCUT"))|grepl("enoxaparin|dalteparin|fondaparinux",d))]; apl<-rx[grepl("aspirin|clopidogrel|ticagrelor|prasugrel|cangrelor",d)]
m[,ther0:=as.integer(stay_id%in%ther[h>=-24&h<=24]$stay_id)]; m<-m[ther0==0]; fl$m_ther<-nrow(m)
m[,`:=`(ther96=as.integer(stay_id%in%ther[h>24&h<=LM]$stay_id),proph=as.integer(stay_id%in%proph[h>=-12&h<=LM]$stay_id),ufh=as.integer(stay_id%in%proph[grepl("heparin",d)&!grepl("enoxaparin|dalteparin",d)&h>=-12&h<=LM]$stay_id),antipl=as.integer(stay_id%in%apl[h>=-24&h<=LM]$stay_id))]
m[,`:=`(chf=fx("^(428|I50)"),copd=fx("^(49[0-2]|496|J4[0-4])"),ckd=fx("^(585|N18)"),dm=fx("^(250|E1[0-4])"),liver=fx("^(571|K70|K7[2-6])"),cancer=fx("^(1[4-9][0-9]|19[0-9]|C[0-7][0-9]|C80)"),sepsis=fx("^(99591|99592|78552|A40|A41|R652)"),ards=fx("^(51882|J80)"),af=fx("^(42731|I48)"),dic=fx("^(2866|D65)"),hit=fx("^(28981|D7582)"),bleed=fx("^(578|K92[0-2]|430|431|432|I6[0-2]|5997|R31|7863|R042)"),asp=fx("^(5070|J69)"),male=as.integer(gender=="M"))]
L<-Q("select stay_id,specimen_id,itemid,charttime,v,value from w.lab"); L<-merge(L[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); L[,h:=H(charttime,intime)]
gv<-function(it,lo,hi,fun="first",a=-6,z=24){x<-L[itemid%in%it&v>=lo&v<=hi&h>=a&h<=z][order(h)]; x<-if(fun=="first") x[,.(v=v[1]),by=stay_id] else x[,.(v=if(fun=="max") max(v) else min(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(lact=gv(50813,.1,30),bun=gv(51006,1,250),creat=gv(50912,.1,25),hb=gv(51222,3,20),wbc=gv(51301,.1,200),alb=gv(50862,1,6),bili=gv(50885,.1,60,"max"),crmax=gv(50912,.1,25,"max",0,24))]
sp<-L[itemid==52033,.(specimen_id,st=value)]; g<-dcast(L[itemid%in%c(50818,50820,50821,50816)&h>=-6&h<=24],stay_id+specimen_id+h~itemid,value.var="v",fun.aggregate=function(x)x[1]); setnames(g,c("50818","50820","50821","50816"),c("paco2","ph","pao2","fl"),skip_absent=TRUE)
g[,st:=sp$st[match(specimen_id,sp$specimen_id)]]; g<-g[(is.na(st)|st=="ART.")&!is.na(pao2)&pao2>=20&pao2<=700]
ce<-Q("select stay_id,itemid,charttime,v from w.ce"); ce<-merge(ce[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); ce[,h:=H(charttime,intime)]
f2<-ce[itemid==223835&v>=21&v<=100]; pc2<-merge(g[,.(stay_id,specimen_id,ha=h)],f2[,.(stay_id,hf=h,f=v)],by="stay_id",allow.cartesian=TRUE)[abs(ha-hf)<=2][order(abs(ha-hf))][,.SD[1],by=specimen_id]
g[,fc:=pc2$f[match(specimen_id,pc2$specimen_id)]]; g[,f:=fifelse(!is.na(fl)&fl>=21&fl<=100,fl,fc)]; pp<-g[!is.na(f)][order(h)][,.SD[1],by=stay_id]
m<-merge(m,pp[,.(stay_id,pao2,paco2,ph,fio2=f)],by="stay_id",all.x=TRUE); m[,pf:=pao2/(fio2/100)]
cm<-function(its,lo,hi,fun=mean,a=0,z=24){x<-ce[itemid%in%its&v>=lo&v<=hi&h>=a&h<=z,.(v=fun(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(hr=cm(220045,20,250),map=cm(c(220052,220181,225312),20,200),mapmin=cm(c(220052,220181,225312),20,200,min),rr=cm(c(220210,224690),3,70))]
ip<-as.data.table(readRDS("h_minp.rds")); ip<-ip[!(itemid%in%c(225170,225168,220970))|(!is.na(amount)&amount>0&!statusdescription%in%"Rewritten")]; ip<-merge(ip[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); ip[,h:=H(starttime,intime)]
m[,`:=`(vaso=as.integer(stay_id%in%ip[itemid%in%c(221906,222315,221749,221662,221289)&h>=-6&h<=24]$stay_id),pltx96=as.integer(stay_id%in%ip[itemid==225170&h>=0&h<=LM]$stay_id),pltx=as.integer(stay_id%in%ip[itemid==225170&h>LM&h<=672]$stay_id),rbcx=as.integer(stay_id%in%ip[itemid==225168&h>LM&h<=672]$stay_id))]
fl2<-function(it,lo,hi,f){x<-ml[itemid==it&v>=lo&v<=hi&h>48&h<=LM,.(v=f(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(cr96=fl2(50912,.1,25,max),bili96=fl2(50885,.1,60,max),hb96=fl2(51222,3,20,min))]
rr<-as.data.table(readRDS("h_mrrt.rds")); rr[,stay_id:=as.numeric(stay_id)]; rr<-merge(rr,m[,.(stay_id=as.numeric(stay_id),intime)],by="stay_id"); rr[,h:=H(starttime,intime)]
v96<-ip[itemid%in%c(221906,222315,221749,221662,221289)&h>24&h<=LM]$stay_id; r96<-rr[h>=0&h<=LM]$stay_id
m[,vaso96:=as.integer(stay_id%in%v96)]; m[,rrt96:=as.integer(as.numeric(stay_id)%in%r96)]
md<-as.data.table(readRDS("h_mdx.rds")); m[,pna_top3:=as.integer(hadm_id%in%md[seq<=3&grepl("^(48[0-6]|J1[2-8]|J69|5070)",icd_code)]$hadm_id)]
sc<-function(x,br) fifelse(is.na(x),0L,as.integer(findInterval(x,br)))
m[,sev:=as.numeric(sc(bili,c(1.2,2,6,12))+fifelse(vaso==1,3L,fifelse(!is.na(mapmin)&mapmin<70,1L,0L))+sc(crmax,c(1.2,2,3.5,5)))]
v28<-pe[stay_id%in%m$stay_id&s<672,.(vh=sum(pmin(e,672)-pmax(s,0))),by=stay_id]; m[,vh:=v28$vh[match(stay_id,v28$stay_id)]]
m[,`:=`(vfd=fifelse(!is.na(death_h)&death_h<=672,0,pmax(0,28-fifelse(is.na(vh),0,vh)/24)),los=H(outtime,intime)/24,hd=hospital_expire_flag,db="MIMIC-IV",id=stay_id,hosp=NA_character_)]
# ================= eICU-CRD (diagnosis strings are stored in lower case) =================
e<-Q("select cast(patientunitstayid as bigint) s,uniquepid,gender,age,try_cast(hospitaldischargeoffset as int) hdo,hospitaldischargestatus hds,try_cast(unitdischargeoffset as int) udo,try_cast(unitvisitnumber as int) uvn,hospitaldischargeyear yr,hospitalid hosp from w.ept")
e<-e[uvn==1][order(yr,s)][,.SD[1],by=uniquepid]; e[,age:=as.numeric(ifelse(age=="> 89","90",age))]; e<-e[!is.na(age)&age>=18&hds%in%c("Alive","Expired")&!is.na(hdo)]; fl$e0<-nrow(e)
edx<-Q("select s,d from w.edx"); eph<-Q("select s,d from w.eph"); fe<-function(p) as.integer(e$s%in%c(edx[grepl(p,d,ignore.case=TRUE)]$s,eph[grepl(p,d,ignore.case=TRUE)]$s)); fd<-function(p) as.integer(e$s%in%edx[grepl(p,d,ignore.case=TRUE)]$s); fdc<-function(p) as.integer(e$s%in%edx[grepl(p,d,perl=TRUE)]$s)
e<-e[s%in%edx[grepl("pneumonia",d,ignore.case=TRUE)]$s]; fl$e_pna<-nrow(e)
rc<-Q("select s,vs,ve from w.erc"); e<-e[s%in%rc[vs>=-360&vs<=1440]$s]; fl$e_imv<-nrow(e)
e[,death_h:=fifelse(hds=="Expired",pmax(hdo/60,0),NA_real_)]; e<-e[is.na(death_h)|death_h>LM]; fl$e_lm<-nrow(e)
el<-as.data.table(readRDS("h_elab.rds")); el<-el[s%in%e$s&!is.na(v)]; ep<-el[n=="platelets x 1000"&v>=1&v<=1500]
pb<-ep[o>=-360&o<=1440][order(o)][,.(plt0=v[1],pltmin24=min(v)),by=s]; pn<-ep[o>2880&o<=LM*60,.(nadir=min(v)),by=s]; pn72<-ep[o>2880&o<=4320,.(nadir72=min(v)),by=s]; pd7<-ep[o>LM*60&o<=10080,.(plt7=min(v)),by=s]
e<-merge(e,pb,by="s"); fl$e_pb<-nrow(e); e<-merge(e,pn,by="s"); fl$e_pn<-nrow(e); e[,`:=`(nadir72=pn72$nadir72[match(s,pn72$s)],plt7=pd7$plt7[match(s,pd7$s)])]
ge2<-function(nm,lo,hi){x<-el[n==nm&v>=lo&v<=hi&o>=-360&o<=1440,.(v=max(v)),by=s]; x$v[match(e$s,x$s)]}; e[,`:=`(inr=ge2("PT - INR",.5,20),ptt=ge2("PTT",10,200))]
e[,`:=`(hem=fe("leukemia|lymphoma|myeloma|myelodysplas"),cirr=fe("cirrhosis"))]; tr<-as.data.table(readRDS("h_etr.rds")); e[,ecmo:=as.integer(s%in%tr[grepl("ecmo",d,ignore.case=TRUE)]$s)]
e<-e[hem==0]; fl$e_hem<-nrow(e); e<-e[ecmo==0]; fl$e_ecmo<-nrow(e)
em<-as.data.table(readRDS("h_emed.rds")); em[,d:=tolower(drugname)]; em[,r:=tolower(routeadmin)]; em<-em[!grepl("flush|lock",d)&s%in%e$s]
num1<-function(z){z[is.na(z)]<-""; x<-regmatches(z,regexpr("[0-9]+(\\.[0-9]+)?(?=\\s*mg)",z,perl=TRUE)); out<-rep(NA_real_,length(z)); out[regexpr("[0-9]+(\\.[0-9]+)?(?=\\s*mg)",z,perl=TRUE)>0]<-as.numeric(x); out}
em[,dv:=num1(tolower(dosage))]; em[is.na(dv),dv:=num1(d)]; em[,lmwh_ther:=grepl("enoxaparin|lovenox",d)&!is.na(dv)&dv>40]
inf<-Q("select s,o,d from w.einf"); inf<-inf[s%in%e$s]
thr<-rbind(em[(grepl("heparin",d)&(grepl("^iv",r)|grepl("d5w|nacl|1/2 ns|25,000",d)))|grepl("argatroban|bivalirudin|warfarin|coumadin|apixaban|eliquis|rivaroxaban|xarelto|dabigatran",d)|lmwh_ther,.(s,o)],inf[grepl("heparin|argatroban|bivalirudin",d,ignore.case=TRUE),.(s,o)])
e[,ther0:=as.integer(s%in%thr[o>=-1440&o<=1440]$s)]; e<-e[ther0==0]; fl$e_ther<-nrow(e)
pr<-em[!lmwh_ther&((grepl("heparin",d)&grepl("sub|sc|sq",r))|grepl("enoxaparin|lovenox|dalteparin|fragmin|fondaparinux",d))]; apl<-em[grepl("aspirin|clopidogrel|plavix|ticagrelor|brilinta|prasugrel|effient",d)]
e[,`:=`(ther96=as.integer(s%in%thr[o>1440&o<=LM*60]$s),proph=as.integer(s%in%pr[o>=-720&o<=LM*60]$s),ufh=as.integer(s%in%pr[grepl("heparin",d)&o>=-720&o<=LM*60]$s),antipl=as.integer(s%in%apl[o>=-1440&o<=LM*60]$s))]
e[,`:=`(chf=fe("heart failure"),copd=fe("copd|chronic obstructive|emphysema"),ckd=fe("chronic kidney|chronic renal|renal insufficiency"),dm=fe("diabetes"),liver=fe("cirrhosis|chronic liver"),cancer=fe("cancer|malignan|carcinoma|metasta"),sepsis=fd("sepsis|septic"),ards=fdc("\\bards\\b|acute respiratory distress"),af=fe("atrial fibrillation"),dic=fdc("disseminated intravascular|\\bdic\\b"),hit=fdc("heparin.induced|\\bhit\\b"),bleed=fd("hemorrhage|bleeding"),asp=fd("aspiration"),male=as.integer(gender=="Male"))]
EL<-Q("select s,o,n,v from w.elab where o between -360 and 1440"); EL<-EL[s%in%e$s&!is.na(v)][order(o)]
ge<-function(nm,lo,hi,f="first"){x<-EL[n==nm&v>=lo&v<=hi]; x<-if(f=="first") x[,.(v=v[1]),by=s] else x[,.(v=max(v)),by=s]; x$v[match(e$s,x$s)]}
e[,`:=`(lact=ge("lactate",.1,30),bun=ge("BUN",1,250),creat=ge("creatinine",.1,25),hb=ge("Hgb",3,20),wbc=ge("WBC x 1000",.1,200),alb=ge("albumin",1,6),bili=ge("total bilirubin",.1,60,"max"),paco2=ge("paCO2",10,150),ph=ge("pH",6.5,8))]
po<-EL[n=="paO2"&v>=20&v<=700][,.SD[1],by=s]; r2<-Q("select s,o,l,v from w.erc2 where l in ('FiO2','FIO2 (%)') and o between -720 and 1440"); r2[,v:=suppressWarnings(as.numeric(v))]
fi2<-rbind(r2[!is.na(v)&v>=21&v<=100,.(s,of=o,f=v)],EL[n=="FiO2"&v>=21&v<=100,.(s,of=o,f=v)]); pp<-merge(po[,.(s,o,pao2=v)],fi2,by="s",allow.cartesian=TRUE)[of<=o+15&o-of<=240][order(s,o-of)][,.SD[1],by=s]
e[,`:=`(pao2=pp$pao2[match(s,pp$s)],fio2=pp$f[match(s,pp$s)])]; e[,pf:=pao2/(fio2/100)]
aps<-Q("select cast(patientunitstayid as bigint) s,try_cast(heartrate as double) hr,try_cast(meanbp as double) map,try_cast(respiratoryrate as double) rr from w.eaps"); a<-aps[match(e$s,aps$s)]
e[,`:=`(hr=fifelse(a$hr>0,a$hr,NA_real_),map=fifelse(a$map>0,a$map,NA_real_),rr=fifelse(a$rr>0,a$rr,NA_real_))]
ap<-Q("select * from w.eapr"); e[,sev:=as.numeric(ap$apache[match(s,ap$s)])]; e[sev<0,sev:=NA]
e[,vaso:=as.integer(s%in%inf[grepl("norepi|levophed|vasopressin|phenyleph|neo-syn|dopamine|epineph",d,ignore.case=TRUE)&o>=-360&o<=1440]$s)]
io<-as.data.table(readRDS("h_eio.rds")); io<-io[!is.na(v)&v>0]; tx<-rbind(tr[grepl("transfusion|blood product",d,ignore.case=TRUE),.(s,o,d)],io[,.(s,o,d)]); tx[,d:=tolower(d)]
e[,`:=`(pltx96=as.integer(s%in%tx[grepl("platelet",d)&o>=0&o<=LM*60]$s),pltx=as.integer(s%in%tx[grepl("platelet",d)&o>LM*60&o<=40320]$s),rbcx=as.integer(s%in%tx[grepl("red blood|prbc|packed red",d)&o>LM*60&o<=40320]$s))]
e3<-function(nm,lo,hi,f){x<-el[n==nm&v>=lo&v<=hi&o>2880&o<=LM*60,.(v=f(v)),by=s]; x$v[match(e$s,x$s)]}
e[,`:=`(cr96=e3("creatinine",.1,25,max),bili96=e3("total bilirubin",.1,60,max),hb96=e3("Hgb",3,20,min))]
einf<-as.data.table(readRDS("h_einf.rds")); ev96<-einf[o>1440&o<=LM*60]$s; er96<-tr[grepl("dialysis|crrt",d,ignore.case=TRUE)&o>=0&o<=LM*60]$s
e[,vaso96:=as.integer(s%in%ev96)]; e[,rrt96:=as.integer(s%in%er96)]
ep2<-as.data.table(readRDS("h_epna.rds")); e[,pna_top3:=as.integer(s%in%ep2[o<=1440]$s)]
rv<-merge(rc[s%in%e$s],e[,.(s)],by="s")[vs<40320,.(vh=sum(pmax(0,pmin(ve,40320)-pmax(vs,0)))/60),by=s]; e[,vh:=rv$vh[match(s,rv$s)]]
e[,`:=`(vfd=fifelse(!is.na(death_h)&death_h<=672,0,pmax(0,28-fifelse(is.na(vh),0,vh)/24)),los=udo/1440,hd=as.integer(hds=="Expired"),db="eICU-CRD",id=s)]
e[,disch_h:=hdo/60]; m[,disch_h:=H(dischtime,intime)]; e[,dday:=fifelse(hds=="Expired",hdo/1440,NA_real_)]
cols<-c("id","db","hosp","age","male","chf","copd","ckd","dm","liver","cirr","cancer","sepsis","ards","af","asp","dic","hit","bleed","plt0","pltmin24","nadir","nadir72","plt7","inr","ptt","lact","bun","creat","hb","wbc","alb","bili","pao2","paco2","ph","fio2","pf","hr","map","rr","sev","vaso","proph","ufh","ther96","pltx96","pltx","rbcx","vfd","los","hd","death_h","disch_h","dday","antipl","cr96","bili96","hb96","vaso96","rrt96","pna_top3")
D<-rbind(m[,..cols],e[,..cols]); D[,dec:=(plt0-nadir)/plt0*100]; D[,x:=as.integer(dec>=30)]
D[,`:=`(d28=as.integer(!is.na(dday)&dday<=28),d90=as.integer(!is.na(dday)&dday<=90),d365=as.integer(!is.na(dday)&dday<=365))]
D[,death_h:=fifelse(!is.na(death_h),pmin(death_h,fifelse(d28==1,672,death_h)),death_h)]
D[,t28:=fifelse(d28==1,death_h,672)-LM]  # eICU-CRD: patients discharged alive are treated as 28-day survivors, consistent with the binary outcome
D[db=="eICU-CRD",`:=`(d90=NA,d365=NA)]
saveRDS(D,"D.rds"); saveRDS(fl,"fl.rds"); print(unlist(fl))
print(D[,.(n=.N,dec30=mean(x),d28=mean(d28),hd=mean(hd),proph=mean(proph),pltx=mean(pltx)),by=.(db,x)])
print(D[,lapply(.SD,function(z) round(mean(is.na(z)),2)),by=db,.SDcols=c("inr","ptt","lact","alb","bili","pf","sev","hr","map","creat")])
