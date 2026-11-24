# 00 extraction: builds the local DuckDB (w.duckdb) used by 01_cohort.R from raw MIMIC-IV v3.1 and eICU-CRD v2.0 csv.gz files
# Set MIMIC_DIR (folder containing hosp/ and icu/) and EICU_DIR before running.
suppressMessages({library(DBI)})
M<-paste0(Sys.getenv("MIMIC_DIR"),"/"); E<-paste0(Sys.getenv("EICU_DIR"),"/")
con<-dbConnect(duckdb::duckdb(shared_home=FALSE),"w.duckdb"); X<-function(s){t<-Sys.time();dbExecute(con,s);cat(round(difftime(Sys.time(),t,units="mins"),1),"min:",substr(s,1,50),"\n")}
X("PRAGMA threads=10"); X("SET memory_limit='18GB'")
rm<-function(p) sprintf("read_csv('%s%s.csv.gz',header=true,auto_detect=true,sample_size=-1)",M,p)
re<-function(p) sprintf("read_csv('%s%s.csv.gz',all_varchar=true,strict_mode=false,ignore_errors=true)",E,p)
# ---------- MIMIC ----------
X(paste0("create or replace table st as select i.subject_id,i.hadm_id,i.stay_id,i.intime,i.outtime,p.gender,p.anchor_age+(year(i.intime)-p.anchor_year) age,p.dod,a.deathtime,a.dischtime,a.hospital_expire_flag,a.race from ",rm("icu/icustays")," i join ",rm("hosp/patients")," p using(subject_id) join ",rm("hosp/admissions")," a using(hadm_id) qualify row_number() over(partition by i.subject_id order by i.intime)=1"))
X("delete from st where age<18")
X(paste0("create or replace table lab as select l.subject_id,s.stay_id,l.specimen_id,l.itemid,l.charttime,l.valuenum v,l.value from ",rm("hosp/labevents")," l join st s on l.subject_id=s.subject_id where (l.itemid=50912 and l.charttime between s.intime-interval 7 day and s.intime+interval 30 day) or (l.itemid in (50818,50820,50803,50882,50813,50862,51006,50971,50983,50902,50931,51222,51301,51265,50868,50821,50885,52033,50816) and l.charttime between s.intime-interval 24 hour and s.intime+interval 30 hour)"))
X(paste0("create or replace table ce as select c.stay_id,c.itemid,c.charttime,c.valuenum v from ",rm("icu/chartevents")," c join st s on c.stay_id=s.stay_id where c.itemid in (226512,224639,220045,220210,224690,220052,220181,225312,220277,223762,223761,220739,223900,223901,223835) and c.charttime between s.intime-interval 6 hour and s.intime+interval 24 hour"))
X(paste0("create or replace table pe as select p.stay_id,p.itemid,p.starttime,p.endtime from ",rm("icu/procedureevents")," p join st s using(stay_id) where p.itemid in (225792,225802,225803,225805,225809,225955,225441,225436)"))
X(paste0("create or replace table inp as select n.stay_id,n.itemid,n.starttime,n.endtime from ",rm("icu/inputevents")," n join st s using(stay_id) where n.itemid in (221906,221289,221662,222315,221749) and n.starttime<=s.intime+interval 24 hour and n.endtime>=s.intime"))
X(paste0("create or replace table dx as select d.hadm_id,d.icd_code,d.icd_version from ",rm("hosp/diagnoses_icd")," d join st s using(hadm_id)"))
# ---------- eICU ----------
X(paste0("create or replace table ept as select * from ",re("patient")))
X(paste0("create or replace table elab as select cast(patientunitstayid as bigint) s,cast(labresultoffset as int) o,labname n,try_cast(labresult as double) v from ",re("lab")," where labname in ('creatinine','paCO2','pH','HCO3','bicarbonate','lactate','albumin','BUN','potassium','sodium','chloride','glucose','Hgb','WBC x 1000','platelets x 1000','anion gap','total bilirubin','paO2','FiO2') and try_cast(labresultoffset as int) between -10080 and 43200"))
X(paste0("create or replace table eaps as select * from ",re("apacheApsVar")))
X(paste0("create or replace table eapr as select cast(patientunitstayid as bigint) s,max(try_cast(apachescore as int)) apache from ",re("apachePatientResult")," where apacheversion='IVa' group by 1"))
X(paste0("create or replace table edx as select cast(patientunitstayid as bigint) s,cast(diagnosisoffset as int) o,lower(diagnosisstring) d from ",re("diagnosis")))
X(paste0("create or replace table eph as select cast(patientunitstayid as bigint) s,lower(pasthistorypath) d from ",re("pastHistory")))
X(paste0("create or replace table etr as select cast(patientunitstayid as bigint) s,cast(treatmentoffset as int) o,lower(treatmentstring) d from ",re("treatment")," where lower(treatmentstring) similar to '.*(dialysis|crrt|cvvh|ultrafiltration|mechanical ventilation|intubation|endotracheal).*'"))
X(paste0("create or replace table erc as select cast(patientunitstayid as bigint) s,try_cast(ventstartoffset as int) vs,try_cast(ventendoffset as int) ve from ",re("respiratoryCare")," where try_cast(ventstartoffset as int)<>0"))
X(paste0("create or replace table einf as select cast(patientunitstayid as bigint) s,cast(infusionoffset as int) o,lower(drugname) d from ",re("infusionDrug")," where cast(infusionoffset as int) between -360 and 1440 and lower(drugname) similar to '.*(norepi|levophed|epineph|dopamine|vasopressin|phenyleph).*'"))
print(dbGetQuery(con,"select table_name,estimated_size from duckdb_tables()")); 
X("create or replace table cm as select distinct stay_id from lab where itemid=50818 and stay_id is not null"); X("create or replace table cei as select distinct s from elab where n='paCO2'")
rm<-function(p) sprintf("read_csv('%s%s.csv.gz',header=true,auto_detect=true,sample_size=200000,ignore_errors=true)",M,p)
re<-function(p) sprintf("read_csv('%s%s.csv.gz',all_varchar=true,strict_mode=false,ignore_errors=true)",E,p)
X(paste0("create or replace table ce2 as select c.stay_id,c.itemid,c.charttime,c.value,c.valuenum from ",rm("icu/chartevents")," c join st s on c.stay_id=s.stay_id where c.stay_id in (select stay_id from cm) and c.itemid in (223900,226732,223849,220339,224684,224685,224695,224696,224686,229314,227287) and c.charttime between s.intime-interval 30 hour and s.intime+interval 30 hour"))
X(paste0("create or replace table pe2 as select p.stay_id,p.itemid,p.starttime,p.endtime from ",rm("icu/procedureevents")," p where p.stay_id in (select stay_id from cm) and p.itemid in (225794,225792)"))
X(paste0("create or replace table erc2 as select cast(patientunitstayid as bigint) s,try_cast(respchartoffset as int) o,respchartvaluelabel l,respchartvalue v from ",re("respiratoryCharting")," where cast(patientunitstayid as bigint) in (select s from cei) and try_cast(respchartoffset as int) between -1800 and 1800"))
X(paste0("create or replace table evp as select cast(patientunitstayid as bigint) s,cast(observationoffset as int) o,try_cast(heartrate as double) hr,try_cast(systemicmean as double) map from ",re("vitalPeriodic")," where cast(patientunitstayid as bigint) in (select s from cei) and try_cast(observationoffset as int) between -1800 and 1800"))
X(paste0("create or replace table evap as select cast(patientunitstayid as bigint) s,cast(observationoffset as int) o,try_cast(noninvasivemean as double) map from ",re("vitalAperiodic")," where cast(patientunitstayid as bigint) in (select s from cei) and try_cast(observationoffset as int) between -1800 and 1800"))
X(paste0("create or replace table enc2 as select cast(patientunitstayid as bigint) s,cast(nursingchartoffset as int) o,nursingchartcelltypevalname n,nursingchartvalue v from ",re("nurseCharting")," where cast(patientunitstayid as bigint) in (select s from cei) and try_cast(nursingchartoffset as int) between -1800 and 1800 and nursingchartcelltypevalname in ('GCS Total','Verbal','O2 Admin Device')"))
print(dbGetQuery(con,"select table_name,estimated_size from duckdb_tables() where table_name in ('ce2','pe2','erc2','evp','evap','enc2')")); dbDisconnect(con,shutdown=TRUE)

dbDisconnect(con,shutdown=TRUE)
