# This file contains helpful functions for my R code

# Function Definitions for paper analysis reproduction:

library(dplyr)

process.monkey.vars = function(pheno, covars, num.vars, fac.vars) {
  
  # Select only variables we are interested in
  pheno = pheno[,covars]
  
  # Remove field names from data and convert categorical variables to factors, numeric to numeric
  for (var in union(num.vars, fac.vars)) {
    old.data = as.character(pheno[,var])
    new.name = strsplit(old.data[1], split=":")[[1]][1]
    new.data = as.character(sapply(old.data, function(x) {return(trimws(strsplit(x, split=":")[[1]][2]))}))
    if (var %in% num.vars)
      new.data = as.numeric(new.data)
    else
      new.data = as.factor(new.data)
    
    colnames(pheno)[colnames(pheno) == var] = new.name
    pheno[,new.name] = new.data
  }
  
  # Replace spaces with periods in variable names to facilitate R subsetting
  colnames(pheno) = gsub(" ", ".", colnames(pheno), fixed=T)
  
  # The clinical status of monkey 18 at time point 0 is incorrectly labeled as latent. Should be active. This is an error in the data submission to GEO.
  pheno$clinical.status[pheno$monkeyid == "M18" & pheno$time.point==0] = "Active"
  
  # Create an extra numeric variable to combine the baseline into one time point
  pheno$time.point.comb = ifelse(pheno$time.point==1, 0, pheno$time.point)
  
  return(pheno)
}

filter.monkey.pheno = function(pheno) {
    covars = c("title",
  "characteristics_ch1.1",
  "characteristics_ch1.2",
   "characteristics_ch1.3",
  "characteristics_ch1.4",
  "characteristics_ch1.6",
  "characteristics_ch1.7",
  "characteristics_ch1.8",
  "characteristics_ch1.9",
  "description",
  "description.1")

  num.vars = c("characteristics_ch1.7")

  fac.vars = c(
  "characteristics_ch1.1",
  "characteristics_ch1.2",
   "characteristics_ch1.3",
  "characteristics_ch1.4",
  "characteristics_ch1.6",
  "characteristics_ch1.8",
  "characteristics_ch1.9")

  pheno.f = process.monkey.vars(pheno, covars, num.vars, fac.vars)
  print(pheno.f)
    return(pheno.f)
}

get.monkey.expressed.genes = function(raw.expres) {
    det_pval_thresh = 0.01

  # 13 Time points
  # 100/13 = 7.69 %
  # Choose 5% to get equivalent of significant part of one time point
  percent_samples= 0.05
  Detection_Pval = raw.expres[,grepl("Detection", colnames(raw.expres))]
  
  PAL.5 = raw.expres$ID_REF[apply(Detection_Pval <= det_pval_thresh, 1, mean) >= percent_samples]
  
  print(paste("Genes expressed in at least", percent_samples*100, "% of samples:", length(PAL.5)))
    
  return(PAL.5)
}

process.monkey.exprs = function(expres, PAL) {
  
  expres = expres[PAL,]
  # All signal intensity values less than 10 are set equal to 10
  expres[expres < 10.0] = 10
  
  # Log 2 transform the data
  expres = log2(expres)
  
  
  return(expres)
}


if (!require("lubridate")) {
  install.packages("lubridate")
  library("lubridate")
}


make.train.test = function(data, pheno) {
    for (p in 1:dim(data)[1]) {
        for (d in 2:dim(data)[2])
            if (data[p,d] %in% c("Training", "Test")) {
                idx = which(pheno$subjectid == as.character(data[p, 1]))
                if (length(idx) == 0) {
                    print(as.character(data[p,1]))
                    print(data[p,d])
                    }
                pheno$dataset[idx] = as.character(data[p,d])
                }
        }
    return(pheno)
}

filter.human.pheno = function(data) {
    cols = c("age:ch1",
             "code:ch1",
             "gender:ch1",
             "group:ch1",
             "site:ch1",
             "subjectid:ch1",
             "time.from.exposure.months:ch1",
             "time.to.tb.months:ch1")
     
    new.names = c("age",
                 "code",
                 "gender",
                 "group",
                 "site",
                 "subjectid",
                 "time.from.exposure.months",
                 "time.to.tb.months")
    
    factors = c("code","gender",
                 "group",
                 "site",
                 "subjectid")
    
    numbers = c("age",
                 "time.from.exposure.months",
                 "time.to.tb.months")
    

    new.data = data[,cols]
    colnames(new.data) = new.names
    for (col in factors)
        new.data[,col] = as.factor(new.data[,col])
    for (col in numbers)
        new.data[,col] = as.numeric(new.data[,col])
    
    # Six (6) subjects have NA on their clinical data in the GEO series of Suliman et al, but the paper's supplementary tables have their information. I have added them here below.
    
    addit.row.names = c("GSM2475704", "GSM2475705", "GSM2475706", "GSM2475722", "GSM2475742", "GSM2475748")
    add.age = c(32, 34, 28, 40, 24, 18)
    add.code = c(672, 694, 695, 994, 1061, 1194)
    add.gender = c("F", "F", "F", "F", "F", "F")
    add.group = c("case (TB)", "Control", "Control", "Control", "Control", "Control")
    add.site = c("UGA", "UGA", "UGA", "AHRI", "AHRI", "AHRI")
    add.subjectid = c("92245", "91420103",  "91451104", "KHHC87", "DZHHC06", "KAZHHC23")
    add.time.from.exposure.months = c(0, 0, 0, 18, 6, 18)
    add.time.to.tb.months = c(6, NA, NA, NA, NA, NA)
    
    add.data = data.frame(age = add.age, code = as.factor(add.code), gender = as.factor(add.gender),
                          group = as.factor(add.group), site = as.factor(add.site), subjectid = as.factor(add.subjectid),
                          time.from.exposure.months = add.time.from.exposure.months,
                          time.to.tb.months = add.time.to.tb.months)
    rownames(add.data) = addit.row.names
    
    print(new.data[new.data$gender == "NA",])
    new.data = new.data[new.data$gender != "NA",]
    new.data = droplevels(new.data)
    new.data = rbind(new.data, add.data)
    

    # Filter out the duplicated GEO Access Numbers (16 samples, have same code and patient id and clinical info)
    # These are simply from sequencing files being split into 2 files for these samples (correspondance with Gerhard Walzl and Daniel Zak)
    
    dup.codes = sort(as.numeric(names(table(new.data$code)[table(new.data$code) == 2])))
    dup.pheno = new.data[new.data$code %in% dup.codes, ]
    dup.pheno <- dup.pheno[order(dup.pheno$code),] 
    print(rownames(dup.pheno))
    print(dup.pheno$code)

    gsm.to.remove = rownames(dup.pheno)[seq(1, length(dup.pheno$code), by=2)]
    
    new.data = new.data[!(rownames(new.data) %in% gsm.to.remove),]
    
    new.data = new.data[order(as.numeric(as.character(new.data$code))),]
    
    # Add training and test split from Suliman et al. I collected these data from the Supplemental tables of Suliman et al.
    
    control_data = read.csv("data/Suliman_et_al_Additional_Data/Additional_data_table_trainingtest_controls.csv", header=F)
    prog_data  = read.csv("data/Suliman_et_al_Additional_Data/Additional_data_table_trainingtest_progressors.csv", header=F)

    # Match subjectids to how they are in pheno table
    control_data$V1 =gsub("20([0178923]*)/", "\\1/", control_data$V1)
    prog_data$V1 =gsub("20([0178923]*)/", "\\1/", prog_data$V1)
    
    new.data$dataset = NA
    
    new.data = make.train.test(control_data, new.data)
    new.data = make.train.test(prog_data, new.data)

    # Make those not assigned a set in Suliman et al Training. These are only 2-3 subjects
    print("Subjects not assigned a training-test set that are now training")
    print(table(droplevels(new.data$subjectid[is.na(new.data$dataset)])))
    new.data$dataset[is.na(new.data$dataset)] = "Training"
    
    # the one AHRI labeled training should be test, according to Suliman et al.
    new.data$dataset[new.data$site == "AHRI"] = "Test"
    
    new.data$dataset = as.factor(new.data$dataset)
    
    # Go ahead and remove the few Ugandan samples excluded from the Suliman et al study:
    new.data = droplevels(new.data[new.data$site != "UGA",])
    
    
    print("about to return new data frame")
    return(new.data)
}

filter.human.exprs = function(exprs, pheno, splice=F) {
    sample.start = 2
    if (splice)
        sample.start = 7
    exprs.cols = colnames(exprs)[sample.start:dim(exprs)[2]]
    exprs.cols = gsub("X", "", exprs.cols)
   
    # Remove the gene symbol
    exprs = exprs[,-c(1:(sample.start-1))]
    colnames(exprs) = exprs.cols
    
    # remove any sample codes not in the pheno data.
    exprs = exprs[,exprs.cols %in% pheno$code]
    
    # Sort columns by code
    exprs = exprs[, order(as.numeric(as.character(colnames(exprs))))]
    
    print("Are all the codes in the expression data in the same order as the codes in pheno?")
    print(identical(as.character(colnames(exprs)), as.character(pheno$code)))

    return(exprs)
}


detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

easy.print = function(data, qte=T) {
    sepr = ifelse(qte, "\"", "") 
    
    cat(sapply(colnames(data), function(x) {paste(sepr, x, sepr ,",\n", sep="")}))
   }

easy.v.print = function(x) {cat(sapply(x, function(x) {paste("\"", x, "\",\n", sep="")}))}

calculate.distribution = function(data) {
    data.tadd = data
    data.tadd$relative.time.from.conversion.days[data.tadd$relative.time.from.conversion.days == 0] = 1
    table(data$relative.time.from.conversion.days)
    print(table(tapply(data.tadd$relative.time.from.conversion.days, data.tadd$patientid, sum)))

}

filter.Garra.pheno = function(contact.pheno) {
   
    colnames(contact.pheno) = gsub(" ",".", gsub(":ch1|\\(|\\)","",colnames(contact.pheno)))
    
    contact.pheno = contact.pheno[,c("title", "source_name_ch1",
                                 "age_at_baseline_visit",
 "birth_place",
 "ethnicity",
 "gender",
 "group",
 "outlier",
 "patient_id",
 "smear_result",
 "tb_disease_type",
 "timepoint_months",
 "uk_arrival_year",
 "visit_date")]

fact.f = c(
 "birth_place",
 "ethnicity",
 "gender",
 "group",
 "outlier",
 "patient_id",
 "smear_result",
 "tb_disease_type",
 "timepoint_months")

num.f = c  ("age_at_baseline_visit",
 "uk_arrival_year")
        
    for (f in fact.f) {
    contact.pheno[,f] = as.factor(contact.pheno[,f])
}

for (f in num.f) {
    contact.pheno[,f] = as.numeric(contact.pheno[,f])
}
    
    
    # Add in a time since exposure variable, which is time since baseline
    time.since.exposure.days = rep(0, dim(contact.pheno)[1])

    for (sample in 1:dim(contact.pheno)[1]) {
        time.since.exposure.days[sample] = ymd(contact.pheno[sample, "visit_date"]) - 
                                    ymd(dplyr::filter(contact.pheno, 
                                               patient_id == contact.pheno[sample, "patient_id"], 
                                               timepoint_months == "Baseline")$visit_date[1])
    }
    contact.pheno$time.since.exposure.days = as.numeric(time.since.exposure.days)
    
    print(summary(contact.pheno))
    return(contact.pheno)
          
}

rmse <- function(truth, pred) 
    sqrt(mean((truth - pred)^2))

mea <- function(truth, pred)
    mean(abs(truth-pred))

mea.med <- function(truth, pred)
    median(abs(truth-pred))

all.metrics = function(truth, pred) {
    
    print("Root Mean Squared Error (RMSE)")
    print(rmse(truth, pred))
    print("Mean Absolute Error")
    print(mea(truth, pred))
    print("Median Absolute Error")
    print(mea.med(truth, pred))
    print("Pearson Correlation Coefficient")
    print(cor(truth, pred))
    print("Pearson Correlation Signifcance Test")
    print(cor.test(truth, pred))
    print("Spearman Correlation Signifcance Test")
    print(cor.test(truth, pred, method = 'spearman'))
    print("R squared")
    print(caret::R2(pred, truth))
}

if (!require("verification")) {
  install.packages("verification")
  library("verification")
}


roc.pvalue = function(obs, prob.pred, pos) {
    return(roc.area(ifelse(as.factor(obs)==pos, 1, 0), prob.pred)$p.value)
}
    

# pred.prob is the probability of the positive class
my.roc = function(pred.prob, obs, pos, title="ROC Curve") {
    
    levs = levels(as.factor(obs))
    if (length(levs) != 2) {
        print("roc requires only two classes!")
        return(-1)
    }
    
    if (levs[2] != pos)
        levs = rev(levs)
    
    the.roc = roc(predictor=pred.prob,
              response = as.factor(obs),
              levels=levs)
    print("This is the AUC:")
    print(the.roc$auc)
    
    print("This is the AUC p-value:")
    print(roc.pvalue(obs, pred.prob, pos ))
    print("This is the AUC 95% Confidence Interval")
    print(ci(the.roc))
    plot.roc(the.roc, main=title, legacy.axes=T)
    return(the.roc)
}

my.roc.auc = function(pred.prob, obs, pos, title="ROC Curve") {
    
    levs = levels(as.factor(obs))
    if (length(levs) != 2) {
        print("roc requires only two classes!")
        return(-1)
    }
    
    if (levs[2] != pos)
        levs = rev(levs)
    
    the.roc = roc(predictor=pred.prob,
              response = as.factor(obs),
              levels=levs)

    return(the.roc$auc)
}

if (!require("akima")) {
  install.packages("akima")
  library("akima")
}

graph.hyper = function(x, y, z) {
   interpdf <-interp2xyz(interp(x=x, y=y, z=z, duplicate="mean"), data.frame=TRUE)

interpdf %>%
  dplyr::filter(!is.na(z)) %>%
  tbl_df() %>%
  ggplot(aes(x = x, y = y, z = z, fill = z)) + 
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.05) + 
  scale_fill_distiller(palette="Spectral", na.value="white") + 
  theme_bw() 
}

filter.ACS.pheno = function(data) {
    cols = c("title",
             "age:ch1",
             "bin:ch1",
             "gender:ch1",
             "ethnicity:ch1",
             "previous diagnosis of tb:ch1",
             "group:ch1",
             "qft:ch1",
             "tst:ch1")
    
    
    new.names  = c("title",
                 "age",
             "bin",
             "gender",
             "ethnicity",
             "previous.diagnosis.of.tb",
             "group",
             "qft",
             "tst")
    factors = c("title",
             "bin",
             "gender",
             "ethnicity",
             "previous.diagnosis.of.tb",
             "group",
             "qft")
    
    numbers = c("age",
                 "tst")
    
    
    new.data = data[,cols]
    colnames(new.data) = new.names
    for (col in factors)
        new.data[,col] = as.factor(new.data[,col])
    for (col in numbers)
        new.data[,col] = as.numeric(new.data[,col])
    
    print("about to return new data frame")
    return(new.data)
}

filter.HUMAN.exprs = function(exprs, pheno) {
    exprs = exprs[, colnames(exprs) %in% row.names(pheno)]
    exprs = exprs[, match(row.names(pheno), colnames(exprs) )]
    
    print("Identical column and rownames between exprs and pheno tables?")
    print(identical(colnames(exprs), row.names(pheno)))
    return(exprs)
}   

generate.regres.graph = function(data, label, log = F, break.90 = F) {
    if (log) {
        data$obs = 2 ^ data$obs
        data$pred = 2 ^ data$pred
    }
    
    label.y = 180
    label.x = 75
    
    if (break.90) {
    label.x = 42
    label.y = 90
    }
    
    
    break_points = c(0, 30, 60, 90 , 120, 150, 180)
    if (break.90)
        break_points = c(0, 20, 30, 42, 56, 90)
    q = ggplot(data, aes(obs, pred)) + #, group=obs)) +
        
        scale_y_continuous(breaks=break_points, limits = c(min(break_points), max(break_points)+15)) +
        scale_x_continuous(breaks=break_points, limits = c(min(break_points), max(break_points)+15)) +
        geom_point() + 
        geom_boxplot(aes(obs, pred, group=obs), outlier.shape=NA, colour = "red", fill="white", width=10) +
        geom_smooth(method = "lm", col= "blue", se = FALSE, size=0.5) + 
        stat_cor(method = "pearson", label.x = label.x, label.y = label.y) + 
        theme_bw() +
        labs(x="Days Post Infection", y="Predicted Days Post Infection") + 
        ggtitle(paste(label)) + 
        theme(plot.title = element_text(hjust = 0.5))
    
    return(q)
}

posGenes = c("RP11-552F3.12", "PYURF" ,        "TRIM7"  ,       "TUBGCP4"  )
negGenes = c("ZNF608", "BEAN1")

# Calculates arithmetic mean of log2 gene values of upregulated genes - downregulated genes
# exprs is an expression matrix of genes with nsamples as rows and ngenes as columns
calculateScore = function(exprs, posGenes, negGenes) {
    pos.exprs = exprs[,colnames(exprs) %in% posGenes]
    neg.exprs = exprs[,colnames(exprs) %in% negGenes]
    
    if (!is.null(dim(neg.exprs))) {
        scores = apply(pos.exprs, 1, mean) - apply(neg.exprs, 1, mean)
        } else {
        scores = apply(pos.exprs, 1, mean) - neg.exprs
    }
    return (scores)
    }