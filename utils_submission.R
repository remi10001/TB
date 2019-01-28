# This file contains helpful functions for my R code

# Function Definitions for paper analysis reproduction

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




library(lubridate)


save.progress = function(file="Monkey-TimesinceTB-caret-bulmer", dir="/master/rault/TB/data") {
    save.image(file=paste(dir, 
                      paste(file, Sys.Date(), ".RData", sep="-"),
                      sep="/"))
    }


# --------------------------------------------------------------------------------------------------

# Function Definitions for Time Since Infection analysis (the first set of functions is for the mouse data)




# I want to reproduce what I created
graph.PCA = function(PCA, strat, pheno) {
  # pheno and PCA may have different samples. I want to remove items from pheno
  # It may be the case that I have errors with retained factors, could deal with these by converting pheno
  # to all character at beginning if I have to. That would be the easiest fix if I have this problem
  # print(paste("Old dimensions of pheno table:", dim(pheno)))
  # retain = overlap(rownames(PCA), pheno$chip.name)
  # pheno = pheno[retain,]
  # print(paste("New dimensions of pheno table:", dim(pheno)))
  
  # How many colors will I need
  cols = rainbow(length(table(pheno[,strat])))
  # plot.col = rep("blue", length(rownames(PCA)))
  # I need to assign the different variables in pheno stat
  strat.names = names(table(pheno[,strat]))
  # can just use which on the names to get the index of the color
  plot.col = rep("blue", length(rownames(PCA)))
  
  # I need to go through this plot.col with each sample, find the samples category and change color if necessary
  i = 1
  for (sample in rownames(PCA)) {
    plot.col[i] = cols[strat.names == pheno[rownames(pheno) == sample,strat]]
    i = i + 1
  }
  for (col in plot.col)
    print(col)
  
  # should modify the code to choose PCA axes
  x <- PCA[,4]
  y <- PCA[,5]
  z <- PCA[,6]
  plot3d(x, y, z, col=plot.col, type='p', size=10.0, lwd=1,
         xlab="PC4", ylab="PC5", zlab="PC6",
         
         
         main="Principle Component Analysis of Mouse Microarray Data")
  legend3d("topright", legend = strat.names, pch = 16, col = cols, cex=1, inset=c(0.02))
  
  # print("Here are the plot colors", plot.col)
  # we have the new pheno, it is all in order
  
}

# Here are functions for the monkey data analysis
process.vars = function(pheno, covars, num.vars, fac.vars) {
  
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



process.data = function(expres, PAL) {
  
  expres = expres[PAL,]
  # All signal intensity values less than 10 are set equal to 10
  expres[expres < 10.0] = 10
  
  # Log 2 transform the data
  expres = log2(expres)
  
  
  return(expres)
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
    
    # Depending on the analysis, I may want time.from.exposure.months to be a factor, but for now is a numeric
    numbers = c("age",
                 "time.from.exposure.months",
                 "time.to.tb.months")
    

    
    new.data = data[,cols]
    colnames(new.data) = new.names
    for (col in factors)
        new.data[,col] = as.factor(new.data[,col])
    for (col in numbers)
        new.data[,col] = as.numeric(new.data[,col])
    
 
    # Filter out subjects with NA on clinical data. There are 6 of them.
    # Actually, we need to add them back in! It's a typo.
    
    # I need to email them again about all the samples that are missing... Just go off expression codes.
    # Also tell them that they have errors in their uploaded pheno table matrix file.
    # I should just check the totals too.
    

    
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
    # I have emailed Gerhard Werlzl and Daniel Zak, the corresponding authors, asking about this
    
    dup.codes = sort(as.numeric(names(table(new.data$code)[table(new.data$code) == 2])))
    dup.pheno = new.data[new.data$code %in% dup.codes, ]
    dup.pheno <- dup.pheno[order(dup.pheno$code),] 
    print(rownames(dup.pheno))
    print(dup.pheno$code)

    gsm.to.remove = rownames(dup.pheno)[seq(1, length(dup.pheno$code), by=2)]
    
    new.data = new.data[!(rownames(new.data) %in% gsm.to.remove),]
    
    new.data = new.data[order(as.numeric(as.character(new.data$code))),]
    
    # Add training and test split from Suliman et al
    
    control_data = read.csv("data/Suliman_et_al_Human_Data//Additional_data_table_trainingtest_controls.csv", header=F)
    prog_data  = read.csv("data/Suliman_et_al_Human_Data//Additional_data_table_trainingtest_progressors.csv", header=F)

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
    
    # the one AHRI should be training
    new.data$dataset[new.data$site == "AHRI"] = "Test"
    
    new.data$dataset = as.factor(new.data$dataset)
    
    # Go ahead and remove the few UGandan samples:
    new.data = droplevels(new.data[new.data$site != "UGA",])
    
    #new.data$tb.status = as.factor(ifelse(is.na(new.data$time.to.tb.months), "control", "prog"))
    
    print("about to return new data frame")
    return(new.data)
}

filter.human.exprs = function(exprs, pheno, splice=F) {
    sample.start = 2
    if (splice)
        sample.start = 7
    exprs.cols = colnames(exprs)[sample.start:dim(exprs)[2]]
    exprs.cols = gsub("X", "", exprs.cols)
   
    
    # Some of the double genes are not labeled correctly. For example, ENSG00000124191 is TOX2, not KLRD1. Also ENSG00000260539 no longer maps to a gene. I am going to go ahead and stick with the ENSEMBL identifier. I'll go back to the genes on the interesting genes, etc.
    
    # I will just remove the gene symbol for now.
    exprs = exprs[,-c(1:(sample.start-1))]
    colnames(exprs) = exprs.cols
    
    # I also need to remove codes not in the pheno data, those 6 samples that didn't have clinical data.
    exprs = exprs[,exprs.cols %in% pheno$code]
    
    # Sort columns by code
    exprs = exprs[, order(as.numeric(as.character(colnames(exprs))))]
    
    
     print("Are all the codes in the expression data in the same order as the codes in pheno?")
    print(identical(as.character(colnames(exprs)), as.character(pheno$code)))
    #print(data.frame(exprs=colnames(exprs)[1:50], pheno=pheno$code[1:50]))
    # Remove gene symbol
    return(exprs)
}

graph_PCA = function(PCA, n_dim, var, pheno) {
    gg.PCA = as.data.frame(PCA$x[,1:n_dim])
    gg.PCA$time.period = pheno[,var]
    gg.PCA$time.period = as.character(gg.PCA$time.period)
    options(repr.plot.width=8, repr.plot.height=5)
    print(ggpairs(gg.PCA, aes(colour = time.period)))
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
    
    
    # Add in a time since exposure variable, which is time since baseline until Anne O'Garra provides me with further information.
    time.since.exposure.days = rep(0, dim(contact.pheno)[1])

    for (sample in 1:dim(contact.pheno)[1]) {
        time.since.exposure.days[sample] = ymd(contact.pheno[sample, "visit_date"]) - 
                                    ymd(filter(contact.pheno, 
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

library(verification)

roc.pvalue = function(obs, prob.pred, pos) {
    return(roc.area(ifelse(as.factor(obs)==pos, 1, 0), prob.pred)$p.value)
}

debug.my.roc = function(pred.prob, obs, pos, title="ROC Curve") {
    print("WHAT FOLLOWS HEREAFTER IS THE CORRECT ORDER OF FACTORS")
    levs = levels(as.factor(obs))
    print("Levels before reversing")
    print(levs)
    if (length(levs) != 2) {
        print("roc requires only two classes!")
        return(-1)
    }
    
    if (levs[2] != pos)
        levs = rev(levs)
    
    print("Levels after reversing")
    print(levs)
    
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
    
    print("WHAT FOLLOWS HEREAFTER IS THE OPPOSITE, INCORRECT ORDER OF FACTORS")
    pos = levs[which(levs != pos)]
    levs = levels(as.factor(obs))
    print("Levels before reversing")
    print(levs)
    if (length(levs) != 2) {
        print("roc requires only two classes!")
        return(-1)
    }
    
    # Here is the changed line for debug
    if (levs[2] != pos)
        levs = rev(levs)
    print("Levels after reversing")
    print(levs)
    
    the.roc.wrong = roc(predictor=pred.prob,
              response = as.factor(obs),
              levels=levs)
    print("This is the AUC:")
    print(the.roc.wrong$auc)
    
    print("This is the AUC p-value:")
    print(roc.pvalue(obs, pred.prob, pos ))
    print("This is the AUC 95% Confidence Interval")
    print(ci(the.roc.wrong))
    plot.roc(the.roc.wrong, main=title, legacy.axes=T)
    
    return(the.roc) # this is the correct roc
}

    
    
# pred.prob is the probability of the pos class
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

library(akima)

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
    
    # Depending on the analysis, I may want time.from.exposure.months to be a factor, but for now is a numeric
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
    #q = qplot(obs, pred, data= data, geom = c("point", "smooth"), method="loess")
    if (log) {
        data$obs = 2 ^ data$obs
        data$pred = 2 ^ data$pred
    }
    
    break_points = c(0, 30, 60, 90 , 120, 150, 180)
    if (break.90)
        break_points = c(0, 20, 30, 42, 56, 90)
    q = ggplot(data, aes(obs, pred, group=obs)) +
        
        scale_y_continuous(breaks=break_points, limits = c(min(break_points), max(break_points)+15)) +
        scale_x_continuous(breaks=break_points, limits = c(min(break_points), max(break_points)+15)) +
        geom_point() + geom_smooth(method="loess") +
        geom_boxplot(outlier.shape=NA, colour = "red", fill="white", width=10) +
        theme_bw() +
        labs(x="Days Post Infection", y="Predicted Days Post Infection") + 
        ggtitle(paste(label)) +  #, "regression in cynomolgus macaques")) +
        theme(plot.title = element_text(hjust = 0.5))
    
    return(q)
}

get.GEO.eset = function(old_array, c_array, pheno, conv, raw=F) {
    # Several Steps
    # (1) Separate signal and p-value, then remove from all arrays not in the pheno (i.e. not C57BL/6, since only B6 is being published)
    # (2) Order all arrays by their order in the pheno (reorder, then convert to new name)
    # (3) Reorder the rows according to the order in the old array.
    # (4) Put together the ID_REF, signal and p-value pieces.
    
    
    # (1) Separate signal and p-value, then remove from all arrays not in the pheno (i.e. not C57BL/6, since only B6 is being published)
    ID_REF = as.character(c_array$ID_REF)
    p_val = c_array[,grepl("Detection", colnames(c_array))]
    if (raw) p_val = 1 - p_val
    if (!raw) {
         expres = c_array[,!grepl("Detection", colnames(c_array))]
           expres = expres[,2:dim(expres)[2]] # remove ID_REF
        } else {
        expres = c_array[,grepl("AVG_Signal", colnames(c_array))]
        colnames(expres) = gsub("AVG_Signal.", "", colnames(expres), fixed=T)
        
        }
   
    colnames(expres) = convert.barcode(expres, conv)
    select_samples = colnames(expres) %in% rownames(pheno)
    p_val = p_val[,select_samples]
    expres = expres[,select_samples]
    
    
    # (2) Order all arrays by their order in the pheno (reorder, then convert to new name)
    reorder = match(rownames(pheno), colnames(expres))
    p_val = p_val[,reorder]
    expres = expres[,reorder] 
    print('Are pheno rows the same as expression columns?')
    print(identical(rownames(pheno), colnames(expres)))
    
    colnames(expres) = pheno$Mouse.ID
    
    # (3) Reorder the rows according to the order in the old array.
    old_ID_REF = as.character(old_array$PROBE_ID)
    row_reorder = match(old_ID_REF, ID_REF)
    ID_REF = ID_REF[row_reorder]
    p_val = p_val[row_reorder,]
    expres = expres[row_reorder,]
    
    # (4) Put together the ID_REF, signal and p-value pieces.
    new_array = data.frame(ID_REF = ID_REF)
    name_c = 2
    for (i in 1:dim(expres)[2]) {
        new_array = cbind(new_array, expres[,i])
        colnames(new_array)[name_c] = colnames(expres)[i]
        name_c = name_c + 1
        new_array = cbind(new_array, p_val[,i])
        colnames(new_array)[name_c]= "Detection Pval"
        name_c = name_c + 1
    }
    
    return(new_array)
    
}

get.GEO.pheno = function(pheno, conv) {
    pdat = pheno[,c("Mouse.ID", "Race", "Condition..Tx.Group", "Time.Point", "Tissue.Type", "chip.name")]

    rownames(pdat) = pdat$chip.name

    # trimws(
    pdat$Time.Point = as.numeric(  unlist(lapply(strsplit(as.character(pdat$Time.Point), split=" "), function(x) {return(x[[2]])})))
    pdat = pdat[, !(colnames(pdat) == "chip.name")]
    colnames(pdat) = c("Mouse.ID", "Strain", "Infect.Status", "Time.point.days", "Tissue")
    pdat$Infect.Status = as.character(pdat$Infect.Status)
    #print(pdat$Infect.Status)
    pdat$Infect.Status[pdat$Infect.Status != "Mtb"] = "Naive"
    #print(pdat$Infect.Status)
    pdat$Infect.Status = as.factor(pdat$Infect.Status)
    pdat$Gender = "female"


    pdat = pdat[pdat$Strain == "C57BL/6",]
    pdat$Mouse.ID = paste(pdat$Infect.Status, pdat$Time.point.days, pdat$Mouse.ID, sep="_")
    return(pdat)
}

# --------------------------------------------------------------------------------------------------
