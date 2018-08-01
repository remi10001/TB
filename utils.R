# This file contains helpful functions for my R code

save.progress = function(file="Monkey-TimesinceTB-caret-bulmer", dir="/master/rault/TB/data") {
    save.image(file=paste(dir, 
                      paste(file, Sys.Date(), ".RData", sep="-"),
                      sep="/"))
    }


# --------------------------------------------------------------------------------------------------

# Function Definitions for Time Since Infection analysis (the first set of functions is for the mouse data)

overlap = function (x,y) {
  z = rep(F, length(y))
  z[match(x, y)] = T
  return(z)
}

convert.barcode = function(array, conv) {
  conv$Sample_name = as.character(conv$Sample_name)
  for (name in colnames(array)) {
    for (code in conv$Barcode) {
      if (grepl(code, name))
        colnames(array)[colnames(array) == name] = conv[conv$Barcode == code, "Sample_name"]
      # print(paste("Sample name is:", conv[conv$Barcode == code, "Sample_name"]))
    }
  }
  return(colnames(array))
}

# Puts expression and phenotype data from csv format into a Bioconductor ExpressionSet Object
# The FeatureData and PhenoAnnotation Data are not utilized (phenodataAnnotation is currently just transpose of actual phenodata)
# Initial Genes are seleccted based on detection in the microarray data

make.expression.set = function(array_data, pheno, conv) {
  # present in 10% of samples means signal precision <0.01 in 10% of samples
  
  sample_n = 50
  det_pval_thresh = 0.01
  
  # get all the pval columns
  colnames(array_data)
  Detection_Pval = array_data[,grepl("Detection", colnames(array_data))]
  PAL.10  = apply(Detection_Pval <= det_pval_thresh, 1, mean) >= 0.10
  array.Pal.10 = array_data[PAL.10,]
  
  # Want to just remove the detection pVal now that they are irrelevant
  array.Pal.10 = array.Pal.10[,  !grepl("Detection", colnames(array.Pal.10))]
  
  # Changes name of microarray samples from barcode to the sample name in the phenotype data.
  colnames(array.Pal.10) = convert.barcode(array.Pal.10, mousechip_conv)
  expres = array.Pal.10
  rownames(expres) = expres$PROBE_ID
  expres = expres[,grepl("OSU", colnames(expres))]
  reorder = sapply(pheno_data$chip.name, function (x) {which(colnames(expres)==x)})
  expres= expres[,reorder]
  # Grabs only the microarray signal data from the array data frame
  # gene_data = array.Pal.10[,grepl("OSU", colnames(array.Pal.10))]
  pdat = pheno_data[,c("Race", "Condition..Tx.Group", "Time.Point", "chip.name")]
    
  rownames(pdat) = pdat$chip.name
  
  # trimws(
  pdat$Time.Point = as.numeric(  unlist(lapply(strsplit(as.character(pdat$Time.Point), split=" "), function(x) {return(x[[2]])})))
  pdat = pdat[, !(colnames(pdat) == "chip.name")]
  colnames(pdat) = c("Strain", "Infect.Status", "Time.point.days")
    pdat$Infect.Status = as.character(pdat$Infect.Status)
    #print(pdat$Infect.Status)
    pdat$Infect.Status[pdat$Infect.Status != "Mtb"] = "Naive"
    #print(pdat$Infect.Status)
    pdat$Infect.Status = as.factor(pdat$Infect.Status)
    #print(pdat$Infect.Status)
  pdat$chip.name
  colnames(pdat)
  rownames(pdat)
  # What pheno data do I really want for the expression set?
  colnames(pheno_data)
  
  phenoDat = new("AnnotatedDataFrame", data=pdat, varMetadata=as.data.frame(t(pdat)))
  exampleSet = ExpressionSet(assayData=as.matrix(expres), phenoData=phenoDat)
  return(exampleSet)
}



#   Per-gene normalization:
#     Signal intensity of each probe in each sample is divided by 
#     median intensity for that probe across all samples
med.normal = function (gset) {
  meds = apply(exprs(gset), 1, median)
  exprs(gset) = exprs(gset) / meds         # This line of code was inspired by the above website consulted
  return(gset)
}

# Performs the normalization procedure described in Berry et al. 2010 and Mejias et al. 2013
# However, median normalized genes are not saved to the new gene set.
normalization = function(gset, filter = T, med = FALSE) {

  # (1) All signal intensity values less than 10 were set equal to 10
  exprs(gset)[exprs(gset) < 10.0] = 10
  
  # (2) Per-gene normalization by probe median intensity
  med.gset = med.normal(gset)
  
  # (3) Most variable probes:
  #     Select probes that have a minimum of twofold expression change
  #     compared with the median intensity across all samples,
  #     in greater than 10% of all samples
  
  diff.genes = (apply((exprs(med.gset) > 2.00) , 1, mean) > 0.10) | (apply((exprs(med.gset) < 0.50) , 1, mean) > 0.10)
  if (filter)
    gset = gset[diff.genes,]
  
  print(paste(dim(exprs(gset))[1], "genes passed the 2 fold change from median in 10% of samples filter"))
  
# Return gene set with 10 to 10 setting, not median normalized, with only filtered genes

  return(gset)
}

# When less than 10 to 10 was performed:
# "3558 genes passed the 2 fold change from median in 10% of samples filter"
# When less than 10 to 10 not performed
# "4437 genes passed the 2 fold change from median in 10% of samples filter"

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

easy.print = function(data) {
    cat(sapply(colnames(data), function(x) {paste("\"", x, "\",\n", sep="")}))
   }

# --------------------------------------------------------------------------------------------------
