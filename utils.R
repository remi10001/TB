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
normalization = function(gset, med = FALSE) {

  # (1) All signal intensity values less than 10 were set equal to 10
  exprs(gset)[exprs(gset) < 10.0] = 10
  
  # (2) Per-gene normalization by probe median intensity
  med.gset = med.normal(gset)
  
  # (3) Most variable probes:
  #     Select probes that have a minimum of twofold expression change
  #     compared with the median intensity across all samples,
  #     in greater than 10% of all samples
  
  diff.genes = (apply((exprs(med.gset) > 2.00) , 1, mean) > 0.10) | (apply((exprs(med.gset) < 0.50) , 1, mean) > 0.10)
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

# --------------------------------------------------------------------------------------------------
