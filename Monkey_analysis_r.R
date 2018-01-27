# This post answers my question on using Jupyter notebook with R
#https://discuss.analyticsvidhya.com/t/how-to-run-r-on-jupyter-ipython-notebooks/5512

# Run my personal libraries
source("/master/rault/R-packages/Rconfigure_default.R")

# Set-up environment

if (!require(randomForest)) {
  install.packages("randomForest")
  library(randomForest)
}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require("glmnet")) {
  install.packages("glmnet")
  library("glmnet")
}

if (!require("caret")) {
  install.packages("caret")
  library("caret")
}

data.dir = "/master/rault/TB/data"
script.dir = "/master/rault/TB"

# Read in the data
pheno = read.table(file = paste(data.dir, "Monkey_PhenoData_middle-late.txt", sep="/"), sep="\t", header=T)
expres = read.table(file = paste(data.dir, "Monkey_Processed_ExpressionData_middle-late.txt", sep="/"), header=T, sep="\t")

