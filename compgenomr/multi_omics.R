# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_gex.csv", 
                       package="compGenomRData")
x1 <- read.csv(csvfile, row.names=1)
colnames(x1)
sum(grepl("\\|", rownames(x1)))
knitr::kable(head(t(head(x1))), caption="Example gene expression data (head)")


# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_muts.csv", 
                       package="compGenomRData")
x2 <- read.csv(csvfile, row.names=1)
x2[x2>0] = 1
knitr::kable(head(t(head(x2))), caption="Example mutation data (head)")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_cnv.csv", 
                       package="compGenomRData")
x3 <- read.csv(csvfile, row.names=1)
x3
knitr::kable(head(t(head(x3))), 
             caption="Example copy number data for CRC samples")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_subtypes.csv",
                       package="compGenomRData")
covariates <- read.csv(csvfile, row.names=1)
rownames(covariates) <- gsub(pattern = '-', replacement = '\\.',
                             rownames(covariates))
covariates

covariates = covariates[colnames(x1),]
anno_col <- data.frame(cms=as.factor(covariates$cms_label))
rownames(anno_col) <- rownames(covariates)

library(pheatmap)
pheatmap(x1, show_rownames = F, show_colnames = F, annotation_col = anno_col)
pheatmap(x2, show_rownames = F, show_colnames = F, annotation_col = anno_col)
pheatmap(x3, show_rownames = F, show_colnames = F, annotation_col = anno_col)

##############
## Factorise via PCA
# run the MFA function from the FactoMineR package
r.mfa <- FactoMineR::MFA(
  t(rbind(x1,x2,x3)), # binding the omics types together
  c(dim(x1)[1], dim(x2)[1], dim(x3)[1]), # specifying the dimensions of each
  graph=FALSE)

mfa.h = r.mfa$global.pca$ind$coord
mfa.w <- r.mfa$quanti.var$coord
# create a dataframe with the H matrix and the CMS label
mfa_df <- as.data.frame(mfa.h)
mfa_df$subtype <- factor(covariates[rownames(mfa_df),]$cms_label)

# create the plot
ggplot2::ggplot(mfa_df, ggplot2::aes(x=Dim.1, y=Dim.2, color=subtype)) +
  ggplot2::geom_point() + ggplot2::ggtitle("Scatter plot of MFA")

pheatmap(t(mfa_df[,1:4]), show_colnames = F,
         annotation_col = anno_col)

#####################################
## NMF
x1.featnorm <- x1 / rowSums(x1)
x2.featnorm <- x2 / rowSums(x2)
x3.featnorm <- x3 / rowSums(x3)

x1.featnorm.frobnorm = x1.featnorm / norm(as.matrix(x1.featnorm), type = 'F')
x2.featnorm.frobnorm = x2.featnorm / norm(as.matrix(x2.featnorm), type = 'F')
x3.featnorm.frobnorm = x3.featnorm / norm(as.matrix(x3.featnorm), type = 'F')
compGenomRData::split_neg_columns(t(x3.featnorm.frobnorm))

split_neg_columns <- function(mat){
  df = as.data.frame(mat)
  split_cols= apply(df, 2, function(x){
    pos = x > 0
    neg = x < 0
    x_pos = x
    x_neg = x
    x_pos[neg] = 0
    x_neg[pos] = 0
    x_neg = abs(x_neg)
    col_name_pos = paste0(names(x), "_pos")
    col_name_neg = paste0(names(x), "_neg")
    x_df = data.frame(col_name_pos = x_pos,
                     col_name_neg = x_neg)
    return(x_df)
  })
  
  x_df <- do.call(cbind, split_cols)
  return(x_df)
}

test_split_neg_columns <- function() {
  input <- as.data.frame(cbind(c(1,2,1),c(0,1,-2)))
  output <- as.data.frame(cbind(c(1,2,1), c(0,0,0), c(0,1,0), c(0,0,2)))
  stopifnot(all(output == split_neg_columns(input)))
}
test_split_neg_columns()

x3.featnorm.frobnorm.nonneg <- t(split_neg_columns(t(x3.featnorm.frobnorm)))

require(NMF)
r.nmf <- nmf(t(rbind(x1.featnorm.frobnorm,
                     x2.featnorm.frobnorm,
                     x3.featnorm.frobnorm.nonneg)),
             2,
             method='Frobenius')
r.nmf

# exctract the H and W matrices from the nmf run result
nmf.h <- NMF::basis(r.nmf)
nmf.w <- NMF::coef(r.nmf)
nmfw <- t(nmf.w)

# create a dataframe with the H matrix and the CMS label (subtype)
nmf_df <- as.data.frame(nmf.h)
colnames(nmf_df) <- c("dim1", "dim2")
nmf_df$subtype <- factor(covariates[rownames(nmf_df),]$cms_label)

# create the scatter plot
ggplot2::ggplot(nmf_df, ggplot2::aes(x=dim1, y=dim2, color=subtype)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("Scatter plot of 2-component NMF for multi-omics integration")

pheatmap::pheatmap(t(nmf_df[,1:2]),
                   annotation_col = anno_col,
                   show_colnames=FALSE,
                   main="Heatmap of 2-component NMF")
