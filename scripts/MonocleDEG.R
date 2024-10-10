#### Revised version of how Monocle runs DE, to produce Log2FC values and also shrinkage of Log2FC using apeglm
monocle_diff <- function(cds, cells = NULL, subtype = NULL, num_cells = 5, min.pct = 0.5, full_formula, reduced_formula, foldChange = 'mle', apeglm_method = 'nbinomC'){
  cds0 <- cds
  
  if(!is.null(cells)){
    cds0 <- cds0[,cells]
  }
  
  cds0 <- detectGenes(cds0, 0)
  
  cds0 <- cds0[fData(cds0)$num_cells_expressed >= num_cells,]
  genes <-  c()
  if(!is.null(subtype)){
    if(length(subtype) == 1){
      pData(cds0)[pData(cds0)$cellType != subtype,]$cellType <- 'other'
      subtype <- c(subtype, 'other')
    }else{
      cds0 <- cds0[,row.names(subset(pData(cds0), cellType %in% subtype))]
    }
    for(s in subtype){
      ct <- exprs(cds0[,row.names(subset(pData(cds0), cellType == s))])
      print(dim(ct))
      print(round(ncol(ct)*min.pct))
      genes <- union(genes, row.names(ct[rowSums(ct > 1)>= min(num_cells, round(ncol(ct)*min.pct)),]))
      print(length(genes))
    }
    
  }else{
    genes <- row.names(fData(cds0))
  }
  
  cds0 <- cds0[genes,]
  if(grepl('num_genes_expressed', full_formula)){
    pData(cds0)$num_genes_expressed <- scale(pData(cds0)$num_genes_expressed)
  }
  cds0 <- estimateSizeFactors(cds0)
  cds0 <- estimateDispersions(cds0)
  ptm <- proc.time()
  #diff_test <- differentialGeneTest(cds = cds0, fullModelFormulaStr = full_formula, reducedModelFormulaStr = reduced_formula)
  diff_test_res <- smartEsApply(cds0,1,diff_test_helper,
                                convert_to_dense=TRUE,
                                fullModelFormulaStr=full_formula,
                                reducedModelFormulaStr=reduced_formula, 
                                expressionFamily=cds0@expressionFamily, 
                                relative_expr=T,
                                disp_func=cds0@dispFitInfo[["blind"]]$disp_func,
                                verbose=F)
  res <- do.call(rbind.data.frame, lapply(diff_test_res, FUN = function(x){x$res}))
  
  res$qval <- 1
  res$qval[which(res$status == 'OK')] <- p.adjust(subset(res, status == 'OK')[, 'pval'], method="BH")
  
  res <- merge(res, fData(cds0), by="row.names")
  row.names(res) <- res[, 1] #remove the first column and set the row names to the first column
  res[, 1] <- NULL 
  
  res <- res[row.names(cds0), ]
  #cell_models0 <- fit_models(cds0,
  #                model_formula_str = full_formula,
  #               expression_family="negbinomial")
  #null_models0 <- fit_models(cds0,
  #                      model_formula_str = reduced_formula,
  #                     expression_family="negbinomial")
  #diff_test <- compare_models(cell_models0, null_models0)
  #res <- subset(diff_test, q_value < 0.05)
  
  
  #res <- diff_test
  print('diff test done')
  print(proc.time() - ptm)
  if(!is.null(subtype) & length(subtype) <= 2){
    # most direct calculation of log2FC, add pseudocounts of 1 for stability and slight shrunkage of values
    mat <- t(t(exprs(cds0))/sizeFactors(cds0))
    MeanExp <- tapply(colnames(mat), as.character(pData(cds0)[colnames(mat),]$cellType), function(x){rowMeans(mat[,x])})
    MeanExp <- log2((MeanExp[[subtype[2]]]+1)/(MeanExp[[subtype[1]]]+1))
    res$Log2FC_direct <- MeanExp[row.names(res)]
    
    if(foldChange == 'vst'){# using vst values to calculate log2FC, not recommended by authors
      mat <- vstExprs(cds0)
      MeanExp <- tapply(colnames(mat), as.character(pData(cds0)[colnames(mat),]$cellType), function(x){rowMeans(mat[,x])})
      MeanExp <- MeanExp[[subtype[2]]]- MeanExp[[subtype[1]]]
      res$Log2FC <- MeanExp[row.names(res)]
    }else if(foldChange == 'mle'){
      for(l in c("BiocGenerics", "VGAM", "Matrix")){library(l, character.only = T)}
      #ptm <- proc.time()
      #diff_test_res <- smartEsApply(cds0,1,diff_test_helper,
      #                      convert_to_dense=TRUE,
      #                      fullModelFormulaStr=full_formula,
      #                      reducedModelFormulaStr=reduced_formula, 
      #                      expressionFamily=cds0@expressionFamily, 
      #                      relative_expr=T,
      #                      disp_func=cds0@dispFitInfo[["blind"]]$disp_func,
      #                      verbose=F)
      #print('models fitted for effect extraction')
      #print(proc.time() - ptm)
      return(logfc_shrink(list(cds=cds0, res = res, models=diff_test_res), subtype = subtype, apeglm_method = apeglm_method ))
      
    }
    #res <- subset(res, abs(Log2FC) >= 1)
  }
  
  #res <- res[order(res$q_value),]
  #res <- res[order(res$qval),]
  return(res)
  #res <- subset(diff_test, qval < 0.05)
  #return(res)
}



logfc_shrink <- function(obj, subtype, apeglm_method = 'nbinomC'){
  res <- obj$res
  cds0 <- obj$cds
  diff_test_res <- obj$models
  oppose_signs = ifelse(sum(grepl(subtype[2], names(coef(diff_test_res[[1]]$full)))) > 0, F, T)
  ptm = proc.time()
  Log2FC_effect_disp <- do.call(rbind, lapply(row.names(res), FUN = function(x){
    disp = diff_test_res[[x]]$disp
    #print('1')
    coefs <- coef(diff_test_res[[x]]$full)
    coefs_SE <- sqrt(diag(vcov(diff_test_res[[x]]$full)))
    #print('2')
    if(!oppose_signs){
      effect <- coefs[grepl(subtype[2], names(coefs))]
      effect_SE <- coefs_SE[grepl(subtype[2], names(coefs_SE))]
    }else{
      effect <- -coefs[grepl(subtype[1], names(coefs))]
      effect_SE <- coefs_SE[grepl(subtype[1], names(coefs_SE))]
    }
    #print('3')
    return(c(effect, effect_SE, disp))
  }))
  print('Effects and Errors calculated')
  print(proc.time() - ptm)
  res$Log2FC_effect <- Log2FC_effect_disp[,1]
  res$Log2FC_effect_SE <- Log2FC_effect_disp[,2]
  res$disp <- Log2FC_effect_disp[,3]
  design <- diff_test_res[[1]]$full@x
  mle_prior=res[,c('Log2FC_effect', 'Log2FC_effect_SE')]
  if(oppose_signs){
    mle_prior[,1] <- -mle_prior[,1]
  }
  print('starting apeglm for logfc shrinkage')
  apeglm_mat <- exprs(cds0)[row.names(res),]
  ptm = proc.time()
  
  fit <- apeglm(Y= apeglm_mat, x=design, log.lik=logLikNB, param=as.matrix(res[,'disp', drop = F]),
                coef=2, mle=mle_prior, offset =  matrix(log(sizeFactors(cds0)), nrow=nrow(apeglm_mat), ncol=ncol(apeglm_mat), byrow=TRUE), method = apeglm_method)
  
  print('effects shrunk with apeglm')
  print(proc.time() - ptm)
  
  res$Log2FC_shrunk <- fit$map[row.names(res),2]*log2(exp(1))
  res$lfc_sval <- fit$svalue[row.names(res),]
  res$lfc_fsr <- fit$fsr[row.names(res),]
  
  if(oppose_signs){
    res$Log2FC_shrunk <- -res$Log2FC_shrunk
  }
  res$Log2FC_unshrunk <- Log2FC_effect_disp[,1]*log2(exp(1))
  
  res$emp_disp <- data.frame(cds0@dispFitInfo$blind$disp_table, row.names = 1)[row.names(res),'disp']
  res$est_disp <- cds0@dispFitInfo$blind$disp_func(data.frame(cds0@dispFitInfo$blind$disp_table, row.names = 1)[row.names(res),'mu'])

  return(res)
}

sparseApply <- function(Sp_X, MARGIN, FUN, convert_to_dense, ...){
  
  if (convert_to_dense){
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...) 
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...) 
      }, FUN, ...)
    }
  }else{
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...) 
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...) 
      }, FUN, ...)
    }
  }
  
  return(res)
  
}

smartEsApply <- function(X, MARGIN, FUN, convert_to_dense, ...) {
  isSparseMatrix <- function(x){
    class(x) %in% c("dgCMatrix", "dgTMatrix")
  }
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  Biobase::multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1
  
  if (sum(isSparseMatrix(exprs(X)))){
    res <- sparseApply(exprs(X), MARGIN, FUN, convert_to_dense, ...)
  }else{
    res <- apply(exprs(X), MARGIN, FUN, ...)
  }
  
  if (MARGIN == 1)
  {
    names(res) <- row.names(X)
  }else{
    names(res) <- colnames(X)
  }
  
  res
}


diff_test_helper <- function(x, 
                             fullModelFormulaStr, 
                             reducedModelFormulaStr, 
                             expressionFamily, 
                             relative_expr,
                             weights,
                             disp_func=NULL,
                             verbose=FALSE
){ 
  calculate_NB_dispersion_hint <- function(disp_func, f_expression, expr_selection_func=mean)
  {
    expr_hint <- expr_selection_func(f_expression)
    if (expr_hint > 0 && is.null(expr_hint) == FALSE) {
      disp_guess_fit <- disp_func(expr_hint)
      
      # For NB: Var(Y)=mu*(1+mu/k)
      f_expression_var <- var(f_expression)
      f_expression_mean <- mean(f_expression)
      
      disp_guess_meth_moments <- f_expression_var - f_expression_mean
      disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k
      
      #return (max(disp_guess_fit, disp_guess_meth_moments))
      return (disp_guess_fit)
    }
    return (NULL)
  }
  
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep="")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep="")
  
  x_orig <- x
  disp_guess <- 0
  
  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    if (relative_expr == TRUE)
    {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE){
      disp_guess <- calculate_NB_dispersion_hint(disp_func, round(x))
      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE  ) {
        # FIXME: In theory, we could lose some user-provided parameters here
        # e.g. if users supply zero=NULL or something. 
        if (expressionFamily@vfamily == "negbinomial")
          expressionFamily <- negbinomial(isize=1/disp_guess)
        else
          expressionFamily <- negbinomial.size(size=1/disp_guess)
      }
    }
  }else if (expressionFamily@vfamily %in% c("uninormal")){
    f_expression <- x
  }else if (expressionFamily@vfamily %in% c("binomialff")){
    f_expression <- x
    #f_expression[f_expression > 1] <- 1
  }else{
    f_expression <- log10(x)
  }
  
  test_res <- tryCatch({
    if (expressionFamily@vfamily %in% c("binomialff")){
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
      }
    }else{
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
      }
    }
    
    #list(full=full_model_fit, reduce = reduced_model_fit, disp = disp_guess, res = res)
    #print(coef(reduced_model_fit))
    res <- compareModels(list(full_model_fit), list(reduced_model_fit))
    list(full=full_model_fit, reduce = reduced_model_fit, disp = disp_guess, res = res)
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { 
    if(verbose)
      print (e);
    list(full=NULL, reduce = NULL, disp = disp_guess, res = data.frame(status = "FAIL", family=expressionFamily@vfamily, pval=1.0, qval=1.0))
    #data.frame(status = "FAIL", pval=1.0) 
  }
  )
  test_res
}



BootstrapNPLogFC <- function(mat1, mat2, perm = 100, exprs_fun = mean){
  logfc <- sapply(1:nrow(mat1), FUN =function(x){
    if(x>99 & x%%100 == 0){
      cat('calculated ', x, ' genes\n')
    }
    if(sum(mat1[x,]) == 0 | sum(mat2[x,]) == 0 | sum(mat1[x,] > 0) < 3 | sum(mat2[x,] > 0) < 3){
      return(NA)
    }else{
      mean(sapply(1:perm, FUN = function(m){ 
        m1 = exprs_fun(sample(mat1[x,], size = length(mat1[x,]), replace = T))
        m2 = exprs_fun(sample(mat2[x,], size = length(mat2[x,]), replace = T))
      return(log((m2+0.1)/(m1+0.1)))
      }))
    }
  })
  names(logfc) <- row.names(mat1)
  return(logfc)
}



