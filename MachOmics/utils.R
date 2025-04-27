# utils.R

library(dplyr)
library(ggplot2)
library(ggridges)
library(reshape2)
library(ranger)
library(pROC)
library(rstatix)
library(purrr)
library(patchwork)


# :: Check Data -----------------------------------------------------------

change_names = function(meta = meta,
                        type = "classification",
                        case = "case",
                        control = "control"){
  # change metadata column names
  colnames(meta) = c("sample", "value")
  if(type == "classification"){
    # replace names of classes
    meta$value = ifelse(meta$value == case, "case", "control")
    # and make them a factor
    meta$value = factor(meta$value, levels=c("case", "control"))
    
  }
  return(meta)
}

# :: Distribution Check ---------------------------------------------------


# plot distributions of 50 random genes

distribution_check = function(data = data,
                              scale = F,
                              transform = FALSE){
  if(transform == TRUE){
    data = log_transform(data)
  }
  if(ncol(data)>1000){
  set.seed(25)
  sample((data), size=50, replace=F)%>% as.matrix() %>%
    reshape2::melt() %>%
    ggplot(aes(x=value, y = 1))+
    ggridges::geom_density_ridges(aes(group = Var2), color="black", fill="white", alpha=0.2)+
    theme_minimal()+theme(legend.position="none")+
    labs(x="Value", y="Frequency", title="Sampled Features")
  }else{
    data %>% as.matrix() %>%
      reshape2::melt() %>%
      ggplot(aes(x=value, y = 1))+
      ggridges::geom_density_ridges(aes(group = Var2), color="black", fill="white", alpha=0.2)+
      theme_minimal()+theme(legend.position="none")+
      labs(x="Value", y="Frequency", title="All Features")
  }
}

# :: Log Transform --------------------------------------------------------

log_transform = function(data = data,
                         scale = F,
                         transform = TRUE){
  # replace na with 0
  data[is.na(data)] = 0
  
  if(transform == T){
  # calculate pseudocount from data
  shiny_pseudo = min(data[data!=0])/2
  # add pseudocount
  data = data+shiny_pseudo
  }
  
  # conditionally scale
  if(scale == T){
    data = sweep(data, 1, rowSums(data), FUN = "/")*100
  }
  
  # apply log2 transform
  if(transform == T){
    data = log(data, base=2)
  }else{
    data=data
  }
  # return data
  return(data)
}

# :: Check class imbalances -----------------------------------------------

check_imbalances = function(meta = meta,
                            type = "classification"){
  if(type == "classification"){
    # visualize class distribution
    imbalance_plot = ggplot(meta)+
      geom_bar(stat="count", position="stack", aes(fill=value, x=0), color="black")+
      coord_polar("y", start=0)+
      theme_void()+
      labs(fill="")
    return(imbalance_plot)
  }
  if(type == "regression"){
    imbalance_plot = ggplot(meta)+
      ggridges::geom_density_ridges(aes(x = value, y=0), color="black", fill="white", alpha=0.2)+
      theme_void()+
      labs(title="Outcome")
    return(imbalance_plot)
  }
}


# :: Build model ----------------------------------------------------------

build_model = function(data = data,
                       meta = meta,
                       task = "classification"){

  # merge data
  ml_data = data
  ml_data$sample = rownames(ml_data)
  ml_data = merge(ml_data, meta, by="sample")
  rownames(ml_data) = ml_data$sample
  ml_data$sample = NULL
  # done
  
  # check data type (if numeric, change type to classification)
  if(class(ml_data$value) == "numeric"){
    task = "regression"
  }
  if(class(ml_data$value) == "factor"){
    task = "classification"
  }
  
  if(task == "classification"){
    ml_loop_output = lapply(1:15, function(iter){
      lapply(c("real", "null"), function(type){
        # save true preds
        true.values = ml_data$value
        # if null, scramble data
        if(type == "null"){
          set.seed(iter)
          ml_data$value = sample(ml_data$value)
        }
        # build model
        set.seed(iter)
        ml_output = ranger::ranger(value ~.,
                                   ml_data,
                                   importance = "permutation",
                                   probability = TRUE)
        # extract predictions
        ml_predictions = data.frame(prob = ml_output$predictions[,1],
                                    true = true.values,
                                    sample = rownames(ml_data),
                                    iter = iter,
                                    type = type)
        # calculate sens, spec, thres, AUC
        ml_roc = pROC::roc(ml_predictions$true, ml_predictions$prob) %>% suppressMessages()
        ml_roc = data.frame(sens = ml_roc$sensitivities,
                            spec = ml_roc$specificities,
                            thresh = ml_roc$thresholds,
                            auc = rep(ml_roc$auc, times = length(ml_roc$specificities)))
        # add iter and datatype
        ml_roc$iter = iter
        ml_roc$type = type
        
        # extract importances
        ml_importances = data.frame(ml_output$variable.importance)
        ml_importances$feature = rownames(ml_importances)
        rownames(ml_importances) = NULL
        # add iter and datatype
        ml_importances$iter = iter
        ml_importances$type = type
        
        list(ml_predictions, ml_roc, ml_importances)
      })})
  }
  
  if(task == "regression"){
    ml_loop_output = lapply(1:15, function(iter){
      lapply(c("real", "null"), function(type){
        # save true preds
        true.values = ml_data$value
        # if null, scramble data
        if(type == "null"){
          set.seed(iter)
          # shuffle true preds
          ml_data$value = sample(ml_data$value)
        }
        # build model
        set.seed(iter)
        ml_output = ranger::ranger(value ~.,
                                   ml_data,
                                   importance = "permutation")
        # extract predictions
        ml_predictions = data.frame(pred = ml_output$predictions,
                                    true = true.values,
                                    sample = rownames(ml_data),
                                    iter = iter,
                                    type = type)
        # calculate R2 and Spearman Rho
        # R2
        ml_r2 = lm(true ~ pred, ml_predictions)
        ml_r2 = data.frame(
          metric = "R2",
          coef = summary(ml_r2)$r.squared,
          pval = summary(ml_r2)$coefficients[2,4])
        # Rho
        ml_cor = cor.test(ml_predictions$true,
                          ml_predictions$pred, method="spearman")
        ml_cor = data.frame(
          metric = "Rho",
          coef = ml_cor$estimate,
          pval = ml_cor$p.value)
        
        # merge
        ml_r2_cor = rbind(ml_r2,
                          ml_cor) %>% data.frame()
        
        # add iter and datatype
        ml_r2_cor$iter = iter
        ml_r2_cor$type = type
        
        # extract importances
        ml_importances = data.frame(ml_output$variable.importance)
        ml_importances$feature = rownames(ml_importances)
        rownames(ml_importances) = NULL
        # add iter and datatype
        ml_importances$iter = iter
        ml_importances$type = type
        
        list(ml_predictions, ml_r2_cor, ml_importances)
        
      })})
    
    # return a list of performance metrics (per iter and type) and importances
    
  }
  return(ml_loop_output)
  
}


# :: Visualize Performance ------------------------------------------------

visualize_performance = function(build_output = build_output,
                                 task = "classification"){
  # extract performance data
  predictions_data = do.call(rbind, lapply(1:15, function(x){
    do.call(rbind, lapply(1:2, function(y){
      purrr::pluck(build_output[[x]][[y]][[1]])
    }))}))
  performance_data = do.call(rbind, lapply(1:15, function(x){
    do.call(rbind, lapply(1:2, function(y){
      purrr::pluck(build_output[[x]][[y]][[2]])
    }))}))
  if(task == "classification"){
    # plot confusion
    real_conf = subset(predictions_data, type == "real")
    real_conf$pred = ifelse(real_conf$prob > 0.5, "case", "control")
    real_conf = data.frame(real_conf[,c("iter", "true", "pred")] %>% table()) %>%
      # calculate median
      group_by(true, pred) %>%
      mutate(med.freq = median(Freq))%>%
      # calculate sd
      mutate(sd.freq = round(sd(Freq), digits=2)) %>%
      # clean df
      dplyr::select(-iter, -Freq) %>% distinct() %>% data.frame()
    plot_confusion_matrix = ggplot(real_conf,
                                   aes(x=true, y=pred))+
      geom_tile(aes(fill=ifelse(true == pred, 1, "grey")), color="black")+
      scale_fill_manual(values=c(1, "grey"))+
      geom_text(aes(label=paste0(med.freq,"\n", "(", sd.freq, ")", sep="")), color="white")+
      guides(fill="none", color="none")+
      theme_minimal()+
      labs(x="True", y="Prediction")
    
    # plot AUC
    performance_data_auc = performance_data %>% dplyr::select(-sens, -spec) %>% distinct()
    performance_data_auc_median = performance_data_auc %>% group_by(type) %>% mutate(med.auc = median(auc))%>%
      dplyr::select(med.auc, type) %>% distinct() %>% data.frame()
    # run wilcox
    stat.test.auc <- performance_data_auc %>%
      rstatix::wilcox_test(auc ~ type) %>%
      mutate(label = ifelse(p < 0.05, "*", "ns"))
    
    plot_auc = ggplot(performance_data_auc,
                      aes(x=type, y=auc))+
      geom_violin(data=subset(performance_data_auc, type == "null"),
                  fill="grey", alpha=0.4, color=NA)+
      geom_violin(data=subset(performance_data_auc, type == "real"),
                  fill="red", color="black", alpha=0.4)+
      geom_point(data=subset(performance_data_auc_median, type == "null"), aes(x=type, y=med.auc),
                 shape=21, size=3, fill="grey", alpha=0.8)+
      geom_point(data=subset(performance_data_auc_median, type == "real"), aes(x=type, y=med.auc),
                 shape = 21, size=3, fill="white", alpha=1)+
      geom_segment(x=1, xend=2,
                   y=max(performance_data_auc$auc)+0.02,
                   yend=max(performance_data_auc$auc)+0.02)+
      geom_text(x=1.5, y=max(performance_data_auc$auc)+0.03,
                label = stat.test.auc$label,
                size=ifelse(stat.test.auc$label == "*", 5, 3))+
      expand_limits(y = max(performance_data_auc$auc) + 0.05)+
      theme_minimal()+
      labs(x="", y="AUC")
    plot_auc
    
    # plot ROC
    performance_data_mean = performance_data %>%
      # remove Inf threshold
      mutate(thresh = ifelse(thresh > 1, 1, 
                             ifelse(thresh < 0, 0, thresh))) %>%
      # round threshold to 0.01
      mutate(thresh.round = round(thresh, digits=2)) %>%
      # calculate the mean of 15 iterations
      group_by(type, thresh.round) %>%
      mutate(mean.sens = mean(sens),
             mean.spec = mean(spec)) %>%
      # second, calculate 95% confidence (sens) interval using 15 iterations
      mutate(lower.sens = quantile(sens, 0.05),
             upper.sens = quantile(sens, 0.95)) %>%
      # last, calculate mean AUC
      group_by(type)%>%
      mutate(mean.auc = mean(auc)) %>%
      # clean df
      dplyr::select(type, thresh.round, mean.sens, mean.spec, lower.sens, upper.sens, mean.auc) %>% distinct() %>% data.frame()
    
    # replace type with AUC, too
    performance_data_mean$type.auc = ifelse(performance_data_mean$type == "real",
                                            paste0("Real (AUC: ", round(unique(subset(performance_data_mean, type == "real")$mean.auc), digits=2),")", sep=""),
                                            paste0("Null (AUC: ", round(unique(subset(performance_data_mean, type == "null")$mean.auc), digits=2),")", sep=""))
    
    performance_data_mean$type.auc = factor(performance_data_mean$type.auc, levels=(unique(performance_data_mean$type.auc)))
    
    # force ordering
    performance_data_mean = performance_data_mean %>%
      group_by(type) %>%
      mutate(mean.sens = -sort(-mean.sens),
             mean.spec = sort(mean.spec),
             lower.sens = -sort(-lower.sens),
             upper.sens = -sort(-upper.sens))
    
    # draw plot
    plot_roc = ggplot(subset(performance_data_mean, type %in% c("real", "null")),
                      aes(x=1-mean.spec, y=mean.sens))+
      # ROC curve
      geom_line(aes(group=type.auc, color=type.auc), size=1)+
      # 95% confidence interval
      geom_ribbon(aes(x = 1-(mean.spec), ymin = lower.sens, ymax = upper.sens,
                      fill=type.auc), alpha = 0.2)+
      # Add diagonal reference line (random classifier)
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      # change colors
      scale_color_manual(values = c("red", "grey"))+
      scale_fill_manual(values = c("red", "grey"))+
      # Labels and theme
      labs(x = "1 - Specificity", y = "Sensitivity", fill="", color="") +
      theme_minimal()+theme(legend.position = c(0.75,0.25),
                            legend.text=element_text(size=11))
    return(plot_confusion_matrix | plot_auc | plot_roc)
    
    
  }
  if(task == "regression"){
    # plot correlations
    correlation_data = do.call(rbind, lapply(1:15, function(x){
      do.call(rbind, lapply(1:2, function(y){
        purrr::pluck(build_output[[x]][[y]][[1]])
      }))}))
    performance_data = do.call(rbind, lapply(1:15, function(x){
      do.call(rbind, lapply(1:2, function(y){
        purrr::pluck(build_output[[x]][[y]][[2]])
      }))}))
    
    # plot Spearman correlations
    performance_data_cor = performance_data %>% subset(metric == "Rho")
    performance_data_cor_median = performance_data_cor %>% group_by(type) %>% mutate(med.cor = median(coef))%>%
      dplyr::select(med.cor, type) %>% distinct() %>% data.frame()
    # run wilcox
    stat.test.cor <- performance_data_cor %>%
      rstatix::wilcox_test(coef ~ type) %>%
      mutate(label = ifelse(p < 0.05, "*", "ns"))
    
    plot_spear = ggplot(performance_data_cor,
                        aes(x=type, y=coef))+
      geom_violin(data=subset(performance_data_cor, type == "null"),
                  fill="grey", alpha=0.4, color=NA)+
      geom_violin(data=subset(performance_data_cor, type == "real"),
                  fill="red", color="black", alpha=0.4)+
      geom_point(data=subset(performance_data_cor_median, type == "null"), aes(x=type, y=med.cor),
                 shape=21, size=3, fill="grey", alpha=0.8)+
      geom_point(data=subset(performance_data_cor_median, type == "real"), aes(x=type, y=med.cor),
                 shape = 21, size=3, fill="white", alpha=1)+
      geom_segment(x=1, xend=2,
                   y=max(performance_data_cor$coef)+0.02,
                   yend=max(performance_data_cor$coef)+0.02)+
      geom_text(x=1.5, y=max(performance_data_cor$coef)+0.03,
                label = stat.test.cor$label,
                size=ifelse(stat.test.cor$label == "*", 5, 3))+
      expand_limits(y = max(performance_data_cor$coef) + 0.05)+
      theme_minimal()+
      labs(x="", y="Spearman ρ")
    plot_spear
    
    # plot predictions
    correlation_data_median = correlation_data %>% 
      group_by(type, sample) %>%
      mutate(med.pred = median(pred)) %>%
      dplyr::select(type, sample, med.pred, true) %>% dplyr::distinct() %>% data.frame()
    # add median correlation to title
    correlation_data_median$type.cor = ifelse(correlation_data_median$type == "real",
                                              paste0("Real (Rho: ", round(unique(subset(performance_data_cor_median, type == "real")$med.cor), digits=2),")", sep=""),
                                              paste0("Null (Rho: ", round(unique(subset(performance_data_cor_median, type == "null")$med.cor), digits=2),")", sep=""))
    correlation_data_median$type.cor = factor(correlation_data_median$type.cor, levels=(unique(correlation_data_median$type.cor)))
    
    # plot
    plot_cor = ggplot(correlation_data_median,
                      aes(x=true, y=med.pred))+
      # plot null first
      geom_smooth(data=subset(correlation_data_median, type == "null"), aes(color=type.cor), method="lm")+
      geom_point(data=subset(correlation_data_median, type == "null"), aes(fill=type.cor), shape=21, color="black")+
      # plot real second
      geom_smooth(data=subset(correlation_data_median, type == "real"), aes(color=type.cor), method="lm")+
      geom_point(data=subset(correlation_data_median, type == "real"), aes(fill=type.cor), shape=21, color="black")+
      scale_fill_manual(values=c("grey", "red"))+
      scale_color_manual(values=c("grey", "red"))+
      theme_minimal()+theme(legend.position = c(0.75,0.15),
                            legend.text=element_text(size=11))+
      labs(x="True", y="Prediction", color="", fill="")
    
    return(plot_spear | plot_cor)
    
  }
}

# :: Feature Importances ------------------------------------------------

feature_importances = function(build_output = build_output,
                               data = data,
                               meta = meta){

  # extract feature importance data
  importance_data = do.call(rbind, lapply(1:15, function(x){
    do.call(rbind, lapply(1:2, function(y){
      purrr::pluck(build_output[[x]][[y]][[3]])
    }))}))
  
  # clean importance_data
  colnames(importance_data)[1] = "importance"
  #importance_data = subset(importance_data, importance > 0)
  
  # prepare data for lm
  # merge data
  ml_data = data
  ml_data$sample = rownames(ml_data)
  ml_data = merge(ml_data, meta, by="sample")
  rownames(ml_data) = ml_data$sample
  ml_data$sample = NULL
  # done
  
  # loop through features
  importance_coef = do.call(rbind, lapply(1:(ncol(ml_data)-1), function(x){
    # subset data
    data.subset = ml_data[,c(x, ncol(ml_data))]
    feature.name = colnames(data.subset)[1]
    colnames(data.subset)[1] = "feature"
    # run lm
    lm.output = lm(feature ~ value, data.subset) %>% summary() %>% coef() %>% data.frame()
    # extract coef
    lm.output = lm.output[2,1]
    # save as dataframe
    lm.output = data.frame(feature = feature.name,
                           coef = lm.output)
    lm.output
  }))
  
  # assess significance of features (real vs null, across 15 iterations; apply BH-FDR)
  feature_wilcox = do.call(rbind, lapply(unique(importance_data$feature), function(x){
    # subset to a single feature
    data.subset = subset(importance_data, feature == x)
    # calculate wilcox difference
    wilcox.output = wilcox.test(subset(data.subset, type == "real")$importance,
                                subset(data.subset, type == "null")$importance)
    # save as dataframe
    wilcox.output = data.frame(
      feature = x,
      estimate = wilcox.output$statistic,
      pval = wilcox.output$p.value
    )
    wilcox.output
  }))
  # apply fdr
  feature_wilcox$padj = p.adjust(feature_wilcox$pval, method="BH")
  # sig
  feature_wilcox$sig = ifelse(feature_wilcox$padj >= 0.05 | is.na(feature_wilcox$padj), "", "*")
  # merge with coef
  feature_wilcox = merge(feature_wilcox, importance_coef, by="feature")
  # calculate median value
  feature_median = importance_data %>%
    group_by(type, feature) %>%
    mutate(med.imp = median(importance)) %>% 
    dplyr::select(type, feature, med.imp) %>% distinct() %>% data.frame()
  # merge
  feature_importance_data = merge(feature_wilcox, feature_median, by="feature")
  
  # plot feature association with outcome
  
  # plot
  plot_feature_interpretation = ggplot(subset(feature_importance_data, type=="real") %>% mutate(feature.label = ifelse(sig == "*" & med.imp > 0, as.character(feature), NA)),
                                       aes(x=coef, y=med.imp, fill=coef, size=med.imp) )+
    # 
    geom_point(shape=21, alpha=1)+
    ggrepel::geom_text_repel(aes(label=feature.label), size=3)+
    theme_minimal()+theme(legend.position="none")+
    scale_fill_gradient2(low = scales::muted("red"),
                         high = scales::muted("blue"))+
    labs(x="Coefficient", y="Median Feature Importance")
  
  if(sum(feature_importance_data$sig=="*")>0){
    # subset to sig
    feature_importance_data = subset(feature_importance_data, sig == "*" & med.imp > 0)
    importance_data_sig = subset(importance_data, feature %in% feature_importance_data$feature)
    feature_median = subset(feature_median, feature %in% feature_importance_data$feature)
    
    # order based on median
    feature_importance_data$feature = factor(feature_importance_data$feature, levels=unique(arrange(subset(feature_importance_data, type=="real" & med.imp > 0), med.imp)$feature))
    importance_data_sig$feature = factor(importance_data_sig$feature, levels=unique(arrange(subset(feature_importance_data, type=="real"& med.imp > 0),med.imp)$feature))
    feature_median$feature = factor(feature_median$feature, levels=unique(arrange(subset(feature_importance_data, type=="real"& med.imp > 0),med.imp)$feature))
    
    # plot
    plot_feature_importance = ggplot()+
      # plot null results
      geom_violin(data=subset(importance_data_sig, type =="null" & !is.na(feature)), 
                  aes(x=importance, y=feature), fill="grey", alpha=0.4)+
      geom_point(data=subset(feature_median, type =="null" & !is.na(feature) & feature %in% unique(arrange(feature_importance_data,med.imp)$feature)), 
                 aes(x=med.imp, y=feature), shape=21, fill="grey", alpha=0.4)+
      # plot real results
      geom_violin(data=subset(importance_data_sig, type =="real"& !is.na(feature)), 
                  aes(x=importance, y=feature), fill="red", alpha=0.4)+
      geom_point(data=subset(feature_median, type =="real" & !is.na(feature)& feature %in% unique(arrange(feature_importance_data,med.imp)$feature)), 
                 aes(x=med.imp, y=feature), shape=21, fill="red", alpha=1)+
      theme_minimal()+
      labs(x="Feature Importance", y="")
    
    # plot both
    plot_feature_interpretation | plot_feature_importance
    
  } else {
    
    plot_feature_interpretation
    
    print("No significant features identified")
  }
  
  
}


