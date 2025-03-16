### 2025_03_04  Piloting MLP

BiocManager::install("MLP")
library("MLP")
MLP(sample.16s.data,
         "rdp_train_set_16", 
         "galaxy", 
         "load")
# training datasets are not 16S
# does not work

# :: load model from their 16S data ------------------------------------------

# load XGB model
model.16s.xgb = readRDS("~/Downloads/model.16S_rRNA.rds")
# extract data
model.16s.xgb.training.data = model.16s.xgb$trainingData

# check features
varImp(model.16s.xgb)
 
nrow(model.16s.xgb.training.data)

cor.test(model.16s.xgb$pred$pred,
         model.16s.xgb$pred$obs)
# Pearson cor = 0.104!

github.model.16s.selected.predictions = model.16s.xgb$pred %>%
  subset(eta == model.16s.xgb$bestTune$eta &
           nrounds == model.16s.xgb$bestTune$nrounds &
           max_depth == model.16s.xgb$bestTune$max_depth &
           gamma == model.16s.xgb$bestTune$gamma &
           colsample_bytree == model.16s.xgb$bestTune$colsample_bytree &
           min_child_weight == model.16s.xgb$bestTune$min_child_weight &
           subsample == model.16s.xgb$bestTune$subsample) %>%
  group_by(rowIndex) %>%
  mutate(ave.pred = mean(pred)) %>%
  dplyr::select(ave.pred, rowIndex, obs) %>% distinct() %>% data.frame()

cor.test(github.model.16s.selected.predictions$ave.pred,
         github.model.16s.selected.predictions$obs)
# Cor = 0.789

# ``` plot: training ------------------------------------------------------
original.xgb.model = data.frame(pred = github.model.16s.selected.predictions$ave.pred,
           true = github.model.16s.selected.predictions$obs) %>%
  ggplot(aes(x=true, y=pred))+
  geom_smooth(method="lm", se=F, color="black")+
  geom_point(alpha=0.2)+
  theme_bw()+
  labs(x="Measured", y="Predicted")
# Cor = 0.789

# :: rebuild model from their 16S data ------------------------------------

colnames(model.16s.xgb.training.data)[ncol(model.16s.xgb.training.data)] <- "outcome"

# NOTE: this tune grid (and repeats=1) yields a recapitulation of their model; with DISMAL performance

tune_grid = expand.grid(
  nrounds=c(100, 500, 1000),
  max_depth = c(3:5),
  eta = c(0.01, 0.05, 0.1),
  gamma = c(0.01),
  colsample_bytree = c(0.75),
  subsample = 1,
  min_child_weight = 1
)

# train model with caret
t1 = Sys.time()
model.16s.selected = caret::train(
  outcome ~ .,
  data = model.16s.xgb.training.data,
  trControl = caret::trainControl(method="repeatedcv", 
                                  number=10, 
                                  repeats = 1, 
                                  savePredictions = TRUE),
  #tuneGrid = tune_grid,
  method = "xgbTree",
  preProc = c("center", "scale")
)
t2 = Sys.time()
t2 - t1

model.16s.selected.predictions = model.16s.selected$pred %>%
  subset(eta == model.16s.selected$bestTune$eta &
           nrounds == model.16s.selected$bestTune$nrounds &
           max_depth == model.16s.selected$bestTune$max_depth &
           gamma == model.16s.selected$bestTune$gamma &
           colsample_bytree == model.16s.selected$bestTune$colsample_bytree &
           min_child_weight == model.16s.selected$bestTune$min_child_weight &
           subsample == model.16s.selected$bestTune$subsample) %>%
  group_by(rowIndex) %>%
  mutate(ave.pred = mean(pred)) %>%
  dplyr::select(ave.pred, rowIndex, obs) %>% distinct() %>% data.frame()

cor.test(model.16s.selected.predictions$ave.pred,
         model.16s.selected.predictions$obs)
# Cor = 0.79

varImp(model.16s.xgb)
varImp(model.16s.selected)



# ``` plot: training ------------------------------------------------------
retrained.xgb.model = data.frame(pred = model.16s.selected.predictions$ave.pred,
           true = model.16s.selected.predictions$obs) %>%
  ggplot(aes(x=true, y=pred))+
  geom_smooth(method="lm", se=F, color="black")+
  geom_point(alpha=0.2)+
  theme_bw()+
  labs(x="Measured", y="Predicted")
# Cor = 0.79


original.xgb.model | retrained.xgb.model

# :: vandeputte test ------------------------------------------------------


# load test data
sample.16s.data = read.csv("~/Downloads/Vandeputte_2017_16S.tsv", sep="\t") # supp table 14
dim(sample.16s.data) # 95 samples
rownames(sample.16s.data) <- sample.16s.data$X
sample.16s.data$X <- NULL
# convert to %
sample.16s.data = sweep(sample.16s.data, 1, rowSums(sample.16s.data), FUN = "/")

# make a function to add missing variables

# add missing variables to testing data
missing.taxa = colnames(model.16s.selected$trainingData)[!colnames(model.16s.selected$trainingData) %in%
                                                 colnames(sample.16s.data)]
sample.16s.data = cbind(sample.16s.data, 
                        do.call(cbind, lapply(missing.taxa, function(x){
  new.data = data.frame(x = rep(0, times = nrow(sample.16s.data)))
  colnames(new.data) = x
  new.data
}))) %>% data.frame()

sample.16s.data$X = NULL
pseudo.test.van = min(sample.16s.data[sample.16s.data!=0])/2
sample.16s.data = log(sample.16s.data+pseudo.test.van, base=10)

model.output = predict(model.16s.selected, sample.16s.data)

hist(model.output)
# check real values

test.16s.real = read.csv("~/Downloads/41586_2017_BFnature24460_MOESM14_ESM.csv")

plot(model.output,
         log(test.16s.real$Average.cell.count..per.gram.of.frozen.feces.))
cor.test(model.output,
     log(test.16s.real$Average.cell.count..per.gram.of.frozen.feces., base=10),
     method="spearman")
cor.test(model.output,
         log(test.16s.real$Average.cell.count..per.gram.of.frozen.feces., base=10),
         method="pearson")
# Hmm. Pearson R = 0.610; Spearman = 0.548


# :: loop through train + test --------------------------------------------

sample.16s.data = read.csv("~/Downloads/Vandeputte_2017_16S.tsv", sep="\t")
test.16s.real = read.csv("~/Downloads/41586_2017_BFnature24460_MOESM14_ESM.csv")

# note: rpart and nnet do not work; all predict the same value

ml.models = c("ranger", "gbm", "xgbTree", "lm", "glmnet", "svmLinear", "svmRadial")

t0 <- Sys.time()
test.16s.iter = do.call(rbind, lapply(1:15, function(iter) {
  do.call(rbind, lapply(ml.models, function(model){
    # iter = 1
    # model = "ranger"
    print(paste0(iter, model))
    # train model with caret
    t1 = Sys.time()
    set.seed(iter)
    model.16s.ml = caret::train(
      outcome ~ .,
      data = model.16s.xgb.training.data,
      trControl = caret::trainControl(method="cv", 
                                      number=5, 
                                      savePredictions = TRUE),
      method = model,
      verbose=F,
      preProc = c("center", "scale")) %>% suppressWarnings()
    t2 = Sys.time()
    t2 - t1
    
    # fix test data
    sample.16s.data.processed = sample.16s.data
    rownames(sample.16s.data.processed) <- sample.16s.data.processed$X
    sample.16s.data.processed$X <- NULL
    # convert to %
    sample.16s.data.processed = sweep(sample.16s.data.processed, 1, rowSums(sample.16s.data.processed), FUN = "/")
    
    # add missing variables to testing data
    missing.taxa = colnames(model.16s.ml$trainingData)[!colnames(model.16s.ml$trainingData) %in%
                                                                    colnames(sample.16s.data.processed)]
    if(length(missing.taxa)>0){
    sample.16s.data.processed = cbind(sample.16s.data.processed, 
                            do.call(cbind, lapply(missing.taxa, function(x){
                              new.data = data.frame(x = rep(0, times = nrow(sample.16s.data.processed)))
                              colnames(new.data) = x
                              new.data
                            }))) %>% data.frame()
    }
    pseudo.test.van.loop = min(sample.16s.data.processed[sample.16s.data.processed!=0])/2
    sample.16s.data.processed = log(sample.16s.data.processed+pseudo.test.van.loop, base=10)
    
    # predict test data
    
    model.output = predict(model.16s.ml, sample.16s.data.processed)

    # save
    data.frame(true = log(test.16s.real$Average.cell.count..per.gram.of.frozen.feces., base=10),
               pred = model.output,
               iter = iter,
               model = model,
               time = as.numeric(t2 - t1))
}))}))
t3 <- Sys.time()
t3 - t0
#saveRDS(test.16s.iter, "~/Documents/PhD/For others/2025_03_05_mlp_16s_test.Rds")
test.16s.iter.df <- test.16s.iter


# plot raw vals
ggplot(test.16s.iter.df,
       aes(x=true, y=pred))+
  geom_smooth(aes(group = iter), method="lm")+
  geom_point()+
  facet_wrap( ~ model, scales="free")

# evaluate using Spearman

test.16s.iter.eval = do.call(rbind, lapply(1:15, function(seed) {
  do.call(rbind, lapply(ml.models, function(ml){
    print(paste0(ml, seed))
    data.subset = subset(test.16s.iter.df, iter == seed & model == ml) 
    #
    spear = cor.test(data.subset$true,
             data.subset$pred, method="spearman")$estimate
    pear = cor.test(data.subset$true,
                     data.subset$pred, method="spearman")$estimate
    if(sd(data.subset$pred)==0){
      r2 = NA
    }else{
    r2 = summary(lm(scale(pred) ~ scale(true), data.subset))$r.squared
    }
    # save
    data.frame(model = ml,
               iter = seed,
               spearman = spear,
               pearson = pear,
               r2 = r2)
  }))}))


# ```plot: model evaluation -----------------------------------------------

# plot time
ggplot(test.16s.iter.df[,c("iter", "model", "time")] %>% distinct())+
  geom_boxplot(aes(x = reorder(model, time), y=time))+
  geom_point(aes(x = reorder(model, time), y=time))+
  theme_bw()+
  labs(x="Model", y="Time elapsed (sec)")
# RF is the fastest
# glmnet is not bad; much better than svmLinear

# plot R2
ggplot(test.16s.iter.eval)+
  geom_boxplot(aes(x = model, y=r2))+
  geom_point(aes(x = model, y=r2))+
  theme_bw()+
  labs(x="Model", y="R-Squared")

# plot pearson
test.pearson.plot = ggplot(test.16s.iter.eval)+
  geom_boxplot(aes(x = model, y=pearson))+
  geom_point(aes(x = model, y=pearson), size=2.5, color="black")+
  geom_point(aes(x = model, y=pearson), size=1.5, color="white")+
  geom_point(aes(x = model, y=pearson, color=model), size=1.5, alpha=0.6)+
  theme_bw()+theme(legend.position="none")+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x="Model", y="Pearson Correlation")
# why is there no variation with glmnet, etc?
# it is deterministic

test.time.plot = ggplot(subset(test.16s.iter.df, model %in% ml.models)[,c("iter", "model", "time")] %>% distinct())+
  geom_boxplot(aes(x = model, y=time))+
  geom_point(aes(x = model, y=time), size=2.5, color="black")+
  geom_point(aes(x = model, y=time), size=1.5, color="white")+
  geom_point(aes(x = model, y=time, color=model), size=1.5, alpha=0.6)+
  theme_bw()+theme(legend.position="none")+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x="Model", y="Seconds Elapsed")

test.pearson.plot/test.time.plot

# :: REBUILT MODEL ----------------------------------------------


# train model with caret
t1 = Sys.time()
model.16s.rf = caret::train(
  outcome ~ .,
  data = model.16s.xgb.training.data,
  trControl = caret::trainControl(method="repeatedcv", 
                                  number=5, # no increase in internal validation performance if set to 5
                                  repeats = 1, # no increase in internal validation performance if set to 5 
                                  savePredictions = TRUE),
  method = "ranger",
  preProc = c("center", "scale")
)
t2 = Sys.time()
t2 - t1

# don't use caret
model.16s.rf = ranger::ranger(outcome ~., model.16s.xgb.training.data)
cor.test(model.16s.rf$predictions,
         model.16s.xgb.training.data$outcome)
# Cor = 0.784


cor.test(model.16s.rf$finalModel$predictions,
           model.16s.xgb.training.data$outcome)
# Cor = 0.788

retrained.rf.model = data.frame(pred = model.16s.rf$finalModel$predictions,
                                true = model.16s.xgb.training.data$outcome) %>%
  ggplot(aes(x=true, y=pred))+
  geom_smooth(method="lm", se=F, color="black")+
  geom_point(alpha=0.2)+
  theme_bw()+
  labs(x="Measured", y="Predicted")


# load test data
sample.16s.data = read.csv("~/Downloads/Vandeputte_2017_16S.tsv", sep="\t") # supp table 14
dim(sample.16s.data) # 95 samples
rownames(sample.16s.data) <- sample.16s.data$X
sample.16s.data$X <- NULL
# convert to %
sample.16s.data = sweep(sample.16s.data, 1, rowSums(sample.16s.data), FUN = "/")

# add missing variables to testing data
missing.taxa = colnames(model.16s.rf$trainingData)[!colnames(model.16s.rf$trainingData) %in%
                                                                colnames(sample.16s.data)]
sample.16s.data = cbind(sample.16s.data, 
                        do.call(cbind, lapply(missing.taxa, function(x){
                          new.data = data.frame(x = rep(0, times = nrow(sample.16s.data)))
                          colnames(new.data) = x
                          new.data
                        }))) %>% data.frame()
pseud.test.test = min(sample.16s.data[sample.16s.data !=0])/2
sample.16s.data = log(sample.16s.data[,colnames(sample.16s.data) != "outcome"]+pseud.test.test, base=10)


model.output = predict(model.16s.rf, sample.16s.data)
saveRDS(model.16s.rf, "~/Documents/PhD/For others/2025_03_06_mlp_rf.Rds")

test.16s.real = read.csv("~/Downloads/41586_2017_BFnature24460_MOESM14_ESM.csv")

cor.test(model.output,
         log(test.16s.real$Average.cell.count..per.gram.of.frozen.feces.),
         method="spearman")
# Hmm. Pearson R = 0.602

retrained.rf.model.test = data.frame(pred = model.output,
                                    true = log(test.16s.real$`Average.cell.count..per.gram.of.frozen.feces.`, base=10)) %>%
  ggplot(aes(x=true, y=pred))+
  geom_smooth(method="lm", se=F, color="black")+
  geom_point(alpha=0.2)+
  geom_text(x=-Inf, y=Inf, label="R = 0.788", hjust=0, vjust = 1, nudge_x=10, nudge_y=-2)+
  theme_bw()+
  labs(x="Measured", y="Predicted", title="External Validation")

retrained.rf.model.test

retrained.rf.model|retrained.rf.model.test
# Cor = 0.788 vs Cor = 0.602


# :: Test on RapidAIM data ------------------------------------------------

# load data
peter.rds = readRDS("~/Documents/PhD/git_pd_archfolder/2025_02_07_asv_mend_data_filtered_rare.Rds")
# make tax table
peter.rds.tax = data.frame(phyloseq::tax_table(peter.rds))
peter.rds.tax$seq = rownames(peter.rds.tax)
# apply RDP
combined.tax.table.rdp <- dada2::assignTaxonomy(unique(as.character(peter.rds.tax$seq)), 
                                                  "~/Downloads/rdp_train_set_16.fa.gz")
# format RDP to be compatible with MLP taxa naming
combined.tax.table.rdp.df = combined.tax.table.rdp %>% data.frame()
# make LCA
combined.tax.table.rdp.df$LCA = ifelse(!is.na(combined.tax.table.rdp.df$Genus), paste("g_", combined.tax.table.rdp.df$Genus, sep=""),
                                         ifelse(!is.na(combined.tax.table.rdp.df$Family), paste("uc_f_", combined.tax.table.rdp.df$Family, sep=""),
                                                ifelse(!is.na(combined.tax.table.rdp.df$Order), paste("uc_o_", combined.tax.table.rdp.df$Order, sep=""),
                                                       ifelse(!is.na(combined.tax.table.rdp.df$Class), paste("uc_c_", combined.tax.table.rdp.df$Class, sep=""),
                                                              ifelse(!is.na(combined.tax.table.rdp.df$Phylum), paste("uc_p_", combined.tax.table.rdp.df$Phylum, sep=""),
                                                                     combined.tax.table.rdp.df$Kingdom)))))
combined.tax.table.rdp.df$seq = rownames(combined.tax.table.rdp.df)
# make OTU table
peter.rds.df = speedyseq::psmelt(peter.rds)
peter.rds.df$LCA = combined.tax.table.rdp.df[match(peter.rds.df$OTU, combined.tax.table.rdp.df$seq),]$LCA

peter.rds.mat = peter.rds.df %>%
  group_by(HM, RS_Name, Replicate, LCA) %>%
  mutate(sum.abun = sum(Abundance)/50000) %>%
  mutate(sample = paste(HM, RS_Name, Replicate, sep="_")) %>%
  dplyr::select(sample, LCA, sum.abun) %>% distinct() %>%
  reshape2::acast(sample ~ LCA, value.var="sum.abun") %>% data.frame()
# done


# add missing variables to testing data
missing.taxa.peter = colnames(model.16s.glmnet$trainingData)[!colnames(model.16s.glmnet$trainingData) %in%
                                                                colnames(peter.rds.mat)] %>% sort()

peter.rds.mat.plus = cbind(peter.rds.mat, 
                        do.call(cbind, lapply(missing.taxa.peter, function(x){
                          new.data = data.frame(x = rep(0, times = nrow(peter.rds.mat)))
                          colnames(new.data) = x
                          new.data
                        }))) %>% data.frame()
pseud.test.peter = min(peter.rds.mat.plus[peter.rds.mat.plus !=0])/2
peter.rds.mat.plus = log(peter.rds.mat.plus[,colnames(peter.rds.mat.plus) != "outcome"]+pseud.test.peter, base=10)

## APPLY MODEL
model.output.peter = predict(model.16s.rf, peter.rds.mat.plus)

hist(model.output.peter)

# apply algorithm
peter.rds.predicted <- data.frame(sample = rownames(peter.rds.mat.plus),
                                  load = model.output.peter)
# merge with biomass
peter.rapidaim.predicted = phyloseq::sample_data(peter.rds) %>% data.frame()
peter.rapidaim.predicted$sample = paste(peter.rapidaim.predicted$HM,
                                        peter.rapidaim.predicted$RS_Name,
                                        peter.rapidaim.predicted$Replicate, sep="_")

peter.rapidaim.predicted <- merge(peter.rapidaim.predicted,
                                  peter.rds.predicted, by="sample")
peter.rapidaim.predicted$DNA = as.numeric(peter.rapidaim.predicted$Qubit)
peter.rapidaim.predicted = subset(peter.rapidaim.predicted, !is.na(DNA))

# subset to RapidAIM samples
peter.rapidaim.predicted = subset(peter.rapidaim.predicted, RS_Name %in% rs.names.pbs)
peter.rapidaim.predicted$RS_Name = factor(peter.rapidaim.predicted$RS_Name, levels=rs.names.pbs)
hist(peter.rapidaim.predicted$load)

peter.rapidaim.predicted.plot = ggplot(peter.rapidaim.predicted,
       aes(x=load, y=log(DNA,base=10)))+
  geom_smooth(method="lm", se=F, color="black")+
  geom_point(aes(color = RS_Name), alpha=0.3)+
  scale_color_manual(values=labelcolors$cols[c(10,1:9)])+
  theme_bw()+theme(legend.position="none")+
  labs(x="Predicted Microbial Load", y="Log10(DNA ng/ul)")
cor.test(peter.rapidaim.predicted$load,
         log(peter.rapidaim.predicted$DNA, base=10))
# Pearson Cor = 0.257, p < 0.001
cor.test(peter.rapidaim.predicted$load,
         log(peter.rapidaim.predicted$DNA, base=10), method="spearman")
# Spearman Rho = 0.327, p < 0.001

lmerTest::lmer(load ~ log(DNA+0.05) + (1|HM) + (1|RS_Name), peter.rapidaim.predicted) %>% summary()
lmerTest::lmer(log(DNA+0.05) ~ load + (1|HM) + (1|RS_Name), peter.rapidaim.predicted) %>% summary()
# very sig


peter.rapidaim.predicted.rs <- subset(peter.rapidaim.predicted, RS_Name %in% rs.names.pbs)
peter.rapidaim.predicted.rs$RS_Name = factor(peter.rapidaim.predicted.rs$RS_Name, levels=rs.names.pbs)

# per DNA
peter.rapidaim.predicted.dna.lm = lmerTest::lmer(scale(log(DNA, base=10)) ~ RS_Name + (1|HM), peter.rapidaim.predicted.rs) %>% summary() %>% coef()%>% data.frame() %>% mutate(sig = ifelse(`Pr...t..` < 0.05, "*", ""))
peter.rapidaim.predicted.dna.lm$RS_Name = factor(gsub("RS_Name", "", row.names(peter.rapidaim.predicted.dna.lm)), levels=rs.names)
peter.rapidaim.predicted.dna.lm = peter.rapidaim.predicted.dna.lm %>%
  subset(!is.na(RS_Name)) %>%
  mutate(padj = p.adjust(`Pr...t..`, method="bonferroni")) %>%
  mutate(sig = ifelse(padj < 0.05, "*", ""))

peter.rapidaim.predicted.plot.dna = ggplot(peter.rapidaim.predicted.rs,
                                       aes(x=RS_Name, y=log(DNA,base=10)))+
  geom_point(size=2, color="black",
             position=position_jitter(seed=25, width=0.2))+
  geom_point(size=1.5, color="white",
             position=position_jitter(seed=25, width=0.2))+  
  geom_point(aes(color = RS_Name), alpha=0.3, size=1.5, 
             position=position_jitter(seed=25, width=0.2))+
  geom_boxplot(notch=T, outlier.shape=NA, width=0.5, aes(fill=RS_Name), alpha=0.6)+
  geom_text(data=peter.rapidaim.predicted.dna.lm, 
            aes(x = RS_Name, y = Inf, label=sig), vjust=1.2, size=7)+
  scale_color_manual(values=labelcolors$cols[c(10,1:9)])+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  theme_bw()+theme(legend.position="none",
                   axis.text.x=element_text(angle=45, hjust=1))+
  labs(x="", y="Log10(DNA ng/ul)")
peter.rapidaim.predicted.plot.dna

# per RS
peter.rapidaim.predicted.load.lm = lmerTest::lmer(scale(load) ~ RS_Name + (1|HM), peter.rapidaim.predicted.rs) %>% summary() %>% coef()%>% data.frame() %>% mutate(sig = ifelse(`Pr...t..` < 0.05, "*", ""))
peter.rapidaim.predicted.load.lm$RS_Name = factor(gsub("RS_Name", "", row.names(peter.rapidaim.predicted.load.lm)), levels=rs.names)
peter.rapidaim.predicted.load.lm = peter.rapidaim.predicted.load.lm %>%
  subset(!is.na(RS_Name)) %>%
  mutate(padj = p.adjust(`Pr...t..`, method="bonferroni")) %>%
  mutate(sig = ifelse(padj < 0.05, "*", ""))

peter.rapidaim.predicted.plot.load = ggplot(peter.rapidaim.predicted.rs,
                                           aes(x=RS_Name, y=log(DNA,base=10)))+
  geom_point(size=2, color="black",
             position=position_jitter(seed=25, width=0.2))+
  geom_point(size=1.5, color="white",
             position=position_jitter(seed=25, width=0.2))+  
  geom_point(aes(color = RS_Name), alpha=0.3, size=1.5, 
             position=position_jitter(seed=25, width=0.2))+
  geom_boxplot(notch=T, outlier.shape=NA, width=0.5, aes(fill=RS_Name), alpha=0.6)+
  geom_text(data=peter.rapidaim.predicted.load.lm, 
            aes(x = RS_Name, y = Inf, label=sig), vjust=1.2, size=7)+
  scale_color_manual(values=labelcolors$cols[c(10,1:9)])+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  theme_bw()+theme(legend.position="none",
                   axis.text.x=element_text(angle=45, hjust=1))+
  labs(x="", y="Predicted Microbial Load")
peter.rapidaim.predicted.plot.load

# plot
peter.rapidaim.predicted.plot | peter.rapidaim.predicted.plot.dna | peter.rapidaim.predicted.plot.load
# Cor = 0.257, p < 0.001


# :: Test on Stool data - Butyrogens ---------------------------------------------------

amplicon.data.gg <- readRDS("~/Documents/PhD/Dissertation/RS_Paper_4 (Longitudinal)/2023_06_26_longitudinal_analysis/all_trials_16S_min50k_mappingFile_240830_physeq_pooled_wTree_240830.Rds")
# subset to usable samples
amplicon.data.gg = phyloseq::subset_samples(amplicon.data.gg, SampleID %in% lsarp.phyloseq.df$SampleID)
# remove taxa not present
amplicon.data.gg = phyloseq::prune_taxa(phyloseq::taxa_sums(amplicon.data.gg) > 0, amplicon.data.gg)
# tax table
amplicon.data.gg.tax = data.frame(seq = phyloseq::refseq(amplicon.data.gg),
                                  ASV = rownames(data.frame(phyloseq::refseq(amplicon.data.gg))))
# apply RDP
amplicon.data.gg.tax.rdp <- dada2::assignTaxonomy(unique(as.character(amplicon.data.gg.tax$seq)), 
                                                "~/Downloads/rdp_train_set_16.fa.gz")
# format RDP to be compatible with MLP taxa naming
amplicon.data.gg.tax.df = amplicon.data.gg.tax.rdp %>% data.frame()
# make LCA
amplicon.data.gg.tax.df$LCA = ifelse(!is.na(amplicon.data.gg.tax.df$Genus), paste("g_", amplicon.data.gg.tax.df$Genus, sep=""),
                                       ifelse(!is.na(amplicon.data.gg.tax.df$Family), paste("uc_f_", amplicon.data.gg.tax.df$Family, sep=""),
                                              ifelse(!is.na(amplicon.data.gg.tax.df$Order), paste("uc_o_", amplicon.data.gg.tax.df$Order, sep=""),
                                                     ifelse(!is.na(amplicon.data.gg.tax.df$Class), paste("uc_c_", amplicon.data.gg.tax.df$Class, sep=""),
                                                            ifelse(!is.na(amplicon.data.gg.tax.df$Phylum), paste("uc_p_", amplicon.data.gg.tax.df$Phylum, sep=""),
                                                                   amplicon.data.gg.tax.df$Kingdom)))))
amplicon.data.gg.tax.df$seq = rownames(amplicon.data.gg.tax.df)
# add ASV
amplicon.data.gg.tax.df = merge(amplicon.data.gg.tax.df,
                                amplicon.data.gg.tax, by="seq")
# make OTU table
amplicon.data.gg.df = speedyseq::psmelt(amplicon.data.gg)
amplicon.data.gg.df$LCA = amplicon.data.gg.tax.df[match(amplicon.data.gg.df$OTU, amplicon.data.gg.tax.df$ASV),]$LCA

amplicon.data.gg.mat = amplicon.data.gg.df %>%
  group_by(Sample, LCA) %>%
  mutate(sum.abun = sum(Abundance)/50000) %>%
  group_by(standard.name) %>%
  mutate(med.abun = median(sum.abun)) %>%
  dplyr::select(standard.name, LCA, med.abun) %>% distinct() %>%
  reshape2::acast(standard.name ~ LCA, value.var="med.abun") %>% data.frame()
# done

# apply prediction
# add missing variables to testing data
missing.taxa.stool = colnames(model.16s.glmnet$trainingData)[!colnames(model.16s.glmnet$trainingData) %in%
                                                                      colnames(amplicon.data.gg.mat)] %>% sort()

amplicon.data.gg.mat.plus = cbind(amplicon.data.gg.mat, 
                           do.call(cbind, lapply(missing.taxa.stool, function(x){
                             new.data = data.frame(x = rep(0, times = nrow(amplicon.data.gg.mat)))
                             colnames(new.data) = x
                             new.data
                           }))) %>% data.frame()
pseudo.test.lsarp = min(amplicon.data.gg.mat.plus[amplicon.data.gg.mat.plus!=0])/2
amplicon.data.gg.mat.plus = log(amplicon.data.gg.mat.plus+pseudo.test.lsarp, base=10)

## APPLY MODEL
model.output.stool = predict(model.16s.rf, amplicon.data.gg.mat.plus) %>% data.frame()
colnames(model.output.stool)[1] = "load"
model.output.stool$standard.name = rownames(amplicon.data.gg.mat.plus)
hist(model.output.stool$load)

# get butyrogens
lsarp.phyloseq.median.butyrogens.load = merge(lsarp.phyloseq.median.butyrogens,
                                         model.output.stool, by="standard.name")


# Compare Placebo to Treatment
lmerTest::lmer(scale(load) ~ rs.days*group + (1|study_id), lsarp.phyloseq.median.butyrogens.load) %>% summary()
lmerTest::lmer(scale(log2(med.but)) ~ rs.days*group + (1|study_id), lsarp.phyloseq.median.butyrogens.load) %>% summary()
# p = 0.766, without adjusting for microbial load
lmerTest::lmer(scale(log2(med.but*10^load))~ rs.days*group + (1|study_id), lsarp.phyloseq.median.butyrogens.load) %>% summary()
# p = 0.641, WITH adjusting for microbial load

# Compare responders to non-responders
lmerTest::lmer(scale(load) ~ rs.days*responder + (1|study_id),
               subset(lsarp.phyloseq.median.butyrogens.load, group=="RS")) %>% summary()
lmerTest::lmer(scale(log2(med.but)) ~ rs.days*responder + (1|study_id),
               subset(lsarp.phyloseq.median.butyrogens.load, group=="RS")) %>% summary()
# p = 0.0211, without adjusting for microbial load
lmerTest::lmer(scale(log2(med.but*10^load)) ~ rs.days*responder + (1|study_id),
               subset(lsarp.phyloseq.median.butyrogens.load, group=="RS")) %>% summary()
# p = 0.0449, WITH adjusting for microbial load

lsarp.load.responders = ggplot(subset(lsarp.phyloseq.median.butyrogens.load, group=="RS"),
       #aes(x=rs.days, y=(residuals)))+
       aes(x=rs.days, y=log2(med.but*10^load)))+
  geom_line(aes(group=study_id), alpha=0.4)+
  geom_smooth(method="lm", color="black")+
  geom_point(color="black", size=2.5)+
  geom_point(color="white", size=1.5)+
  geom_point(aes(color=group), size=1.5)+
  scale_color_manual(values=labelcolors$cols[8])+
  theme_bw()+theme(legend.position="none")+
  facet_wrap(~responder)+
  labs(x="Days on RS", y="Butyrogens (Predicted Absolute)")
lsarp.load.responders

# no association with load
lmerTest::lmer(scale(log2(med.but*10^load)) ~ rs.days*group + (1|study_id),lsarp.phyloseq.median.butyrogens.load) %>% summary()

lsarp.load.placebo.treatment = ggplot(lsarp.phyloseq.median.butyrogens.load,
       #aes(x=rs.days, y=(residuals)))+
       aes(x=rs.days, y=log2(med.but*10^load)))+
  geom_line(aes(group=study_id), alpha=0.4)+
  geom_smooth(method="lm", color="black")+
  geom_point(color="black", size=2.5)+
  geom_point(color="white", size=1.5)+
  geom_point(aes(color=group), size=1.5, alpha=0.6)+
  scale_color_manual(values=labelcolors$cols[c(1,8)])+
  #ggrepel::geom_text_repel(aes(label=study_id))+
  theme_bw()+theme(legend.position="none")+
  facet_wrap(~group)+
  labs(x="Days on RS", y="Butyrogens (Predicted Absolute)")
lsarp.load.placebo.treatment


lsarp.load.placebo.treatment / lsarp.load.responders




# :: REVERSE MODEL --------------------------------------------------------

# TRAINING DATA
reverse.16s.train = read.csv("~/Downloads/Vandeputte_2017_16S.tsv", sep="\t") # supp table 14
dim(reverse.16s.train) # 95 samples
rownames(reverse.16s.train) <- reverse.16s.train$X
reverse.16s.train$X <- NULL
# convert to %
reverse.16s.train = sweep(reverse.16s.train, 1, rowSums(reverse.16s.train), FUN = "/")
reverse.16s.train
# add outcome
reverse.16s.train.outcome = read.csv("~/Downloads/41586_2017_BFnature24460_MOESM14_ESM.csv")
# add pseudocount and log transform
reverse.pseud = min(reverse.16s.train[reverse.16s.train !=0])/2
reverse.16s.train = log(reverse.16s.train[,colnames(reverse.16s.train) != "outcome"]+reverse.pseud, base=10)
reverse.16s.train$outcome = log(reverse.16s.train.outcome$Average.cell.count..per.gram.of.frozen.feces., base=10)

# train model

# train model with caret
t1 = Sys.time()
reverse.16s.xgb = caret::train(
  outcome ~ .,
  data = reverse.16s.train,
  trControl = caret::trainControl(method="repeatedcv", 
                                  number=10, 
                                  repeats = 1, 
                                  savePredictions = TRUE),
  #tuneGrid = tune_grid,
  method = "xgbTree",
  preProc = c("center", "scale")
)
t2 = Sys.time()
t2 - t1

# check performance
best <- reverse.16s.xgb$bestTune
pred <- reverse.16s.xgb$pred %>% filter(nrounds == best$nrounds & max_depth == best$max_depth & eta == best$eta & gamma == best$gamma & colsample_bytree == best$colsample_bytree & min_child_weight == best$min_child_weight & subsample == best$subsample)
cor.test(pred$pred,
         pred$obs, method = "spearman")
# Rho = 0.616

# :: reconstruct testing data ---------------------------------------------

reverse.16s.test.raw = read.csv("~/Downloads/2025_03_05_vandeputte_otu_unrarefied.csv")
reverse.16s.test.load = read.csv("~/Downloads/2025_03_05_vandeputte_load.csv")[,c("X", "Cell_count_per_gram", "ID_Number", "Day_Number")]
unique(reverse.16s.test.load$ID_Number) %>% length()
# n = 20 unique subjects
# remove NA
reverse.16s.test.raw = reverse.16s.test.raw[!is.na(reverse.16s.test.load$Cell_count_per_gram),]
reverse.16s.test.load = reverse.16s.test.load[!is.na(reverse.16s.test.load$Cell_count_per_gram),]
dim(reverse.16s.test.raw)
# rarefy taxa
reverse.16s.test.raw$X = NULL
set.seed(25)
reverse.16s.test.raw = rarefy_even_depth(phyloseq::otu_table(reverse.16s.test.raw, taxa_are_rows=F), 
                                     sample.size=10000,
                                     replace=F) %>% data.frame()
# convert to %
reverse.16s.test.raw <- sweep(reverse.16s.test.raw, 1, rowSums(reverse.16s.test.raw), FUN = "/")

# make predictions
# and lastly, remove samples in load that were removed in OTU table
reverse.16s.test.load = subset(reverse.16s.test.load, X %in% rownames(reverse.16s.test.raw))
nrow(reverse.16s.test.raw)
nrow(reverse.16s.test.load)

# add missing variables to testing data
missing.taxa = colnames(reverse.16s.xgb$trainingData)[!colnames(reverse.16s.xgb$trainingData) %in%
                                                                colnames(reverse.16s.test.raw)]
reverse.16s.test = cbind(reverse.16s.test.raw, 
                        do.call(cbind, lapply(missing.taxa, function(x){
                          new.data = data.frame(x = rep(0, times = nrow(reverse.16s.test.raw)))
                          colnames(new.data) = x
                          new.data
                        }))) %>% data.frame()
# pseud.test = min(reverse.16s.test[reverse.16s.test !=0])/2
# use training data's pseudocount
reverse.16s.test = log(reverse.16s.test[,colnames(reverse.16s.test) != "outcome"]+reverse.pseud, base=10)

plot(predict(reverse.16s.xgb, reverse.16s.test),
         log(reverse.16s.test.load$Cell_count_per_gram, base=10))
cor.test(predict(reverse.16s.xgb, reverse.16s.test),
         log(reverse.16s.test.load$Cell_count_per_gram, base=10))
# Cor = 0.42
# recall: xgboost; need to try other models


## data

reverse.16s.test

# :: loop through train + test --------------------------------------------


ml.models = c("ranger", "gbm", "xgbTree", "lm", "glmnet", "svmLinear", "svmRadial")

t0 <- Sys.time()
reverse.16s.iter = do.call(rbind, lapply(1:15, function(iter) {
  do.call(rbind, lapply(ml.models, function(model){
    #iter = 1
    #model = "ranger"
    print(paste0(iter, model))
    # train model with caret
    t1 = Sys.time()
    set.seed(iter)
    model.16s.ml = caret::train(
      outcome ~ .,
      data = reverse.16s.train,
      trControl = caret::trainControl(method="cv", 
                                      number=5, 
                                      savePredictions = TRUE),
      method = model,
      preProc = c("center", "scale")
    )
    t2 = Sys.time()
    t2 - t1
    
    # fix test data
    reverse.16s.test.processed = reverse.16s.test.raw
    rownames(reverse.16s.test.processed) <- reverse.16s.test.processed$X
    reverse.16s.test.processed$X <- NULL
    # convert to %
    reverse.16s.test.processed = sweep(reverse.16s.test.processed, 1, rowSums(reverse.16s.test.processed), FUN = "/")
    
    # add missing variables to testing data
    missing.taxa = colnames(model.16s.ml$trainingData)[!colnames(model.16s.ml$trainingData) %in%
                                                         colnames(reverse.16s.test.processed)]
    reverse.16s.test.processed = cbind(reverse.16s.test.processed, 
                                      do.call(cbind, lapply(missing.taxa, function(x){
                                        new.data = data.frame(x = rep(0, times = nrow(reverse.16s.test.processed)))
                                        colnames(new.data) = x
                                        new.data
                                      }))) %>% data.frame()
    #pseudo.test.reverse.loop = min(reverse.16s.test.processed[reverse.16s.test.processed!=0])/2
    reverse.16s.test.processed = log(reverse.16s.test.processed+reverse.pseud, base=10)
    
    # predict test data
    
    model.output = predict(model.16s.ml, reverse.16s.test.processed)
    
    # save
    data.frame(true = log(reverse.16s.test.load$Cell_count_per_gram, base=10),
               pred = model.output,
               iter = iter,
               model = model,
               time = as.numeric(t2 - t1))
  }))}))
t3 <- Sys.time() # 22 min
t3 - t0
saveRDS(reverse.16s.iter, "~/Documents/PhD/For others/2025_03_05_mlp_16s_test_reverse.Rds")
reverse.16s.iter.df <- reverse.16s.iter


# plot raw vals
ggplot(reverse.16s.iter.df,
       aes(x=true, y=pred))+
  geom_smooth(aes(group = iter), method="lm")+
  geom_point()+
  facet_wrap( ~ model, scales="free")

# evaluate using Spearman

reverse.16s.iter.eval = do.call(rbind, lapply(1:15, function(seed) {
  do.call(rbind, lapply(ml.models, function(ml){
    print(paste0(ml, seed))
    data.subset = subset(reverse.16s.iter.df, iter == seed & model == ml) 
    #
    spear = cor.test(data.subset$true,
                     data.subset$pred, method="spearman")$estimate
    pear = cor.test(data.subset$true,
                    data.subset$pred, method="spearman")$estimate
    if(sd(data.subset$pred)==0){
      r2 = NA
    }else{
      r2 = summary(lm(scale(pred) ~ scale(true), data.subset))$r.squared
    }
    # save
    data.frame(model = ml,
               iter = seed,
               spearman = spear,
               pearson = pear,
               r2 = r2)
  }))}))


# ```plot: model evaluation -----------------------------------------------

# plot time
ggplot(reverse.16s.iter.df[,c("iter", "model", "time")] %>% distinct())+
  geom_boxplot(aes(x = reorder(model, time), y=time))+
  geom_point(aes(x = reorder(model, time), y=time))+
  theme_bw()+
  labs(x="Model", y="Time elapsed (sec)")
# RF is the fastest
# glmnet is not bad; much better than svmLinear

# plot R2
ggplot(reverse.16s.iter.eval)+
  geom_boxplot(aes(x = model, y=r2))+
  geom_point(aes(x = model, y=r2))+
  theme_bw()+
  labs(x="Model", y="R-Squared")

# plot pearson
reverse.pearson.plot = ggplot(reverse.16s.iter.eval)+
  geom_boxplot(aes(x = model, y=pearson))+
  geom_point(aes(x = model, y=pearson), size=2.5, color="black")+
  geom_point(aes(x = model, y=pearson), size=1.5, color="white")+
  geom_point(aes(x = model, y=pearson, color=model), size=1.5, alpha=0.6)+
  theme_bw()+theme(legend.position="none")+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x="Model", y="Pearson Correlation")
# why is there no variation with glmnet, etc?
# it is deterministic

reverse.time.plot = ggplot(subset(reverse.16s.iter.df, model %in% ml.models)[,c("iter", "model", "time")] %>% distinct())+
  geom_boxplot(aes(x = model, y=time))+
  geom_point(aes(x = model, y=time), size=2.5, color="black")+
  geom_point(aes(x = model, y=time), size=1.5, color="white")+
  geom_point(aes(x = model, y=time, color=model), size=1.5, alpha=0.6)+
  theme_bw()+theme(legend.position="none")+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x="Model", y="Seconds Elapsed")

reverse.pearson.plot/reverse.time.plot


# WRAPPER -----------------------------------------------------------------

mlp <- function(
    input = phyloseq_object, # rarefied or %
    tax = NULL,
    assign = T,
    predict = T){
  # assign tax
  if(assign == T){
    if(is.null(tax)){
      print(paste("Checking tax table"))
    ps.tax = data.frame(phyloseq::tax_table(input))
    ps.tax$seq = rownames(ps.tax)
    if(sum(grepl("ASV", rownames(ps.tax)))>1){ # this detects whether tax_table or refseq table were used
      ps.tax = data.frame(seq = phyloseq::refseq(input),
                          ASV = rownames(ps.tax))
    }
    # apply RDP
    print(paste("Re-assigning taxonomy"))
    ps.tax.rdp <- dada2::assignTaxonomy(unique(as.character(ps.tax$seq)), 
                                                    "./rdp_train_set_16.fa.gz")
    print(paste("Collapsing to lowest common annotated ancestor"))
    # format RDP to be compatible with MLP taxa naming
    ps.tax.rdp.df = ps.tax.rdp %>% data.frame()
    # make LCA
    ps.tax.rdp.df$LCA = ifelse(!is.na(ps.tax.rdp.df$Genus), paste("g_", ps.tax.rdp.df$Genus, sep=""),
                                           ifelse(!is.na(ps.tax.rdp.df$Family), paste("uc_f_", ps.tax.rdp.df$Family, sep=""),
                                                  ifelse(!is.na(ps.tax.rdp.df$Order), paste("uc_o_", ps.tax.rdp.df$Order, sep=""),
                                                         ifelse(!is.na(ps.tax.rdp.df$Class), paste("uc_c_", ps.tax.rdp.df$Class, sep=""),
                                                                ifelse(!is.na(ps.tax.rdp.df$Phylum), paste("uc_p_", ps.tax.rdp.df$Phylum, sep=""),
                                                                       ps.tax.rdp.df$Kingdom)))))
    ps.tax.rdp.df$seq = rownames(ps.tax.rdp.df)
    # make OTU table
    input.df = speedyseq::psmelt(input)
    # replace ASV name with OTU seq, if refseq was used
    if(sum(grepl("ASV", rownames(ps.tax)))>1){ # this detects whether tax_table or refseq table were used
      input.df$OTU = ps.tax[match(input.df$OTU, ps.tax$ASV),]$seq
    }
    input.df$LCA = ps.tax.rdp.df[match(input.df$OTU, ps.tax.rdp.df$seq),]$LCA
    # collapse to LCA
    input.mat = input.df %>%
      # mutate(Sample = row.names(.)) %>%
      group_by(Sample, LCA) %>%
      mutate(sum.abun = sum(Abundance)) %>%
      dplyr::select(Sample, LCA, sum.abun) %>% distinct() %>%
      reshape2::acast(Sample ~ LCA, value.var="sum.abun") %>% data.frame()
    # normalize to %
    input.mat = input.mat / rowSums(input.mat)
    # return
    return(input.mat)
    }
  }
  if(predict == TRUE){
  # load model
    print(paste("Processing data"))
  model.16s.rf = readRDS("./2025_03_06_mlp_rf.Rds")
  # extract pseudo
  pseudo = min(10^model.16s.rf$trainingData)
  # add missing variables
  missing.taxa = colnames(model.16s.rf$trainingData)[!colnames(model.16s.rf$trainingData) %in%
                                                             colnames(input.mat)]
  input.mat.processed = cbind(input.mat, 
                          do.call(cbind, lapply(missing.taxa, function(x){
                            new.data = data.frame(x = rep(0, times = nrow(input.mat)))
                            colnames(new.data) = x
                            new.data
                          }))) %>% data.frame()
  
  input.mat.processed$X = NULL
  input.mat.processed = log(input.mat.processed+pseudo, base=10)
  print(paste("Making predictions"))
  model.output = data.frame(
    sample = rownames(input.mat),
    load = predict(model.16s.rf, input.mat.processed))
  print(paste("Done"))
  return(model.output)
  }
}
    
setwd("~/Documents/PhD/For others/mlp_pd")

test = mlp(
  input = amplicon.data.gg, 
  tax = NULL,
  assign = F,
  predict = T)



# :: VALIDATE -------------------------------------------------------------

# use my downloaded data to reconstruct and evaluate model

reverse.16s.test.raw
reverse.16s.test.load

reverse.16s.validate = reverse.16s.test.raw
pseud.validate = min(reverse.16s.validate[reverse.16s.validate!=0])/2
reverse.16s.validate = log(reverse.16s.validate+pseud.validate, base=10)
reverse.16s.validate$load = log(reverse.16s.test.load$Cell_count_per_gram, base=10)

# train model with caret
t1 = Sys.time()
model.16s.rf.validate = caret::train(
  load ~ .,
  data = reverse.16s.validate,
  trControl = caret::trainControl(method="cv", 
                                  number=5, 
                                  savePredictions = TRUE),
  #tuneGrid = tune_grid,
  method = "ranger",
  preProc = c("center", "scale")
)
t2 = Sys.time()
t2 - t1

# TRAINING DATA
reverse.16s.validate.test = read.csv("~/Downloads/Vandeputte_2017_16S.tsv", sep="\t") # supp table 14
rownames(reverse.16s.validate.test) <- reverse.16s.validate.test$X
reverse.16s.validate.test$X <- NULL
# convert to %
reverse.16s.validate.test = sweep(reverse.16s.validate.test, 1, rowSums(reverse.16s.validate.test), FUN = "/")
# add missing variables to testing data
missing.taxa.validate = colnames(model.16s.rf.validate$trainingData)[!colnames(model.16s.rf.validate$trainingData) %in%
                                                     colnames(reverse.16s.validate.test)]
reverse.16s.validate.test = cbind(reverse.16s.validate.test, 
                        do.call(cbind, lapply(missing.taxa.validate, function(x){
                          new.data = data.frame(x = rep(0, times = nrow(reverse.16s.validate.test)))
                          colnames(new.data) = x
                          new.data
                        }))) %>% data.frame()
# add pseudocount and log transform
reverse.16s.validate.test = log(reverse.16s.validate.test+pseud.validate, base=10)

# apply model
model.16s.validate.predictions = predict(model.16s.rf.validate, reverse.16s.validate.test)
# add outcome
reverse.16s.train.outcome = read.csv("~/Downloads/41586_2017_BFnature24460_MOESM14_ESM.csv")
# correlate
cor.test(model.16s.validate.predictions,
         log(reverse.16s.train.outcome$Average.cell.count..per.gram.of.frozen.feces., base=10))
# Cor = 0.608

# consider this validated



# :: Feature Importance ---------------------------------------------------

# rebuilt ranger
set.seed(25)
model.16s.ranger = ranger::ranger(outcome ~ .,
                                 data = model.16s.xgb.training.data,
                                 importance = "impurity")

# feature importances
model.16s.ranger.importance = data.frame(imp = ranger::importance(model.16s.ranger)) %>% arrange(imp)

# check TimbR
class(model.16s.ranger)
timbr.output = timbR::select_trees(
  rf = model.16s.ranger,
  distance.matrix = timbR::measure_distances(rf = model.16s.ranger, metric = "splitting variables"),
  num.trees = 1)

tree_df <- ranger::treeInfo(timbr.output, tree = 1)
