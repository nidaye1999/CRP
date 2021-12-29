# 0:female/no; 1:male/yes
# https://rpubs.com/aryn999/LinearRegressionAssumptionsAndDiagnosticsInR
# https://data.library.virginia.edu/diagnostic-plots/
# log-transformed data; residual-vs-leverage; QQ plot; VIF; independence


suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(VIM))
suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(lmtest))

set.seed(60)

# load data
data=read.csv("data.csv", check.names=F, row.names=1, stringsAsFactors=TRUE)


# plot missing value distribution
temp=data[,colnames(data)[colSums(is.na(data)) > 0]]
png("results/linear model/plots/missing.png",width=8, height=10, units = "in", res = 300)
aggr(temp, col=c('navyblue','yellow'),numbers=TRUE, sortVars=TRUE,labels=colnames(temp), cex.axis=.8,cex.numbers=0.5,gap=3, ylab=c("Proportion of Missings for Single Attribute","Proportion of Missings for Attribute Combination"), oma = c(11.7, 4.5, 1, 0))
dev.off()
rm(temp)

data=na.omit(data)

# binary coding
binary.coding.variables=which(apply(data, 2, function(x) length(unique(x))) %in% c(2))
binary.variable=colnames(data)[binary.coding.variables]
data$Gender=factor(data$Gender, exclude = NULL, levels=c(0,1),labels=c("Female", "Male"))
data$'C-reactive Protein Detection Range'=factor(data$'C-reactive Protein Detection Range', exclude = NULL, levels=c(0,1),labels=c("LOD (0.5-23)", "LOD (0.8-20)"))
for (item in binary.variable[!binary.variable %in% c("Gender", "C-reactive Protein Detection Range")]){
  data[,item]=factor(data[,item], exclude = NULL, levels=c(0,1),labels=c("No", "Yes"))
}
rm(item)

sum(round(10^(data$`lg C-reactive Protein`),2)==0.50)
sum(round(10^(data$`lg C-reactive Protein`),2)==0.80)


# log transformation
for (item in c("lg Income", "lg Alcohol Use", "lg Caffeine Use")){
  data[[item]]=log10(data[[item]]+1)
  # colnames(data)[colnames(data) == item] = paste("log", item, sep=".")
}
rm(item)

# get continuous variable
continuous.variable=colnames(data)[-binary.coding.variables]
rm(binary.coding.variables)

# pair plot
# p=ggpairs(data, diag=list(continuous=wrap("barDiag", bins=50)), upper = list(continuous = wrap("cor", size = 3), combo=wrap("box")), lower = list(continuous = wrap("points", size = 0.1), combo=wrap("facethist", bins=50)))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave("results/linear model/plots/original.png", p, dpi=400, height = 25, width = 25, limitsize = F)
# rm(p)
temp=aggregate(data$`lg C-reactive Protein`, list(data$Gender), FUN=mean)
p<-ggplot(data, aes(x=`lg C-reactive Protein`, fill=Gender, color=Gender)) +
  geom_histogram(position="identity", alpha=0.5, bins = 50)+
  geom_vline(data=temp, aes(xintercept=x, color=Group.1),linetype="dashed")
ggsave(file="results/crp_gender_hist.png", p, width=4, height=2.5, units = "in", dpi=300)

rm(temp)


# generate table1
suppressPackageStartupMessages(library(table1))
data.table1=data
data.table1$Education=factor(data.table1$Education, exclude = NULL, levels=c(1,2,3,4,5,6,7),labels=c("Less than seven years of school", "Junior high school (7th, 8th, 9th)", "Some high school (10th, 11th)", "High school graduate (including equivalency exam)", "Some college or technical school (at least one year)", "College graduate", "Graduate professional training (Masters or above)"))
table1(~ .|Gender, data=data.table1[,colnames(data.table1)[c(1:2,4:21, 3)]])
table1(~ .|Gender, data=data.table1[,colnames(data.table1)[c(22:37,3)]])
rm(data.table1)

# t-test for numerical variables
t.test.results=data.frame(row.names = c("Female", "Male", "t-statistics", "df", "p-value"))
for (item in continuous.variable){
  temp=t.test(data[[item]]~Gender, data=data)
  t.test.results[,item]=c(temp$estimate[[1]],temp$estimate[[2]],temp$statistic[[1]], temp$parameter[[1]], temp$p.value)
}
write.csv(t.test.results, "results/linear model/excels/t_test_gender.csv")

# chi-squared test for binary variables
chisq.results=data.frame(row.names = c("Chi-squared", "df", "p-value"))
for (item in binary.variable[-c(1)]){
  temp=chisq.test(data[[item]], data$Gender)
  chisq.results[,item]=c(temp$statistic[[1]], temp$parameter[[1]], temp$p.value)
}
write.csv(chisq.results, "results/linear model/excels/chi_squared_gender.csv")
rm(t.test.results, chisq.results, item, temp)

# define linear regression diagnostics function
lm.diagnostics=function(name, data, continuous.variable, id.n=3){
  print(name)
  data[,continuous.variable]=scale(data[,continuous.variable])
  # print(describe(data))
  model=lm(data$'lg C-reactive Protein'~., data=data[2:dim(data)[2]])
  print(summary(model))
  write.csv(summary(model)$coefficients, paste("results/linear model/excels/standardized_beta_",name,".csv", sep=""))
  # Linearity of the data
  png(paste("results/linear model/plots/fit_",name,".png",sep=""),width=5.5, height=5.5, units = "in", res = 150)
  par(mgp=c(2.5, 1, 0))
  plot(model, 1, cex.caption=1.75, cex.axis=1.25, cex.lab=1.45, sub.caption=" ")
  dev.off()
  # Normality of Residuals
  png(paste("results/linear model/plots/QQ_plot_",name,".png",sep=""),width=5.5, height=5.5, units = "in", res = 150)
  par(mgp=c(2.5, 1, 0))
  plot(model, 2, cex.caption=1.75, cex.axis=1.25, cex.lab=1.45, sub.caption=" ")
  dev.off()
  # Homoscedasticity Assumption
  png(paste("results/linear model/plots/Sca_Loc_",name,".png",sep=""),width=5.5, height=5.5, units = "in", res = 150)
  par(mgp=c(2.5, 1, 0))
  plot(model, 3, cex.caption=1.75, cex.axis=1.25, cex.lab=1.45, sub.caption=" ")
  dev.off()
  print(bptest(model))
  # High leverage points
  png(paste("results/linear model/plots/Res_Lev_",name,".png",sep=""),width=5.5, height=5.5, units = "in", res = 150)
  par(mgp=c(2.5, 1, 0))
  plot(model, 5, id.n = id.n, cex.caption=1.75, cex.axis=1.25, cex.lab=1.45, sub.caption=" ")
  dev.off()
  # Durbin Watson Test for Autocorrelation
  print(durbinWatsonTest(model))
  # detecting multicollinearity
  print(vif(model))
  write.csv(vif(model), paste("results/linear model/excels/vif_",name,".csv", sep=""))
  # plot betas
  coef=tidy(model, conf.int = TRUE)
  coef$term=str_replace(coef$term, "Yes|Male|LOD \\(0.8-20\\)", "")
  p=ggplot(coef, aes(term, estimate)) +
    geom_hline(yintercept = 0, size = I(1.1), color = I("black"), linetype=2, alpha=0.5) +
    geom_point(stat="identity", color="black", size=3) +
    geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2,position=position_dodge(.9)) +
    theme(axis.text.x=element_text(color = "black", size=16, angle=0),
          axis.text.y=element_text(color = "black", size=16, angle=0),
          axis.title=element_text(size=20,face="bold"),
          plot.title = element_text(size = 24, face = "bold"))+
    xlab("Variable")+ylab("Standardized beta")+ggtitle(name)+scale_x_discrete(limits = rev(coef$term), labels = str_replace(str_replace(str_replace(str_replace_all(rev(coef$term), "`", ""), "Yes", "_Yes"), "LOD", "_LOD"), "Male", "_Male"))+coord_flip(ylim=c(-1.5, 1.5))
  ggsave(paste("results/linear model/plots/beta_",name,".png",sep=""), p, width = 10, height = 12, dpi = 200, units = "in")
  # corrplot
  png(paste("results/linear model/plots/corr_", name, ".png", sep=""), width=7.5, height=7.5, units = "in", res = 200)
  corrplot.mixed(cor(data[,continuous.variable]), mar=c(0,0,0,0), tl.pos="lt", tl.cex=.9, lower.col = "black", number.cex = .75, upper = "ellipse")
  mtext(name, at=0, line=0, cex=1.5)
  dev.off()
#
#   png(paste("results/linear model/plots/ellipse_corr_", name, ".png", sep=""), width=7.5, height=7.5, units = "in", res = 200)
#   corrplot(cor(data[,continuous.variable]), mar=c(0,0,0,0), tl.pos="lt", tl.cex=.9, method = "ellipse")
#   mtext(name, at=0, line=0, cex=1.5)
#   dev.off()
#
#   png(paste("results/linear model/plots/number_corr_", name, ".png", sep=""), width=7.5, height=7.5, units = "in", res = 200)
#   corrplot(cor(data[,continuous.variable]), mar=c(0,0,0,0), tl.pos="lt", tl.cex=.9, method = "number", number.cex = .75)
#   mtext(name, at=0, line=0, cex=1.5)
#   dev.off()
  # return(coef)
}

# run lm.diagnostics for all, female and male subjects
temp=data
temp$Contraceptive=NULL
lm.diagnostics("All", temp, continuous.variable)
temp=data[data$Gender=="Female",]
temp$Gender=NULL
lm.diagnostics("Female", temp, continuous.variable)
temp=data[data$Gender=="Male",]
temp$Gender=NULL
temp$Contraceptive=NULL
lm.diagnostics("Male", temp, continuous.variable)
rm(temp, lm.diagnostics)


# investigate confounding effect: CRP, gender, BMI, % body fat and waist-hip ratio
library(ggpubr)
formula = y ~ x
crp.gender.obesity=data[,c(1,3,15:17)]
# crp.gender.obesity[,c(1,3,4,5)]=scale(crp.gender.obesity[,c(1,3,4,5)])
p=ggplot(crp.gender.obesity, aes(x=crp.gender.obesity$'Body Mass Index', y=crp.gender.obesity$'lg C-reactive Protein', color=Gender))+geom_point(size=0.6)+
  geom_smooth(aes(fill = Gender, color = Gender), method=lm, formula = formula)+
  stat_regline_equation(aes(label =  paste(..eq.label..)), formula = formula, show.legend = FALSE)+
  ylim(min(crp.gender.obesity$'lg C-reactive Protein')-0.3, max(crp.gender.obesity$'lg C-reactive Protein')+0.4)+
  xlab('Body Mass Index')+ylab('lg C-reactive Protein')
ggsave("results/linear model/plots/crp_bmi_gender.png", p, width = 4.5, height = 3.5, dpi=300)

p=ggplot(crp.gender.obesity, aes(x=crp.gender.obesity$'Percent Body Fat', y=crp.gender.obesity$'lg C-reactive Protein', color=Gender))+geom_point(size=0.6)+
  geom_smooth(aes(fill = Gender, color = Gender), method=lm, formula = formula)+
  stat_regline_equation(aes(label =  paste(..eq.label..)), formula = formula, show.legend = FALSE)+
  ylim(min(crp.gender.obesity$'lg C-reactive Protein')-0.3, max(crp.gender.obesity$'lg C-reactive Protein')+0.4)+
  xlab('Percent Body Fat')+ylab('lg C-reactive Protein')
ggsave("results/linear model/plots/crp_body_fat_gender.png", p, width = 4.5, height = 3.5, dpi=300)

p=ggplot(crp.gender.obesity, aes(x=crp.gender.obesity$'Waist-to-hip Ratio', y=crp.gender.obesity$'lg C-reactive Protein', color=Gender))+geom_point(size=0.6)+
  geom_smooth(aes(fill = Gender, color = Gender), method=lm, formula = formula)+
  stat_regline_equation(aes(label =  paste(..eq.label..)), formula = formula, show.legend = FALSE)+
  ylim(min(crp.gender.obesity$'lg C-reactive Protein')-0.3, max(crp.gender.obesity$'lg C-reactive Protein')+0.4)+
  xlab('Waist-to-hip Ratio')+ylab('lg C-reactive Protein')
ggsave("results/linear model/plots/crp_waist_hip_gender.png", p, width = 4.5, height = 3.5, dpi=300)
rm(p, formula, crp.gender.obesity)

############################################################################################################
# Machine learning
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(source('caretStack_FUNC.R'))


# parameter for NCV
n_folder=10
tuneL=7
cpu=detectCores()-2


get.pred.grid=function(matrix_list, variable){
  x_temp=list()
  for (i in 1:length(matrix_list)){
    x_temp[[i]]=(matrix_list[[i]][,1])
  }
  x_temp=data.frame(sort(unique(unlist(x_temp))))
  colnames(x_temp)=c(variable)
  return(x_temp)
}

plot_summ=function(preds){
  require(plotrix)
  train=array(0,dim=c(dim(modelPerf.summ(preds[[1]]$prediction)$train.by.fold)[[1]],1,length(preds)))
  test=array(0,dim=c(dim(modelPerf.summ(preds[[1]]$prediction)$test)[[1]],dim(modelPerf.summ(preds[[1]]$prediction)$test)[[2]],length(preds)))
  for (i in 1:length(preds)){
    perf=modelPerf.summ(preds[[i]]$prediction)
    train[,,i]=perf$train.by.fold$m
    test[,,i]=perf$test
  }
  perf$train.by.fold$m=apply(train, c(1,2), mean)
  perf$train.by.fold$se=apply(train, c(1,2), std.error)
  perf$test[,]=apply(test, c(1,2), mean)
  metrics = colnames(perf$test)
  p=list()
  for (j in 1:length(metrics)){
    p[[j]] = print(plot.perf(perf, metrics=metrics[j]))
  }
  fig=grid.arrange(grobs=p, ncol=3, top="Overall")
  return(fig)
}


fitControl = trainControl(method = 'cv',
                           number = n_folder,
                           search = 'grid',
                           summaryFunction = defaultSummary,
                           selectionFunction = 'oneSE',
                           savePredictions = 'final')

# factor to numeric
data.numeric=data
for (item in binary.variable){
  data.numeric[,item]=as.numeric(data.numeric[,item])-1
}
rm(item)


# all subjects
temp=data.numeric
temp$Contraceptive=NULL
temp[,continuous.variable[-c(1)]]=scale(temp[,continuous.variable[-c(1)]])
flds = createFolds(temp$'C-reactive Protein', k = n_folder, list = TRUE, returnTrain = FALSE)
models=list()
preds=list()
fig=list()
cl = makePSOCKcluster(cpu)
registerDoParallel(cl)
for (i in 1:n_folder){
  training = temp[-flds[[i]],]
  testing  = temp[flds[[i]],]
  model = caretModels(training, resp.var='C-reactive Protein', control = fitControl,
                       preProc.opt = NULL, tuneL = tuneL, metric = 'RMSE',
                       methods = c('lm', 'pcr', 'svmRadial', 'ranger'))
  models[[i]] = model
  preds[[i]] = PredVal(model, testing, resp.var='C-reactive Protein')

  perf = modelPerf.summ(preds[[i]]$prediction)
  metrics = colnames(perf$test)
  p = list()
  for (j in 1:length(metrics)){
    p[[j]] = print(plot.perf(perf, metrics=metrics[j]))
  }
  fig[[i]]=grid.arrange(grobs=p, ncol=3, top=paste(toString(i), "th outer loop", sep=""))
  print(i)
}
stopCluster(cl)
fig[[n_folder+1]]=plot_summ(preds)
ggsave(file="results/machine learning/plots/all.png", arrangeGrob(grobs=fig, ncol=1), width=10, height=16*n_folder/5, units = "in")
saveRDS(models, "results/machine learning/models/all.rds")

# female
temp = data.numeric[data.numeric$Gender==0,]
temp$Gender=NULL
temp[,continuous.variable[-c(1)]]=scale(temp[,continuous.variable[-c(1)]])
flds = createFolds(temp$'C-reactive Protein', k = n_folder, list = TRUE, returnTrain = FALSE)
models=list()
preds=list()
fig=list()
cl = makePSOCKcluster(cpu)
registerDoParallel(cl)
for (i in 1:n_folder){
  training = temp[-flds[[i]],]
  testing  = temp[flds[[i]],]
  model = caretModels(training, resp.var='C-reactive Protein', control = fitControl,
                       preProc.opt = NULL, tuneL = tuneL, metric = 'RMSE',
                       methods = c('lm', 'pcr', 'svmRadial', 'ranger'))
  models[[i]] = model
  preds[[i]] = PredVal(model, testing, resp.var='C-reactive Protein')

  perf = modelPerf.summ(preds[[i]]$prediction)
  metrics = colnames(perf$test)
  p = list()
  for (j in 1:length(metrics)){
    p[[j]] = print(plot.perf(perf, metrics=metrics[j]))
  }
  fig[[i]]=grid.arrange(grobs=p, ncol=3, top=paste(toString(i), "th outer loop", sep=""))
  print(i)
}
stopCluster(cl)
fig[[n_folder+1]]=plot_summ(preds)
ggsave(file="results/machine learning/plots/female.png",arrangeGrob(grobs=fig, ncol=1), width=10, height=16*n_folder/5, units = "in")
saveRDS(models, "results/machine learning/models/female.rds")

# male
temp = data.numeric[data.numeric$Gender==1,]
temp$Gender=NULL
temp$Contraceptive=NULL
temp[,continuous.variable[-c(1)]]=scale(temp[,continuous.variable[-c(1)]])
flds = createFolds(temp$'C-reactive Protein', k = n_folder, list = TRUE, returnTrain = FALSE)
models=list()
preds=list()
fig=list()
cl = makePSOCKcluster(cpu)
registerDoParallel(cl)
for (i in 1:n_folder){
  training = temp[-flds[[i]],]
  testing  = temp[flds[[i]],]
  model = caretModels(training, resp.var='C-reactive Protein', control = fitControl,
                       preProc.opt = NULL, tuneL = tuneL, metric = 'RMSE',
                       methods = c('lm', 'pcr', 'svmRadial', 'ranger'))
  models[[i]] = model
  preds[[i]] = PredVal(model, testing, resp.var='C-reactive Protein')

  perf = modelPerf.summ(preds[[i]]$prediction)
  metrics = colnames(perf$test)
  p = list()
  for (j in 1:length(metrics)){
    p[[j]] = print(plot.perf(perf, metrics=metrics[j]))
  }
  fig[[i]]=grid.arrange(grobs=p, ncol=3, top=paste(toString(i), "th outer loop", sep=""))
  print(i)
}
stopCluster(cl)
fig[[n_folder+1]]=plot_summ(preds)
ggsave(file="results/machine learning/plots/male.png",arrangeGrob(grobs=fig, ncol=1), width=10, height=16*n_folder/5, units = "in")
saveRDS(models, "results/machine learning/models/male.rds")
rm(cl, fitControl, i, cpu, fig, j, models, preds, temp)

# partial dependent plot
library(pdp)
data.numeric=data
for (item in binary.variable){
  data.numeric[,item]=as.numeric(data.numeric[,item])-1
}
rm(item)

scale_info=list()

temp=attributes(scale(data.numeric[,continuous.variable[-c(1)]]))
temp=rbind(temp$`scaled:center`,temp$`scaled:scale`)
colnames(temp)[colnames(temp) == 'lg Income'] = "Income"
colnames(temp)[colnames(temp) == 'lg Alcohol Use'] = "Alcohol Use"
colnames(temp)[colnames(temp) == 'lg Caffeine Use'] = "Caffeine Use"
colnames(temp)[colnames(temp) == 'Waist-to-hip Ratio'] = "Waist to hip Ratio"
row.names(temp)=c("mean", "std")
scale_info[["All"]]=temp
temp=attributes(scale(data.numeric[data.numeric$Gender==0,continuous.variable[-c(1)]]))
temp=rbind(temp$`scaled:center`,temp$`scaled:scale`)
colnames(temp)[colnames(temp) == 'lg Income'] = "Income"
colnames(temp)[colnames(temp) == 'lg Alcohol Use'] = "Alcohol Use"
colnames(temp)[colnames(temp) == 'lg Caffeine Use'] = "Caffeine Use"
colnames(temp)[colnames(temp) == 'Waist-to-hip Ratio'] = "Waist to hip Ratio"
row.names(temp)=c("mean", "std")
scale_info[["Female"]]=temp
temp=attributes(scale(data.numeric[data.numeric$Gender==1,continuous.variable[-c(1)]]))
temp=rbind(temp$`scaled:center`,temp$`scaled:scale`)
colnames(temp)[colnames(temp) == 'lg Income'] = "Income"
colnames(temp)[colnames(temp) == 'lg Alcohol Use'] = "Alcohol Use"
colnames(temp)[colnames(temp) == 'lg Caffeine Use'] = "Caffeine Use"
colnames(temp)[colnames(temp) == 'Waist-to-hip Ratio'] = "Waist to hip Ratio"
row.names(temp)=c("mean", "std")
scale_info[["Male"]]=temp


pdp_summ=function(models, title, scale_info){
  dir.create(file.path(paste("results/machine learning/plots/", title, sep="")), showWarnings = FALSE)
  p=list()
  index=1
  for (variable in models[[1]]$ranger$finalModel$xNames){
    print(variable)
    temp=list()
    for (i in 1:length(models)){
      temp[[i]]=partial(models[[i]]$ranger, pred.var = variable, plot = F, rug = TRUE)
      temp[[i]]$group <- rep(i,nrow(temp[[i]]))
      temp[[i]]=as.matrix(temp[[i]])
    }
    x_temp=get.pred.grid(temp, variable)
    temp=list()
    for (i in 1:length(models)){
      temp[[i]]=partial(models[[i]]$ranger,   pred.var = variable, plot = F, rug = TRUE, pred.grid = x_temp)
      temp[[i]]$group <- rep(i,nrow(temp[[i]]))
      temp[[i]]=as.matrix(temp[[i]])
    }
    temp_summ=data.frame(matrix(ncol = 3, nrow = dim(temp[[1]])[[1]]))
    colnames(temp_summ)=c(variable, "yhat", "se")
    temp_summ[,1]=x_temp
    temp_matrix=data.frame(matrix(ncol = length(models), nrow=dim(temp_summ)[[1]]))
    for (i in 1:length(models)){
      temp_matrix[,i]=temp[[i]][,2]
    }
    # temp_matrix=temp_matrix*scale_info["std", "log.CRP"]+scale_info["mean", "log.CRP"]
    temp_summ[,2]=apply(temp_matrix,1, mean, na.rm = TRUE)
    temp_summ[,3]=apply(temp_matrix,1, sd, na.rm = TRUE)
    if (variable %in% colnames(scale_info)){
      if (variable %in% c("Income", "Alcohol Use", "Caffeine Use")){
        variable_1=paste("lg", variable)
      }
      else if(variable=="Waist to hip Ratio"){
        variable_1="Waist-to-hip Ratio"
      }
      else{
        variable_1=variable
      }
      plt_summ=ggplot(data=temp_summ, aes(x = temp_summ[,1]*scale_info["std", variable]+scale_info["mean", variable], y = yhat))+geom_line(size = 1)+geom_ribbon(aes(ymin=temp_summ$yhat-1.96*temp_summ$se, ymax=temp_summ$yhat+1.96*temp_summ$se), alpha=0.3)+xlab(variable_1)+ylab("lg C-reactive Protein")+ylim(0.00,0.80)+theme(axis.text.x=element_text(color = "black", size=12, angle=0),axis.text.y=element_text(color = "black", size=12, angle=0),axis.title=element_text(size=15,face="bold"))
      temp_matrix[[variable]]=x_temp*scale_info["std", variable]+scale_info["mean", variable]
    }else{
      plt_summ=ggplot(data=temp_summ, aes(x = temp_summ[,1], y = yhat))+geom_line(size = 1)+geom_ribbon(aes(ymin=temp_summ$yhat-1.96*temp_summ$se, ymax=temp_summ$yhat+1.96*temp_summ$se), alpha=0.3)+xlab(variable)+ylab("lg C-reactive Protein")+ylim(0.00,0.80)+theme(axis.text.x=element_text(color = "black", size=12, angle=0),axis.text.y=element_text(color = "black", size=12, angle=0),axis.title=element_text(size=15,face="bold"))
      temp_matrix[[variable]]=x_temp
    }
    p[[index]]=ggplotGrob(plt_summ)

    ggsave(paste("results/machine learning/plots/", title,"/pdp-",str_replace(variable,"/","_"), ".png", sep=""), plt_summ, width=6, height=5, dpi=300)
    index=index+1
  }
  fig=grid.arrange(grobs=p, ncol=index-1, top=title)
  return(fig)
}
pdp.fig=list()
pdp.fig[[1]]=pdp_summ(readRDS("results/machine learning/models/all.rds"), "All", scale_info[["All"]])
pdp.fig[[2]]=pdp_summ(readRDS("results/machine learning/models/female.rds"), "Female", scale_info[["Female"]])
pdp.fig[[3]]=pdp_summ(readRDS("results/machine learning/models/male.rds"), "Male", scale_info[["Male"]])
ggsave(file="results/machine learning/plots/All-pdp.png",pdp.fig[[1]], width=200, height=6, dpi=150, units = "in", limitsize = FALSE)
ggsave(file="results/machine learning/plots/Female-pdp.png",pdp.fig[[2]], width=200, height=6, dpi=150, units = "in", limitsize = FALSE)
ggsave(file="results/machine learning/plots/Male-pdp.png",pdp.fig[[3]], width=200, height=6, dpi=150, units = "in", limitsize = FALSE)
ggsave(file="results/machine learning/plots/pdp.png",arrangeGrob(grobs=pdp.fig, ncol=1), width=200, height=18, dpi=150, units = "in", limitsize = FALSE)

Var.imp=function(models, title){
  temp=c()
  for (i in 1:length(models)){
    temp=c(temp,as.vector(varImp(models[[i]]$ranger)$importance))
  }
  summ=as.data.frame(do.call(cbind, temp))
  rownames(summ)=models[[1]]$ranger$finalModel$xNames
  summ.mean=apply(summ,1, mean, na.rm = TRUE)
  summ.sd=apply(summ,1, sd, na.rm = TRUE)
  summ.name <- factor(rownames(summ), levels = rownames(summ))
  summ.name_1=rev(levels(summ.name))
  summ.name_1[summ.name_1=="Income"]="lg Income"
  summ.name_1[summ.name_1=="Alcohol Use"]="lg Alcohol Use"
  summ.name_1[summ.name_1=="Caffeine Use"]="lg Caffeine Use"
  p<- ggplot() +
    geom_bar(mapping=aes(x=summ.name, y=summ.mean), stat="identity", color="black", position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=summ.mean-1.96*summ.sd, ymax=summ.mean+1.96*summ.sd, x=summ.name), width=.2,position=position_dodge(.9))+
    xlab("Variable")+ylab("Importance")+ggtitle(title)+scale_x_discrete(limits = rev(levels(summ.name)), labels=summ.name_1)+coord_flip()+
    theme(axis.text.x=element_text(color = "black", size=25, angle=0),axis.text.y=element_text(color = "black", size=25, angle=0),axis.title=element_text(size=27,face="bold"),plot.title = element_text(size = 29, face = "bold"))
  return(p)
}
rf.Fig=list()
rf.Fig[[1]]=print(Var.imp(readRDS("results/machine learning/models/all.rds"),"All"))
rf.Fig[[2]]=print(Var.imp(readRDS("results/machine learning/models/female.rds"),"Female"))
rf.Fig[[3]]=print(Var.imp(readRDS("results/machine learning/models/male.rds"),"Male"))
ggsave(file="results/machine learning/plots/All-rf.png",rf.Fig[[1]], width=15, height=20, units = "in", dpi=300)
ggsave(file="results/machine learning/plots/Female-rf.png",rf.Fig[[2]], width=15, height=20, units = "in", dpi=300)
ggsave(file="results/machine learning/plots/Male-rf.png",rf.Fig[[3]], width=15, height=20, units = "in", dpi=300)
ggsave(file="results/machine learning/plots/rf.png",arrangeGrob(grobs=rf.Fig, ncol=3), width=45, height=20, units = "in", dpi=300)



###
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
data=read.csv("data.csv", check.names=F, row.names=1, stringsAsFactors=TRUE)

p1=ggplot(data, aes(x=10^data$'C-reactive Protein'))+geom_histogram(bins=50)+xlab("Raw C-reactive Protein")
p2=ggplot(data, aes(x=data$'C-reactive Protein'))+geom_histogram(bins=50)+xlab("10-based Log-transformed C-reactive Protein")
ggsave(file="results/crp_hist.png", ggarrange(p1, p2, ncol = 2, nrow = 1), width=8, height=4, units = "in", dpi=300)

load("CRP_data.RData")
# data=na.omit(data)
p1=ggplot(data, aes(x=.data$'Body Mass Index', color=Gender, fill=Gender)) +
  geom_histogram(position="identity", alpha=0.5, bins = 50)+xlab('Body Mass Index')+theme(legend.position = "none")
p2=ggplot(data, aes(x=.data$'Percent Body Fat', color=Gender, fill=Gender)) +
  geom_histogram(position="identity", alpha=0.5, bins = 50)+xlab('Percent Body Fat')+theme(legend.position = "none")
p3=ggplot(data, aes(x=.data$'Waist to hips Ratio', color=Gender, fill=Gender)) +
  geom_histogram(position="identity", alpha=0.5, bins = 50)+xlab('Waist to hip Ratio')
ggsave(file="results/obesity_hist.png", ggarrange(p1, p2, p3, ncol = 3, nrow = 1,  widths=c(2.65, 2.65, 3.7)), width=9, height=3, units = "in", dpi=300)
