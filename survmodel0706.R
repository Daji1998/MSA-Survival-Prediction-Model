library(ggpubr)
library(corrplot)
library()
library(caret)
library(CBCgrps)
library(tidyverse)
library(rms)
library(pROC)
library(mice)
library(zoo)
library(VIM)
library(timeROC)
library(survivalROC)
library(dplyr)
library(forestmodel)
library(survival)
library(survminer)
library(autoReg)
library(ggsci)
library(forestplot)
library(forestploter)
library(ggprism)
getwd()
data <- read.csv("all_data_raw.csv")
table(is.na(data))
head(data)
str(data)
# 基线表----
basetab <- twogrps(data,gvar = "statues")
write.csv(basetab$Table,"基线表补充.csv")
missing_data <- sapply(data, function(x) sum(is.na(x)))
print(missing_data)
##缺失值统计表----
na_summary <- data.frame(
  Variable = names(data),
  Missing_Count = sapply(data, function(x) sum(is.na(x))),
  Missing_Percent = round(sapply(data, function(x) mean(is.na(x)) * 100), 2)
)
# 查看结果
na_summary[na_summary$Missing_Count > 0, ]  # 只显示有缺失值的列

# 按 statues 分组，计算 H_Y 的中位数和 IQR

summary_stats <- data %>%
  group_by(statues) %>%
  summarise(
    Q1 = quantile(H_Y, 0.25, na.rm = TRUE),
    Median = median(H_Y, na.rm = TRUE),
    Q3 = quantile(H_Y, 0.75, na.rm = TRUE),
    n = sum(!is.na(H_Y))
  )

print(summary_stats)

# 进行 Mann-Whitney U 检验（Wilcoxon rank-sum test）
wilcox_test_result <- wilcox.test(H_Y ~ statues, data = data, exact = FALSE, na.action = na.omit)

print(wilcox_test_result)

# 基本生存分析----
## 生存时间和KM曲线----
fit <- survfit(Surv(time,event)~1,data)
fit
summary(fit)
summary(fit, times = c(3, 5, 6, 9,10))
ggsurvplot(fit,data = data, pval = T,conf.int = T,conf.int.style="ribbon",
           surv.median.line = "hv",conf.int.alpha=0.2,
           palette = "jce",risk.table = T,family="sans", censor.size = 1.5,        # 控制删失标记的大小（默认约 2.5）
           censor.shape = 3 )
surv_median(fit)
data$Subtype <- factor(data$Subtype, levels = c(0, 1), labels = c("MSA-C", "MSA-P"))
fit <- survfit(Surv(time,event)~Subtype,data)
fit
summary(fit)
ggsurvplot(fit,data = data, pval = T,conf.int = T,conf.int.style="ribbon",
           surv.median.line = "hv",conf.int.alpha=0.1,
           palette = "jco",risk.table = T,family="sans", censor.size = 1.5,        # 控制删失标记的大小（默认约 2.5）
           censor.shape = 3)
surv_median(fit)

#折线图
data <- read.csv("featureselect.csv")
str(data)
ggplot(data, aes(x = features, y = Cindex)) +
  geom_line(color = "#1f77b4", size = 1) +
  geom_point(color = "#1f77b4", size = 2) +
  scale_x_reverse(breaks = seq(35, 1, by = -2)) +  # 横轴从35到1
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "#ec0000", size = 0.8) +  # 添加横线
  labs(x = "Number of Features", y = "C-index", title = "C-index vs. Number of Features") +
  theme_prism(base_size = 14)
point11 <- subset(data, features == 11)
ggplot(data, aes(x = features, y = Cindex)) +
  geom_line(color = "#1f77b4", size = 1) +
  geom_point(color = "#1f77b4", size = 2) +
  scale_x_reverse(breaks = seq(35, 1, by = -2)) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "#ec0000", size = 0.8) +
  # 添加注释点和文字
  geom_point(data = point11, aes(x = features, y = Cindex), color = "#ec0000", size = 3) +
  geom_text(data = point11, aes(label = round(Cindex, 3)), 
            vjust = -1, hjust = 0.5, color = "#ec0000", size = 4) +
  labs(x = "Number of Features", y = "C-index", title = "C-index vs. Number of Features") +
  theme_prism(base_size = 14)
##单因素.多因素Cox回归分析----
str(data)
names(data)
data1 <- data[,-c(1,2,8:10,13:15,21)]
str(data1)
factvalue <- sapply(data1,function(x)length(unique(x)))==2
data1[,factvalue] <- lapply(data1[,factvalue],as.factor)
data1$event <- as.numeric(as.character(data1$event)) 
str(data1)
cox_fit <- coxph(Surv(time,event)~.,data1)
reuslt_cox <- autoReg(cox_fit,uni=TRUE,multi=TRUE,threshold=0.1)
reuslt_cox
write.csv(reuslt_cox,"单因素多因素Cox回归.csv",row.names = F)

##多因素cox回归森林图----
data1$DiagnosticCategory <- factor(data$DiagnosticCategory, levels = c(0, 1), 
                                   labels = c("Clinically established", "Clinically probable"))
data1$TongueTremor <- factor(data1$TongueTremor, levels = c(0, 1), labels = c("without", "with"))
data1$NeurogenicOH <- factor(data1$NeurogenicOH, levels = c(0, 1), labels = c("without", "with"))
data1$PathologicalSigns <- factor(data1$PathologicalSigns, levels = c(0, 1), labels = c("without", "with"))
data1$FrequentFalls <- factor(data1$FrequentFalls, levels = c(0, 1), labels = c("without", "with"))
data1$RBD <- factor(data1$RBD, levels = c(0, 1), labels = c("without", "with"))
cox_fit2 <- coxph(Surv(time, event) ~ 
                    Age + DiagnosticCategory+
                    TongueTremor+
                    NeurogenicOH + ResidualUrine+
                    PathologicalSigns + UMSARS_total+
                    FrequentFalls + AOscore+
                    RBD + 
                    NLR + HCY+
                    VB9,
                  data = data1)
summary(cox_fit2)

panels <- list(
  forest_panel(width = 0.2, display = ~variable),
  forest_panel(width = 0.3, display = ~estimate),
  forest_panel(width = 0.6, item = "forest"),  # 主图更宽
  forest_panel(width = 0.1, item = "vline", line_x = 1)  # 中心线
)

forest_model(cox_fit2,format_options = forest_model_format_options(
  colour = "#ec0000",
  shape = 16,
  text_size = 5,
  point_size = 2,
  banded = F
))
modelPlot(cox_fit2, uni = FALSE, show.ref = FALSE, xlim = c(0.4, 4)) 
ggforest(cox_fit2, data = data1, main = "Multivariable Cox Forest Plot", 
         cpositions = c(0.02, 0.22, 0.4), fontsize = 1, refLabel = "Reference", noDigits = 2)

#预后模型----
##插补和分割后数据统计----
data <- read.csv("all_data_final.csv")
head(data)
str(data)
basedata <- twogrps(data,gvar = "test")
write.csv(basedata$Table,"插补分割后基线表.csv")

#RSF模型构建----
##构建生存数据机器学习模型----
traindata <- read.csv("train_dataset.csv")
testdata <- read.csv("test_dataset.csv")
names(traindata)
traindata <- traindata%>%select(-c("GaitAtaxia","RestTremor","Diabetes","COPD","AidWalking","CerebellarDysarthria",
                          "RBD","OculomotorFeatures","Hypertension","LimbAtaxia","MMSE","DiagnosticCategory","Gender",
                          "UrinaryUrgeIncontinence","Constipation","Rigidity","UMSARS_total","VB12","TongueTremor","Cscore",
                          "BMI","Subtype","AG_Ratio","bradykinesia"))
testdata <- testdata%>%select(-c("GaitAtaxia","RestTremor","Diabetes","COPD","AidWalking","CerebellarDysarthria",
                                   "RBD","OculomotorFeatures","Hypertension","LimbAtaxia","MMSE","DiagnosticCategory","Gender",
                                   "UrinaryUrgeIncontinence","Constipation","Rigidity","UMSARS_total","VB12","TongueTremor","Cscore",
                                   "BMI","Subtype","AG_Ratio","bradykinesia"))
###RSF####
library(randomForestSRC)

# 训练 RSF 模型
RSF <- rfsrc(seed = 42,
  formula = Surv(time,event) ~ . -time -event,
  data = traindata,
  ntree = 1042,              # 等价于 n_estimators
  nodesize = 12,             # 等价于 min_samples_leaf
  nsplit = 26,               # 近似 min_samples_split（分裂时随机选择的分裂点个数）
  mtry = floor(0.2 * (ncol(traindata) - 5)),  # max_features，0.2×有效特征数
  sampsize = floor(0.56587 * nrow(traindata)),  # max_samples，采样样本数
  max.depth = 3,             # 限制最大树深度（randomForestSRC 默认没有，需要用 `max.nodes` 近似）
  nthreads = parallel::detectCores()         # 等价于 n_jobs=-1
)
# 打印模型结果
print(RSF)

train_plot_data_rsf <- data.frame(
  traindata[,c("time","event")],
  predicted=as.numeric(predict(RSF, traindata)$predicted))

test_plot_data_rsf <- data.frame(
  testdata[,c("time","event")],
  predicted=as.numeric(predict(RSF, testdata)$predicted))

##模型预测值----
#整理训练集预测值
trainplot_data=cbind(train_plot_data_rsf[,c(1,2)],RSF=train_plot_data_rsf$predicted)

#整理测试集预测值
testplot_data=cbind(test_plot_data_rsf[,c(1,2)],RSF=test_plot_data_rsf$predicted)

###训练集AUC----
RSF = coxph(Surv(time,event) ~ RSF,trainplot_data,x=TRUE,y=TRUE)   #拟合Cox回归模型
# 计算模型AUC
AUCplot_lc <- riskRegression::Score(list("RSF" = RSF),
                                    formula = Surv(time, event) ~ 1,
                                    data = trainplot_data,
                                    metrics = "AUC",
                                    null.model = FALSE,
                                    times = seq(0, 12, 2))
# 获取 AUC 数据
auc_lc <- plotAUC(AUCplot_lc)
print(auc_lc)
# 绘制训练集模型的 AUC 曲线
ggplot() +
  geom_line(data = auc_lc, aes(x = times, y = AUC, group = model, col = model), linewidth = 1) +
  geom_point(data = auc_lc, aes(x = times, y = AUC, group = model, col = model), shape = 16) +
  geom_hline(yintercept = 0.7, linetype = 2, size = 0.8, colour = "grey") +
  theme_classic() + 
  labs(title = "AUC Curve for Lasso_Cox", x = "Time", y = "AUC") +
  ylim(0.5, 1) +
  theme(legend.title = element_blank())
###测试集AUC----
testRSF = coxph(Surv(time,event) ~ RSF,trainplot_data,x=TRUE,y=TRUE)

testAUCplot =riskRegression::Score(list("RSF"=testRSF),
                                   formula=Surv(time,event)~1,
                                   data=testplot_data,
                                   metrics="AUC",
                                   null.model=F,
                                   times=seq(0,12,1))
testauc = plotAUC(testAUCplot)
ggplot()+geom_line(data=testauc, aes(times,AUC,group=model,col=model),linewidth=1)+
  geom_point(data=testauc, aes(times,AUC,group=model,col=model),shape = 16)+
  #geom_ribbon(data=auc, aes(times,ymin = lower, ymax = upper,fill=model),alpha = 0.1)+
  geom_hline(yintercept=0.7, linetype=2,size=0.8,colour = "grey")+theme_classic()+ 
  labs(title = " ", x="Times", y="AUC")+
  ylim(0.5,1)+
  theme(legend.title = element_blank())
##训练集 C-index----
###绘制训练集C-index----
library(tidyr)  #gather函数
library(rms)    #cph函数

C_index = pec::cindex(list("RSF"=RSF),
                      formula=Surv(time,event)~1,
                      data=trainplot_data,
                      eval.times=seq(0,12,2))
cindex=data.frame(times=seq(0,12,2),C_index$AppCindex)
cindex = gather(cindex,key=Model,value = Cindex,-c(times))

ggplot()+geom_line(data=cindex, aes(times,Cindex,group=Model,col=Model),linewidth=1)+
  geom_point(data=cindex, aes(times,Cindex,group=Model,col=Model),shape = 16)+
  #geom_ribbon(data=cindex, aes(times,ymin = lower, ymax = upper,fill=Model),alpha = 0.1)+
  geom_hline(yintercept=0.7, linetype=2,size=0.8,colour = "grey")+theme_classic()+ 
  labs(title = " ", x="Times", y="C-index")+
  ylim(0.5,1)+
  theme(legend.title = element_blank())
###绘制测试集C-index----
testC_index = pec::cindex(list("RSF"=RSF),
                          formula=Surv(time,event)~1,
                          data=testplot_data,
                          eval.times=seq(4,12,2))
testcindex=data.frame(times=seq(4,12,2),testC_index$AppCindex)
testcindex = gather(testcindex,key=Model,value = Cindex,-c(times))

ggplot()+geom_line(data=testcindex, aes(times,Cindex,group=Model,col=Model),linewidth=1)+
  geom_point(data=testcindex, aes(times,Cindex,group=Model,col=Model),shape = 16)+
  #geom_ribbon(data=cindex, aes(times,ymin = lower, ymax = upper,fill=Model),alpha = 0.1)+
  geom_hline(yintercept=0.7, linetype=2,size=0.8,colour = "grey")+theme_classic()+ 
  labs(title = " ", x="Times", y="C-index")+
  ylim(0.5,1)+
  theme_classic()+
  theme(legend.title = element_blank())



##训练集 校准曲线----
train_cal_data=trainplot_data
colnames(train_cal_data)
time=5

f = cph(Surv(time, event) ~ RSF, x=T, y=T, surv=T, data=train_cal_data, time.inc=time)
cal = calibrate(f, cmethod="KM", method="boot", u=time, m=(nrow(train_cal_data)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="RSF predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="black", sub=F)


library(rms)

# 初始化结果列表
cal_list <- list()
time_points <- c(3, 5, 7)

# 循环每个时间点做校准
for (t in time_points) {
  f = cph(Surv(time, event) ~ RSF, x=T, y=T, surv=T, data=train_cal_data, time.inc=t)
  cal <- calibrate(f, cmethod="KM", method="boot", u=t, m=nrow(train_cal_data)/3, B=1000)
  cal_list[[as.character(t)]] <- cal
}

# 绘图
plot(cal_list[["3"]], lwd=1, xlim=c(0,1), ylim=c(0,1),font.lab=2,font.axis=2,cex.lab=1.2,
     xlab="RSF predicted OS", ylab="Observed OS", col="#00468B")
plot(cal_list[["5"]], add=TRUE, lwd=1, col="#ED0000")
plot(cal_list[["7"]], add=TRUE, lwd=1, col="#42B540")
abline(0,1, lty=1, col="gray")

legend("bottomright", legend=c("3 years", "5 years", "7 years"),
       col=c("#00468B", "#ED0000", "#42B540"), lwd=2,text.font=2)

##测试集 校准曲线----
test_cal_data=testplot_data
colnames(test_cal_data)
time=5

f = cph(Surv(time, event) ~ RSF, x=T, y=T, surv=T, data=test_cal_data, time.inc=time)
cal = calibrate(f, cmethod="KM", method="boot", u=time, m=(nrow(test_cal_data)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="RSF predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="black", sub=F)


library(rms)

# 初始化结果列表
cal_list <- list()
time_points <- c(3, 5, 7)

# 循环每个时间点做校准
for (t in time_points) {
  f = cph(Surv(time, event) ~ RSF, x=T, y=T, surv=T, data=test_cal_data, time.inc=t)
  cal <- calibrate(f, cmethod="KM", method="boot", u=t, m=nrow(test_cal_data)/2, B=100)
  cal_list[[as.character(t)]] <- cal
}

# 绘图
plot(cal_list[["3"]], lwd=1, xlim=c(0,1), ylim=c(0,1),font.lab=2,font.axis=2,cex.lab=1.2,
     xlab="RSF predicted OS", ylab="Observed OS", col="#00468B")
plot(cal_list[["5"]], add=TRUE, lwd=1, col="#ED0000")
plot(cal_list[["7"]], add=TRUE, lwd=1, col="#42B540")
abline(0,1, lty=1, col="gray")

legend("bottomright", legend=c("3 years", "5 years", "7 years"),
       col=c("#00468B", "#ED0000", "#42B540"), lwd=2,text.font=2)


##训练集 DCA----
traindca_data=traindata
traindca_data$RSF = c(1-summary(survfit(RSF,newdata = trainplot_data),times=6)$surv)
traindca=dcurves::dca(Surv(time, event) ~ RSF, data = traindca_data,time = 6)
as_tibble(traindca) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +
  geom_line(lwd=1.2) +
  coord_cartesian(ylim = c(-0.05, 0.5)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  #theme_bw()+
  theme_prism()+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.75),
        legend.text = element_text(size = 16))

##测试集 DCA----
testdca_data=testdata
testdca_data$RSF = c(1-summary(survfit(RSF,newdata = testplot_data),times=5.5)$surv)
testdca=dcurves::dca(Surv(time, event) ~ RSF, data = testdca_data,time = 5.5)
as_tibble(testdca) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +
  geom_line(lwd=1.2) +
  coord_cartesian(ylim = c(-0.05, 0.5)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  theme(legend.position = "inside",  # 启用图内图例定位
        legend.position.inside = c(0.95, 0.95),  # 设置相对位置（右上角）
        legend.justification = c(1, 1)) + 
  #theme_bw()+
  theme_prism()+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.75),
        legend.text = element_text(size = 16))


##训练集ROC----
library(plotROC)
library(pROC)
library(ggplot2)
dead=subset(trainplot_data,trainplot_data$event==1)
survival=subset(trainplot_data,trainplot_data$event==0)
time=5   ##指定绘制ROC的时间节点
dead$Result=ifelse(dead$time>time,0,1)
survival$Result=0
ROC_data=rbind(dead,survival)

str(ROC_data)

# 构建 ROC、AUC 和 CI
RSFROC <- roc(response = ROC_data$Result, predictor = ROC_data$RSF)
RSFAUC=round(pROC::auc(RSFROC),3)
RSFCI <- ci.auc(RSFROC)

# 绘图
pROC::ggroc(list(RSF = RSFROC),
            size=1,
            legacy.axes = T
)+theme_bw()+
  #geom_smooth()+
  labs(title = 'ROC curve')+
  theme(
    plot.title = element_text(hjust = 0.5,size = 15,face="bold"),
    axis.text=element_text(size=12,face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size=9,face="bold"),
    legend.position=c(0.65,0.15),
    legend.background = element_blank(),
    axis.title.y = element_text(size=12,face="bold"),#element_blank()
    axis.title.x = element_text(size=12,face="bold"),
    panel.border = element_rect(color="black",size=1),
    panel.background = element_blank())+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),        # 绘制对角线
               colour='grey', 
               linetype = 'dotdash'
  )+
  scale_colour_discrete(
    breaks = c("RSF"),
    labels = c(
      paste0("RSF (AUC = ", sprintf("%0.3f", RSFAUC),
             ", 95% CI: ", sprintf("%0.3f", RSFCI[1]),
             " - ", sprintf("%0.3f", RSFCI[3]), ")")
    )
  )

str(ROC_data)


# 时间点列表
times <- c(3, 5, 7)

# 存储不同时间点下的ROC对象和AUC
roc_list <- list()
auc_labels <- c()

for (t in times) {
  dead$Result <- ifelse(dead$time > t, 0, 1)
  ROC_data <- rbind(dead, survival)
  roc_obj <- roc(response = ROC_data$Result, predictor = ROC_data$RSF)
  auc_val <- round(pROC::auc(roc_obj), 3)
  ci_val <- ci.auc(roc_obj)
  
  # 加入roc列表（命名以便绘图）
  roc_list[[paste0("RSF_", t, "yr")]] <- roc_obj
  
  # 构造label
  auc_labels <- c(auc_labels,
                  paste0("RSF ", t, " year (AUC = ", sprintf("%0.3f", auc_val),
                         ", 95% CI: ", sprintf("%0.3f", ci_val[1]),
                         " - ", sprintf("%0.3f", ci_val[3]), ")"))
}

# 绘图
pROC::ggroc(roc_list, size=1, legacy.axes = TRUE) +
  theme_bw() +
  labs(title = '') +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    legend.position = c(0.65, 0.15),
    legend.background = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 1),
    panel.background = element_blank()
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
               colour = 'grey',
               linetype = 'dotdash') +
  scale_colour_discrete(
    breaks = names(roc_list),
    labels = auc_labels
  )

##测试集ROC----
dead=subset(testplot_data,testplot_data$event==1)
survival=subset(testplot_data,testplot_data$event==0)
time=5   ##指定绘制ROC的时间节点
dead$Result=ifelse(dead$time>time,0,1)
survival$Result=0
ROC_data=rbind(dead,survival)

str(ROC_data)

# 构建 ROC、AUC 和 CI
RSFROC <- roc(response = ROC_data$Result, predictor = ROC_data$RSF)
RSFAUC=round(pROC::auc(RSFROC),3)
RSFCI <- ci.auc(RSFROC)

# 绘图
pROC::ggroc(list(RSF = RSFROC),
            size=1,
            legacy.axes = T
)+theme_bw()+
  #geom_smooth()+
  labs(title = 'ROC curve')+
  theme(
    plot.title = element_text(hjust = 0.5,size = 15,face="bold"),
    axis.text=element_text(size=12,face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size=9,face="bold"),
    legend.position=c(0.65,0.15),
    legend.background = element_blank(),
    axis.title.y = element_text(size=12,face="bold"),#element_blank()
    axis.title.x = element_text(size=12,face="bold"),
    panel.border = element_rect(color="black",size=1),
    panel.background = element_blank())+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),        # 绘制对角线
               colour='grey', 
               linetype = 'dotdash'
  )+
  scale_colour_discrete(
    breaks = c("RSF"),
    labels = c(
      paste0("RSF (AUC = ", sprintf("%0.3f", RSFAUC),
             ", 95% CI: ", sprintf("%0.3f", RSFCI[1]),
             " - ", sprintf("%0.3f", RSFCI[3]), ")")
    )
  )



# 时间点列表
times <- c(3, 5, 7)

# 存储不同时间点下的ROC对象和AUC
roc_list <- list()
auc_labels <- c()

for (t in times) {
  dead$Result <- ifelse(dead$time > t, 0, 1)
  ROC_data <- rbind(dead, survival)
  roc_obj <- roc(response = ROC_data$Result, predictor = ROC_data$RSF)
  auc_val <- round(pROC::auc(roc_obj), 3)
  ci_val <- ci.auc(roc_obj)
  
  # 加入roc列表（命名以便绘图）
  roc_list[[paste0("RSF_", t, "yr")]] <- roc_obj
  
  # 构造label
  auc_labels <- c(auc_labels,
                  paste0("RSF ", t, " yearr (AUC = ", sprintf("%0.3f", auc_val),
                         ", 95% CI: ", sprintf("%0.3f", ci_val[1]),
                         " - ", sprintf("%0.3f", ci_val[3]), ")"))
}

# 绘图
pROC::ggroc(roc_list, size=1, legacy.axes = TRUE) +
  theme_bw() +
  labs(title = '') +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    legend.position = c(0.65, 0.15),
    legend.background = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 1),
    panel.background = element_blank()
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
               colour = 'grey',
               linetype = 'dotdash') +
  scale_colour_discrete(
    breaks = names(roc_list),
    labels = auc_labels
  )

##训练集 K-M 曲线----

library(survminer)  #ggsurvplot函数
###二分组----
km_data <- trainplot_data
colnames(km_data)
RSF_Median <- median(km_data$RSF)
km_data$RSF <- ifelse(km_data$RSF>RSF_Median,"high","low")

fit <- survfit(Surv(time, event) ~ RSF, data = km_data)
plot <- ggsurvplot(fit, data = km_data,
                   title = "RSF",
                   surv.median.line = "hv", # 添加中位数生存时间线
                   legend.title = "Risk", # 设置图例标题
                   legend.labs = c("high", "low"), # 指定图例分组标签
                   legend=c(0.9, 0.9),
                   pval = TRUE, # 设置添加P值
                   pval.method = TRUE, #设置添加P值计算方法
                   conf.int = TRUE, # 设置添加置信区间
                   risk.table = F, # 设置添加风险因子表
                   tables.height = 0.2, # 设置风险表的高度
                   tables.theme = theme_cleantable(), # 设置风险表的主题
                   palette = c("#ED0000FF",  "#0099B4FF"), # 设置颜色画板
                   ggtheme = theme_bw() # Change ggplot2 theme
)
plot$plot <- plot$plot + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)
print(plot)

###三分组----
km_data <- trainplot_data
# 使用 tertile 分组：low（< 33%）, middle（33%-66%）, high（> 66%）
quantiles <- quantile(km_data$RSF, probs = c(0.33, 0.67))  # 得到33%和66%分位数

km_data$RSF_group <- cut(km_data$RSF,
                         breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
                         labels = c("low", "middle", "high"),
                         right = TRUE)

# 查看每组人数
table(km_data$RSF_group)

# 构建生存对象
fit <- survfit(Surv(time, event) ~ RSF_group, data = km_data)

# 绘图
plot <- ggsurvplot(
  fit, data = km_data,
  title = "RSF",
  surv.median.line = "hv",
  legend.title = "Risk",
  legend.labs = c("low", "middle", "high"),  # 注意顺序要与 cut 分组一致
  legend = c(0.85, 0.85),
  pval = TRUE,
  pval.method = TRUE,
  censor.size = 2,
  conf.int = TRUE,
  conf.int.alpha = 0.08,
  risk.table = FALSE,
  tables.height = 0.2,
  tables.theme = theme_cleantable(),
  palette = c("#1f77b4", "#ff7f0e", "#42B540"),  # Lancet风格调色
  ggtheme = theme_bw()
)
plot$plot <- plot$plot + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)
print(plot)

##测试集 K-M 曲线----

###二分组----
km_data <- testplot_data
colnames(km_data)
RSF_Median <- median(km_data$RSF)
km_data$RSF <- ifelse(km_data$RSF>RSF_Median,"high","low")

fit <- survfit(Surv(time, event) ~ RSF, data = km_data)
plot <- ggsurvplot(fit, data = km_data,
                   title = "RSF",
                   surv.median.line = "hv", # 添加中位数生存时间线
                   legend.title = "Risk", # 设置图例标题
                   legend.labs = c("high", "low"), # 指定图例分组标签
                   legend=c(0.9, 0.9),
                   pval = TRUE, # 设置添加P值
                   pval.method = TRUE, #设置添加P值计算方法
                   conf.int = TRUE, # 设置添加置信区间
                   risk.table = F, # 设置添加风险因子表
                   tables.height = 0.2, # 设置风险表的高度
                   tables.theme = theme_cleantable(), # 设置风险表的主题
                   palette = c("#ED0000FF",  "#0099B4FF"), # 设置颜色画板
                   ggtheme = theme_bw() # Change ggplot2 theme
)
plot$plot <- plot$plot + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)
print(plot)

###三分组----
km_data <- testplot_data
# 使用 tertile 分组：low（< 33%）, middle（33%-66%）, high（> 66%）
quantiles <- quantile(km_data$RSF, probs = c(0.33, 0.67))  # 得到33%和66%分位数

km_data$RSF_group <- cut(km_data$RSF,
                         breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
                         labels = c("low", "middle", "high"),
                         right = TRUE)

# 查看每组人数
table(km_data$RSF_group)

# 构建生存对象
fit <- survfit(Surv(time, event) ~ RSF_group, data = km_data)

# 绘图
plot <- ggsurvplot(
  fit, data = km_data,
  title = "RSF",
  surv.median.line = "hv",
  legend.title = "Risk",
  legend.labs = c("low", "middle", "high"),  # 注意顺序要与 cut 分组一致
  legend = c(0.85, 0.85),
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  censor.size = 2,
  risk.table = FALSE,
  tables.height = 0.2,
  conf.int.alpha = 0.08,
  tables.theme = theme_cleantable(),
  palette = c("#1f77b4", "#ff7f0e", "#42B540"),  # Lancet风格调色
  ggtheme = theme_bw()
)
plot$plot <- plot$plot + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

print(plot)


#风险表KM曲线
library(readxl)
library(scales) 
getwd()
data <- read.csv("riskscore.csv")
data <- read_excel("train_dataset2.xlsx", sheet = 3)
table(is.na(data))
head(data)
data$group <- as.factor(data$group)
str(data)

# 基本生存分析----

## 生存时间和KM曲线----
fit <- survfit(Surv(time,event)~1,data)
fit
summary(fit)
ggsurvplot(fit,data = data, pval = T,conf.int = T,conf.int.style="ribbon",
           surv.median.line = "hv",conf.int.alpha=0.2,
           palette = "jce",risk.table = T,family="sans", censor.size = 1.5,        # 控制删失标记的大小（默认约 2.5）
           censor.shape = 3 )
surv_median(fit)
fit <- survfit(Surv(time,event)~group,data)
fit
summary(fit)
ggsurvplot(fit,data = data, pval = T,conf.int = F,conf.int.style="ribbon",
           surv.median.line = "hv",conf.int.alpha=0.1,
            palette = c("#1f77b4", "#ff7f0e", "#42B540"),risk.table = F,family="sans", censor.size = 0,        # 控制删失标记的大小（默认约 2.5）
           censor.shape = 3)
surv_median(fit)


p <- ggsurvplot(fit, data = data,
                pval = TRUE,
                conf.int = FALSE,
                conf.int.style = "ribbon",
                surv.median.line = "none",
                conf.int.alpha = 0.1,
                palette = c( "#ff7f0e","#42B540","#1f77b4" ),
                risk.table = FALSE,
                family = "sans",
                censor.size = 0,
                censor.shape = 3)

# 修改主图的 Y 轴为 1 - S(t)，即累计死亡率
p$plot <- p$plot +
  scale_y_reverse(labels = percent_format(accuracy = 1)) +  # 百分比格式
  labs(y = "Cumulative Incidence",x = "Time/years")+
  theme(legend.position = "bottom",
    legend.title = element_text(size = 14),   # 图注标题字号
    legend.text = element_text(size = 14)     # 图注内容字号
  )

# 显示图
print(p)

#roc曲线----
times <- c(3, 5, 7)  # 3年、6年、9年时间点，单位和data$time一致

roc_obj <- timeROC(T = data$time,
                   delta = data$event,
                   marker = data$riskscore,
                   cause = 1,
                   weighting = "marginal",
                   times = times,
                   iid = TRUE)
# 基础绘图
plot(roc_obj, time = 3, col = "#1f77b4", lwd = 2, title = FALSE)
plot(roc_obj, time = 5, col =  "#ff7f0e", lwd = 2, add = TRUE)
plot(roc_obj, time = 7, col =  "#42B540", lwd = 2, add = TRUE)
legend("bottomright",
       legend = paste0(times, " years AUC: ", 
                       format(round(roc_obj$AUC, 3), nsmall = 3)),
       col = c("#1f77b4", "#ff7f0e", "#42B540"),
       lwd = 2)
title("Time-dependent ROC curves for riskscore")
# 自定义函数来绘制平滑曲线
plot_smooth_ROC <- function(roc_obj, time, col, lwd = 2, add = FALSE) {
  fpr <- roc_obj$FP[, which(roc_obj$times == time)]
  tpr <- roc_obj$TP[, which(roc_obj$times == time)]
  
  ss <- smooth.spline(fpr, tpr, spar = 0.5)  # spar 越大越平滑（0~1之间调节）
  if (!add) {
    plot(ss$x, ss$y, type = "l", col = col, lwd = lwd,
         xlab = "False positive rate", ylab = "True positive rate",
         xlim = c(0, 1), ylim = c(0, 1))
  } else {
    lines(ss$x, ss$y, col = col, lwd = lwd)
  }
}

# 示例绘图
times <- c(3, 5, 7)
plot_smooth_ROC(roc_obj, time = 3, col = "#1f77b4")
plot_smooth_ROC(roc_obj, time = 5, col = "#ff7f0e", add = TRUE)
plot_smooth_ROC(roc_obj, time = 7, col = "#42B540", add = TRUE)

legend("bottomright",
       legend = paste0(times, " years AUC: ", 
                       format(round(roc_obj$AUC, 3), nsmall = 3)),
       col = c("#1f77b4", "#ff7f0e", "#42B540"),
       lwd = 2)

title("Time-dependent ROC curves for riskscore (smoothed)")


#rcs分析
library(rms)
library(survival)
dd <- datadist(data)
options(datadist = "dd")
fit <- cph(Surv(time, event) ~ rcs(riskscore, 4), data = data, x = TRUE, y = TRUE, surv = TRUE)
plot(Predict(fit, riskscore, fun = exp),
     xlab = "Risk score",
     ylab = "Hazard Ratio (HR)",
     main = "Restricted Cubic Spline of Risk Score",
     lwd = 2,
     col = "blue",
     cex.lab = 1.2,
     cex.axis = 1.1)
grid.lines(x = unit(c(0, 1), "npc"),
           y = unit(c(1, 1), "native"),
           gp = gpar(col = "red", lty = 2))
rcspline.eval(data$riskscore, knots = 4, inclx = TRUE)  # 可查看自动生成的结点位置
fit <- cph(Surv(time, event) ~ rcs(riskscore,5), data = data, x = TRUE, y = TRUE, surv = TRUE)

# 获取预测值
pred <- Predict(fit, riskscore, fun = exp)

# 转为 data.frame
df_pred <- as.data.frame(pred)

# 使用 ggplot 绘制
ggplot(df_pred, aes(x = riskscore, y = yhat)) +
  geom_line(color = "#1f77b4", size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, fill = "#1f77b4") +
  geom_hline(yintercept = 1, lwd=1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(5, 10), linetype = "dashed", color = "darkgray") +
  labs(x = "Risk Score", y = "Hazard Ratio (HR)", 
       title = "Restricted Cubic Spline of Risk Score") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),  # 移除主网格线
        panel.grid.minor = element_blank())  # 移除次网格线

