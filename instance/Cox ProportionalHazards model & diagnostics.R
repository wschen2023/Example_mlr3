## Cox Proportional Hazards Model & Diagnostics
## 2024-05-15
## R 4.3.1
## mlr3verse 0.2.8

## Clean up the working environment
rm(list=ls())
## Set seed
set.seed(2024)


## Load package
library(mlr3verse)
# install.packages("mlr3proba", repos = "https://mlr-org.r-universe.dev")
# install.packages("ooplah")
# install.packages("dictionar6")
library(mlr3proba)
library(survival)

## Import Data
## Task Definition
task = tsk("gbcs")
task
task$head()

# time：观察时间或生存时间 
# status：事件状态
# age：患者的年龄
# estrg_recp：雌激素受体状态 
# grade：肿瘤分级 
# hormone：激素治疗情况 
# menopause：更年期状态 
# nodes：正性淋巴结的数量 
# prog_recp：孕激素受体状态 
# size：肿瘤大小（以毫米为单位）

## Kaplan-Meier plot
# install.packages("GGally")
autoplot(task, type = "target")
# Group by "grade"
autoplot(task, rhs = "grade")
# Group by "hormone"
autoplot(task, rhs = "hormone")

## Learner Definition
learner_cox <- lrn("surv.coxph")

## Feature Selection
afs = auto_fselector(
  fselector = fs("exhaustive_search"), #详尽搜索，尝试所有可能的特征子集
  learner = learner_cox,
  resampling = rsmp("cv", folds = 3),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 100)
)

future::plan("multisession")

rr <- resample(task, afs, rsmp("cv", folds = 4), store_models = TRUE)
k <- c("surv.rcll", "surv.cindex", "surv.dcalib")
# 平均模型性能
rr$aggregate(measures = msrs(k)) 
# 外层4次迭代的结果
rr$score(measures = msrs(k)) 
# 内层每次调参的结果
as.data.table(extract_inner_fselect_results(rr))[, 1:10]
# age, menopause在4次的特征选择中都没有采用
# nodes, prog_recp在4次的特征选择中都被采用

## Performance Evaluation
# 与surv.kaplan模型（基线模型）进行比较
learner_km <- lrns("surv.kaplan")
learners = c(afs, learner_km) 
future::plan("multisession")
bmr = benchmark(benchmark_grid(task, learners,
                               rsmp("cv", folds = 3)))
bmr$aggregate(measures = msrs(k))[, c("learner_id", 
                                      "surv.rcll", 
                                      "surv.cindex", 
                                      "surv.dcalib")]

## Model Feature Extraction
model <- afs$train(task)$model
model$features

# INFO  [16:48:57.753] [bbotk] Result:
# INFO  [16:48:57.756] [bbotk]     age estrg_recp  grade hormone menopause  nodes prog_recp   size        features n_features surv.cindex
# INFO  [16:48:57.756] [bbotk]  <lgcl>     <lgcl> <lgcl>  <lgcl>    <lgcl> <lgcl>    <lgcl> <lgcl>          <list>      <int>       <num>
# INFO  [16:48:57.756] [bbotk]   FALSE      FALSE  FALSE   FALSE     FALSE   TRUE      TRUE  FALSE nodes,prog_recp          2   0.7457936

## Training
## Model training using selected features
task_cox <- task$clone()$select(model$features)
model_cox <- learner_cox$train(task_cox)$model
summary(model_cox)

## Testing proportional Hazards assumption
test_ph <- cox.zph(model_cox)
test_ph

# chisq df       p
# nodes      1.05  1 0.30628
# prog_recp 12.68  1 0.00037
# GLOBAL    14.44  2 0.00073

# prog_recp变量的检验和全局检验的p值<0.05,说明数据违背了比例风险的假定

op <- par(mfrow = c(1, 2))
plot(test_ph)

## Schoenfeld residuals 
#  如果Cox模型的假设成立，那么残差应该随着时间或事件的发生是随机分布的，不应该有明显的趋势或者模式
library(survminer)
ggcoxzph(test_ph)
ggcoxdiagnostics(model_cox, type = "schoenfeld")
# 观察prog_recp变量的数值
temp <- task_cox$data()
op <- par(mfrow = c(1, 1))
hist(temp$prog_recp, breaks = 100)

# prog_recp数据是偏态分布
# 尝试分箱
library(mlr3pipelines)
po <- po("quantilebin", numsplits = 4, affect_columns = selector_name("prog_recp"))
task_cox_2 <- po$train(list(task_cox))[[1]]

# 按prog_recp特征的分组生存曲线
autoplot(task_cox_2, rhs = "prog_recp")

# KM曲线存在交叉，意味prog_recp特征可能还是不符合比你风险假设
# 把prog_recp特征，进行分位数分箱分成4分
# 进行独热编码，最后连接上学习器
po_encode <- po("encode", affect_columns = selector_name("prog_recp"))
gra_cox_bin <- po %>>% 
  po_encode %>>%
  learner_cox
# install.packages("igraph")
library(igraph)
op <- par(mfrow = c(1, 1))
gra_cox_bin$plot()

## Performance Evaluation
learner_cox_bin <- as_learner(gra_cox_bin)

# surv.kaplan模型（基线模型）
gra_km_bin <- po %>>% 
  po_encode %>>%
  learner_km
gra_km_bin$plot()

learner_km_bin <- as_learner(gra_km_bin)

learners = c(learner_cox_bin, learner_km_bin)

bmr = benchmark(benchmark_grid(task_cox, learners,
                               rsmp("cv", folds = 4)))
bmr$aggregate(measures = msrs(k))[, c("learner_id", "surv.cindex", "surv.dcalib", "surv.rcll" )]

## 检验分箱后的模型是否符合cox回归的比例风险假设
model_cox_bin <- learner_cox_bin$train(task_cox)$model$surv.coxph$model
summary(model_cox_bin)
test_hp_bin <- cox.zph(model_cox_bin)
test_hp_bin
# 分箱并没有改变不符合cox回归的比例风险假设的现象
plot(test_hp_bin)
ggcoxdiagnostics(model_cox_bin, type = "schoenfeld")

## 解决COX回归不符合比例风险假设
## 1）对时间分层: 对时间使用分层函数, 根据拐点把时间分成几段
## 2）连续性时依系数变换：prog_recp系数随时间变化的曲线明显不是线性的，可以通过数据变换把它变成类似线性的，比如取log
## 3) 样条技术

##Import Data
library(mlr3proba)
task = tsk("gbcs")
data <- task$data()

fit_cox <- coxph(Surv(time, status) ~ nodes + prog_recp, 
                 data = data)
fit_cox
test_ph <- cox.zph(fit_cox)
test_ph

op <- par(mfrow = c(1, 1))
plot(test_ph[2])
# prog_recp变量的系数随着时间的改变，prog_recp偏离的比较厉害, 系数大概从810天开始增加, 趋近于0；
# 1900天又开始下降, 所以prog_recp的系数是一直在随着时间改变的, 不符合比例风险假设

## 1)对时间分层
## 对时间使用分层函数
stra <- survSplit(Surv(time, status) ~ ., data = data, 
                  cut=c(810, 1900), # 两个拐点把时间分为3层
                  episode= "tgroup", 
                  id = "id")
fit_cox_stra <- coxph(Surv(tstart, time, status) ~ nodes + prog_recp:strata(tgroup), 
                      data = stra)
fit_cox_stra
test_ph <- cox.zph(fit_cox_stra)
test_ph

# K-M plot
# 按prog_recp特征的分组生存曲线
task_stra = TaskSurv$new("gbcs", stra, event = "status")
autoplot(task_stra, rhs = "tgroup")

## 连续性时依系数变换
# prog_recp系数随时间变化的曲线明显不是线性的
# 我们可以通过数据变换把它变成类似线性的，比如取log
# 可通过tt(time transform)函数实现
# 构建时依协变量时，可以选择x * t、x * log(t)、x * log(t + 20)、x * log(t + 200)等等
# 没有明确的规定，要结合结果和图示进行选择
fit_cox_log <- coxph(Surv(time, status) ~ nodes + prog_recp + tt(prog_recp), # 对prog_recp进行变换
              data = data, 
              tt = function(x, t, ...) x * log(t) 
              )
fit_cox_log
# Call:
#   coxph(formula = Surv(time, status) ~ nodes + prog_recp + tt(prog_recp), 
#         data = data, tt = function(x, t, ...) x * log(t))
# 
# coef exp(coef)  se(coef)      z        p
# nodes          0.064031  1.066125  0.008728  7.337 2.19e-13
# prog_recp     -0.076417  0.926430  0.019173 -3.986 6.73e-05
# tt(prog_recp)  0.009967  1.010017  0.002659  3.749 0.000178
# 
# Likelihood ratio test=110.5  on 3 df, p=< 2.2e-16
# n= 686, number of events= 171

# prog_recp的时依系数估计为：0.009967 * log(t)
# 变换后的PH检验
test_ph <- cox.zph(fit_cox, transform = function(time) log(time))
plot(test_ph[2])
abline(0, 0, col = "red") # 0水平线
abline(h = fit_cox$coef[2], col = "green", lwd = 2, lty = 2) # 整体估计
abline(coef(fit_cox_log)[2:3], col = "blue", lwd = 2, lty = 3) # 现在的估计

## 样条技术
# 使用样条技术pspline(time)进行数据变换
fit_cox_ps_1 <- coxph(Surv(time, status) ~ nodes + prog_recp + tt(prog_recp), # 对prog_recp进行变换
              data = data, 
              tt = function(x, t, ...) x * pspline(t) 
              )
fit_cox_ps_1

# Call:
#   coxph(formula = Surv(time, status) ~ nodes + prog_recp + tt(prog_recp), 
#         data = data, tt = function(x, t, ...) x * pspline(t))
# 
# coef  se(coef)       se2     Chisq   DF       p
# nodes                  6.48e-02  8.91e-03  8.91e-03  5.28e+01 1.00 3.6e-13
# prog_recp             -5.46e-02  3.74e-02  2.34e-02  2.13e+00 1.00  0.1441
# tt(prog_recp), linear  6.55e-06  3.06e-06  3.06e-06  4.56e+00 1.00  0.0326
# tt(prog_recp), nonlin                                1.45e+01 3.61  0.0042
# 
# Iterations: 10 outer, 40 Newton-Raphson
# Theta= 1 
# Degrees of freedom for terms= 1.0 0.4 4.6 
# Likelihood ratio test=124  on 6 df, p=<2e-16
# n= 686, number of events= 171 

# prog_recp在线性平滑函数下对风险没有显著影响
# 但在非线性平滑函数下具有显著影响

# 还可以使用pspline()函数在R代码里主要是修改了COX回归的基线风险函数
# 这样，基线风险函数不再是固定的，而是可以随时间变化
# 调整了时间的非线性变化
fit_cox_ps_2 <- coxph(Surv(time, status) ~ nodes + prog_recp + pspline(time),
                data = data
                )
fit_cox_ps_2

# Call:
#   coxph(formula = Surv(time, status) ~ nodes + prog_recp + pspline(time), 
#         data = data)
# 
# coef  se(coef)       se2     Chisq   DF       p
# nodes                  3.62e-02  1.06e-02  1.06e-02  1.17e+01 1.00 0.00063
# prog_recp             -3.75e-03  8.31e-04  8.31e-04  2.04e+01 1.00 6.3e-06
# pspline(time), linear -2.83e-02  1.31e-03  1.31e-03  4.65e+02 1.00 < 2e-16
# pspline(time), nonlin                                2.88e+02 2.99 < 2e-16
# 
# Iterations: 3 outer, 60 Newton-Raphson
# Theta= 0.273 
# Degrees of freedom for terms= 1 1 4 
# Likelihood ratio test=1108  on 5.99 df, p=<2e-16
# n= 686, number of events= 171 

# 时间的多项式样条变换（pspline(time)）在统计上是显著的
# 这表明生存风险是随时间的变化的
# 在调整了时间的非线性变化（即基线风险的波动）之后nodes和prog_recp的影响在统计上依然显著
test_ph <- cox.zph(fit_cox_ps_2)
test_ph
op <- par(mfrow = c(1, 1))
plot(test_ph)
