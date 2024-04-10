cd  /public1/guop/mawx/workspace/186sample/fastastructure
mkdir fastastructure2_20_result
seq 2 20 | parallel -j 20 "structure.py -K {} --input=/public1/guop/mawx/workspace/186sample/186_filtered.LD.pruned.noContig --output=/public1/guop/mawx/workspace/186sample/fastastructure/fastastructure2_20_result/LD_faststructure_K --cv=5 --prior=simple --seed=123 > LD_faststructure_K_{}.log 2>&1" &

## input为前缀  上一步output的后缀
chooseK.py --input=LD_faststructure_K


Model complexity that maximizes marginal likelihood = 2
Model components used to explain structure in data = 5

## 提取cv error
for f in *.log; do
    k=$(echo "$f" | sed -n 's/.*_K\.\([0-9]*\)\.log/\1/p')
    cv_error=$(grep "CV error" "$f")
    echo "$f (K=$k): $cv_error"
done

LD_faststructure_K.10.log (K=10): CV error = 0.4088825, 0.0119023
LD_faststructure_K.11.log (K=11): CV error = 0.4215564, 0.0191627
LD_faststructure_K.12.log (K=12): CV error = 0.4122411, 0.0263402
LD_faststructure_K.13.log (K=13): CV error = 0.4170378, 0.0102155
LD_faststructure_K.14.log (K=14): CV error = 0.4250738, 0.0264135
LD_faststructure_K.15.log (K=15): CV error = 0.4158595, 0.0171911
LD_faststructure_K.16.log (K=16): CV error = 0.4230420, 0.0153854
LD_faststructure_K.17.log (K=17): CV error = 0.4236116, 0.0228504
LD_faststructure_K.18.log (K=18): CV error = 0.4168027, 0.0190206
LD_faststructure_K.19.log (K=19): CV error = 0.4151018, 0.0158840
LD_faststructure_K.20.log (K=20): CV error = 0.4187051, 0.0177595
LD_faststructure_K.2.log (K=2): CV error = 0.4372752, 0.0154550
LD_faststructure_K.3.log (K=3): CV error = 0.4232509, 0.0167209
LD_faststructure_K.4.log (K=4): CV error = 0.4318602, 0.0117640
LD_faststructure_K.5.log (K=5): CV error = 0.4235034, 0.0094378
LD_faststructure_K.6.log (K=6): CV error = 0.4188243, 0.0130600
LD_faststructure_K.7.log (K=7): CV error = 0.4205870, 0.0108605
LD_faststructure_K.8.log (K=8): CV error = 0.4028336, 0.0179259
LD_faststructure_K.9.log (K=9): CV error = 0.4201673, 0.0096367

## 提取LLBO
# 创建一个空文件来保存结果
echo "File,Initialization,MarginalLikelihood" > marginal_likelihoods.csv

# 遍历所有.log文件
for f in *.log; do
    # 对每个文件，提取并处理边缘似然值
    grep "Marginal likelihood with initialization" "$f" | 
    awk -F' = ' -v fname="$f" '{print fname","substr($1, length($1)-1, 1)","$2}' >> marginal_likelihoods.csv
done

# 输出结果文件的路径
echo "Results saved to marginal_likelihoods.csv"

########################################################################  本地Rstudio中绘图   ##################################################################################
# 安装和加载必要的包
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# 读取数据
data <- read.csv("marginal_likelihoods.csv")

# 计算每个初始化的平均边缘似然值和标准差
aggregate_data <- aggregate(MarginalLikelihood ~ Initialization, data, function(x) c(mean = mean(as.numeric(x)), sd = sd(as.numeric(x))))
aggregate_data$Mean <- sapply(aggregate_data$MarginalLikelihood, `[`, 1)
aggregate_data$SD <- sapply(aggregate_data$MarginalLikelihood, `[`, 2)

# 绘制折线图和误差条
ggplot(aggregate_data, aes(x = Initialization, y = Mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = .1) +
  xlab("Initialization") +
  ylab("Marginal Likelihood") +
  ggtitle("Marginal Likelihood by Initialization") +
  theme_minimal()













