#############################  CV error  #######################################
# K值
k <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)

# CV错误值
cv_errors <- c(0.4372752, 0.4232509, 0.4318602, 0.4235034, 0.4188243, 0.4205870,
               0.4028336, 0.4201673, 0.4088825, 0.4215564, 0.4122411, 0.4170378, 
               0.4250738, 0.4158595, 0.4230420, 0.4236116, 0.4168027, 0.4151018,
               0.4187051)

# 标准差
std_devs <- c(0.0154550, 0.0167209, 0.0117640, 0.0094378, 0.0130600, 0.0108605,
              0.0179259, 0.0096367, 0.0119023, 0.0191627, 0.0263402, 0.0102155, 
              0.0264135, 0.0171911, 0.0153854, 0.0228504, 0.0190206, 0.0158840,
              0.0177595)

# 计算y轴范围
y_min <- min(cv_errors - std_devs) - 0.05  # 留出一点额外空间
y_max <- max(cv_errors + std_devs) + 0.05

# 绘制折线图并设置y轴范围
plot(k, cv_errors, type = "o", pch = 19, col = "blue", xlab = "K Value", 
     ylab = "CV Error", main = "CV Error vs K Value", ylim = c(y_min, y_max))

# 添加误差条
arrows(k, cv_errors - std_devs, k, cv_errors + std_devs, angle = 90, code = 3, 
       length = 0.03, col = "red")

library(ggplot2)

# 创建数据框
cv_data <- data.frame(K = k, CV_Error = cv_errors, SD = std_devs)

# 绘制折线图，添加误差条
ggplot(cv_data, aes(x = K, y = CV_Error)) +
  geom_line(colour = "blue") +  # 绘制蓝色折线
  geom_point(colour = "blue") +  # 添加蓝色点
  geom_errorbar(aes(ymin = CV_Error - SD, ymax = CV_Error + SD), width = 0.2, colour = "red") +  # 添加红色误差条
  labs(x = "K Value", y = "CV Error", title = "CV Error across different K values") +
  scale_x_continuous(breaks = k) + 
  scale_y_continuous(limits = c(0.37, 0.47)) +  # 设置y轴范围
  theme_bw()  




###############################  llbo #########################################
library(ggplot2)
library(dplyr)


# 读取数据
data <- read.csv("marginal_likelihoods.csv")

# 从文件名提取K值
data$K <- as.numeric(sub(".*_K\\.(\\d+)\\.log", "\\1", data$File))

# 计算每个K值的边缘似然值的平均值和标准差
stats <- data %>%
  group_by(K) %>%
  summarise(Mean = mean(MarginalLikelihood), SD = sd(MarginalLikelihood))

# 使用ggplot2绘图
library(ggplot2)

ggplot(stats, aes(x = K, y = Mean)) +
  geom_line(colour = "blue") +  # 绘制折线
  geom_point(colour = "blue") +  # 添加点
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, colour = "red") +  # 添加误差条
  labs(x = "K value", y = "Log-marginal likelihood lower bound (LLBO)", title = "Marginal Likelihood across different K values") +
  scale_x_continuous(breaks = k) + 
  theme_bw() 
