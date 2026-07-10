# 2026-07-10 AndyP


# 1. Load required library for multivariate normal simulation
# install.packages("MASS") 
library(MASS)

# 2. Set seed for reproducible results
set.seed(123)

# 3. Define simulation parameters
total_points <- 500
correlation <- 0.6  

# Define the covariance matrix
cov_matrix <- matrix(c(1, correlation, correlation, 1), nrow = 2)

# 4. Generate the base correlated data
raw_data <- mvrnorm(n = total_points, mu = c(0, 0), Sigma = cov_matrix)

# 5. Create a data frame
my_data <- data.frame(
  Var1 = raw_data[, 1],
  Var2 = raw_data[, 2]
)

# 6. Assign group labels matching the 75% vs 25% ratio
my_data$Group <- sample(
  c("Group_A", "Group_B"), 
  size = total_points, 
  replace = TRUE, 
  prob = c(0.75, 0.25)
)

# 7. Add group differences
# Large effect for Var1, small effect for Var2
var1_effect <- 1.5
var2_effect <- 0.4 

my_data$Var1 <- ifelse(my_data$Group == "Group_B", my_data$Var1 + var1_effect, my_data$Var1)
my_data$Var2 <- ifelse(my_data$Group == "Group_B", my_data$Var2 + var2_effect, my_data$Var2)

# 8. Verification checks
cat("--- Group Counts ---\n")
print(table(my_data$Group))

cat("\n--- Overall Correlation Matrix ---\n")
print(cor(my_data[, c("Var1", "Var2")]))

cat("\n--- Mean of Var1 by Group (Large difference) ---\n")
print(aggregate(Var1 ~ Group, data = my_data, FUN = mean))

cat("\n--- Mean of Var2 by Group (Small difference) ---\n")
print(aggregate(Var2 ~ Group, data = my_data, FUN = mean))


df <- my_data



#########

# 2. Set seed for reproducible results
set.seed(123)

# 3. Define simulation parameters
total_points <- 500
correlation <- 0.6  

# Define the covariance matrix
cov_matrix <- matrix(c(1, correlation, correlation, 1), nrow = 2)

# 4. Generate the base correlated data
raw_data <- mvrnorm(n = total_points, mu = c(0, 0), Sigma = cov_matrix)

# 5. Create a data frame
my_data <- data.frame(
  Var1 = raw_data[, 1],
  Var2 = raw_data[, 2]
)

# 6. Assign group labels matching the 75% vs 25% ratio
my_data$Group <- sample(
  c("Group_A", "Group_B"), 
  size = total_points, 
  replace = TRUE, 
  prob = c(0.75, 0.25)
)

# 7. Add group differences
# Large positive effect for Var1, small negative effect for Var2
var1_effect <- 1.5
var2_effect <- -0.4 

my_data$Var1 <- ifelse(my_data$Group == "Group_B", my_data$Var1 + var1_effect, my_data$Var1)
my_data$Var2 <- ifelse(my_data$Group == "Group_B", my_data$Var2 + var2_effect, my_data$Var2)

# 8. Verification checks
cat("--- Group Counts ---\n")
print(table(my_data$Group))

cat("\n--- Overall Correlation Matrix ---\n")
print(cor(my_data[, c("Var1", "Var2")]))

cat("\n--- Mean of Var1 by Group (Positive difference) ---\n")
print(aggregate(Var1 ~ Group, data = my_data, FUN = mean))

cat("\n--- Mean of Var2 by Group (Negative difference) ---\n")
print(aggregate(Var2 ~ Group, data = my_data, FUN = mean))

df1 <- my_data



library(MASS)
set.seed(123)

total_points <- 500
correlation <- 0.6  
cov_matrix <- matrix(c(1, correlation, correlation, 1), nrow = 2)
raw_data <- mvrnorm(n = total_points, mu = c(0, 0), Sigma = cov_matrix)

# 1. SHIFT BASELINE MEANS TO 100
var1_base <- raw_data[, 1] + 100
var2_base <- raw_data[, 2] + 100

group_labels <- sample(c("Group_A", "Group_B"), size = total_points, replace = TRUE, prob = c(0.75, 0.25))

# Create df
df <- data.frame(Var1 = var1_base, Var2 = var2_base, Group = group_labels)
df$Var1 <- ifelse(df$Group == "Group_B", df$Var1 + 1.5, df$Var1)
df$Var2 <- ifelse(df$Group == "Group_B", df$Var2 + 0.4, df$Var2)
df$Ratio <- df$Var2 / df$Var1

# Create df1
df1 <- data.frame(Var1 = var1_base, Var2 = var2_base, Group = group_labels)
df1$Var1 <- ifelse(df1$Group == "Group_B", df1$Var1 + 1.5, df1$Var1)
df1$Var2 <- ifelse(df1$Group == "Group_B", df1$Var2 - 0.4, df1$Var2)
df1$Ratio <- df1$Var2 / df1$Var1

# --- TEST 1: RAW RATIO MODELS (OLS) ---
cat("\n=== OLS RATIO MODEL FOR DF (Both Up) ===\n")
print(summary(lm(Ratio ~ Group, data = df))$coefficients)

cat("\n=== OLS RATIO MODEL FOR DF1 (Opposite) ===\n")
print(summary(lm(Ratio ~ Group, data = df1))$coefficients)


# --- TEST 2: GLM RATIO MODELS (Gamma with Log Link) ---
# We use 'offset(log(Var1))' to force the ratio framework mathematically
glm_df  <- glm(Var2 ~ Group + offset(log(Var1)), family = Gamma(link = "log"), data = df)
glm_df1 <- glm(Var2 ~ Group + offset(log(Var1)), family = Gamma(link = "log"), data = df1)

cat("\n=== GLM RATIO MODEL FOR DF ===\n")
print(summary(glm_df)$coefficients)

cat("\n=== GLM RATIO MODEL FOR DF1 ===\n")
print(summary(glm_df1)$coefficients)


print(df1)