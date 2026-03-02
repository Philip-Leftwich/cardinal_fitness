
#open data
data <- read.csv("data/Cardinal_Rps7_Results_20250611.csv")

#load ggplot2
library(ggplot2)

#remove NTCs
library(dplyr)
tidy_data <- data %>%
  filter(Sample != "NTC")

#convert Ct column into numerical
summarised_data <- tidy_data %>%
  mutate(Cq = as.numeric(Cq)) %>%
  group_by(Sample, Target, Replicate) %>%
  summarise(Cq.mean = mean(Cq))

#Plot the mean_Ct for a quick comparison
ggplot(summarised_data, aes(x = Sample, y = Cq.mean, colour = Target)) +
  geom_point()

##Analyse the data using the delta-delta-Ct method
target_data <- summarised_data %>%
  filter(Target == "Cardinal")

ref_data <- summarised_data %>%
  filter(Target == "rps7") %>%
  rename("Cq.ref" = "Cq.mean")
#use the left.join function to merge both tables according to Sample_Name
combined_data <- left_join(target_data, ref_data, by = c("Sample", "Replicate"))
##Use the mutate function to create a new column containing the delta Ct value
combined_data <- mutate(combined_data, delta_Ct = Cq.mean - Cq.ref)
#plot the deltaCt values
ggplot(combined_data, aes(x = Sample, y = delta_Ct)) +
  geom_point()

##calculate the mean_deltaCt value for each treatment
treatment_summary <- combined_data %>%
  group_by(Sample, Replicate) %>%
  summarise(mean_deltaCt = mean(delta_Ct))
##calculate the delta delta Ct value of each replicate compared to the control sample
mean_control <- filter(treatment_summary, Sample == "SDA-500 fem") %>%
  pull(mean_deltaCt)

combined_data <- combined_data %>%
  mutate(delta_delta_Ct = delta_Ct - mean_control)

ggplot(combined_data, aes(x = Sample, y = delta_delta_Ct)) +
  geom_point()

combined_data <- combined_data %>%
  mutate(Relative_Expression = 2^(-delta_delta_Ct))
                      
#plot the mean
library(ggpubr)
library(tidyr)

#separate sample in sample and sex
combined_data <- combined_data %>%
separate(Sample, into = c("Sample", "Sex"), sep = " ", remove = TRUE)

#Rename fem into female in the sex column
combined_data$Sex <- recode(combined_data$Sex, "fem" = "Female", "female" = "Female", "male" = "Male")

#write list of comparisons for the female sex
comparisons <- list(
  c("2360B5", "SDA-500"),
  c("1759", "SDA-500"),
  c("D251", "SDA-500"),
  c("QA383P", "SDA-500")
)

#change the order of the sample column
combined_data$Sample <- factor(combined_data$Sample, levels = c(
  "2360B5",
  "1759",
  "D251",
  "QA383P",
  "SDA-500"))

#check levels
levels(combined_data$Sample)

sample_name <- c(
  "2360B5" = expression("cd"^"g225"),
  "1759" = expression("cd"^"g384"),
  "D251" = expression("cd"^225),
  "QA383P" = expression("cd"^384),
  "SDA-500" = expression("WT")
)

#Define colors for each sample
id_colors <- c(
  "2360B5" = "#8BABD3",
  "1759" = "#BB8BD3",
  "D251" = "#FFA040",
  "QA383P" = "#FFE699",
  "SDA-500" = "lightgray")
                  
#Plot the graph
ggplot(combined_data, aes(x = Sample, y = Relative_Expression, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +          
  geom_jitter(width = 0.2, alpha = 0.6) +     
  facet_grid(. ~Sex, scales = "free_x", space = "free_x") +
  stat_compare_means(comparisons = comparisons,
                     method = "t.test", label = "p.signif", size = 3) +
  scale_x_discrete(labels = sample_name) +
  labs(x = NULL, y = "Relative Expression") +
    theme_bw() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none") +
  scale_fill_manual(values = id_colors)






#open data
data2 <- read.csv("data/Cardinal_GADPH_Results_20250613.csv")


#remove NTC from Sample
tidy_data2 <- data2 %>%
  filter(Sample !="NTC")

#convert Cq to numeric and group group table by sample, target and replicate
summarised_data2 <- tidy_data2 %>%
  mutate(Cq = as.numeric(Cq)) %>%
  group_by(Sample, Target, Replicate) %>%
  summarise(Cq.mean = mean(Cq))

#filter by target
target_data2 <- summarised_data2 %>%
  filter(Target == "Cardinal")

#filter by reference gene
ref_data2 <- summarised_data2 %>%
  filter(Target == "GAPDH") %>%
  rename("Cq.ref" = "Cq.mean")

#combine both tables
combined_data2 <- left_join(target_data2, ref_data2, by = c("Sample", "Replicate"))
# calculate delta_Ct
combined_data2 <- combined_data2 %>%
  mutate(delta_CT = Cq.mean - Cq.ref)
#calculate mean delta_Ct value
treatment_summary2 <- combined_data2 %>%
  group_by(Sample, Replicate) %>%
  summarise(delta_CT.mean = mean(delta_CT))
#Indicate reference sample
mean_control2 <- filter(treatment_summary2, Sample == "SDA-500 fem") %>%
  pull(delta_CT.mean)
#calculate delta_delta_Ct value
combined_data2 <- combined_data2 %>%
  mutate(delta_delta_Ct = delta_CT - mean_control2)
#calculate relative expression
combined_data2 <- combined_data2 %>%
  mutate(Relative_Expression = 2^(-delta_delta_Ct))

#separate sample from sex
combined_data2 <- combined_data2 %>%
  separate(Sample, into = c("Sample", "Sex"), sep = " ")

#Rename fem, female and male into Female and Male
combined_data2$Sex <- recode(combined_data2$Sex, "fem" = "Female", "female" ="Female", "male" = "Male")

#Change the order of the Sample column
combined_data2$Sample <- factor(combined_data2$Sample, levels = c(
  "2360B5",
  "1759",
  "D251",
  "QA383P",
  "SDA-500"
))

#check the order
levels(combined_data2$Sample)

#Change the Sample names
sample_name <- c(
  "2360B5" = expression("cd"^"g225"),
  "1759" = expression("cd"^"g384"),
  "D251" = expression("cd"^225),
  "QA383P" = expression("cd"^384),
  "SDA-500" = expression("WT")
)

#set up comparisons
comparisons <- list(
  c("2360B5", "SDA-500"),
  c("1759", "SDA-500"),
  c("D251", "SDA-500"),
  c("QA383P", "SDA-500")
)

#Define colors for each sample
id_colors <- c(
  "2360B5" = "#8BABD3",
  "1759" = "#BB8BD3",
  "D251" = "#FFA040",
  "QA383P" = "#FFE699",
  "SDA-500" = "lightgray")

#plot graph
ggplot(combined_data2, aes(x = Sample, y = Relative_Expression, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_grid(.~Sex, scales = "free_x", space = "free_x") +
  stat_compare_means(comparisons = comparisons,
                     method = "t.test", label = "p.signif", size = 3) +
  scale_x_discrete(labels = sample_name) +
  labs(x = NULL, y = "Relative Expression") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none") +
  scale_fill_manual(values = id_colors)


##
model <- glm(delta_delta_Ct ~ Sample * Sex, data = combined_data2, family = gaussian)

emm <- emmeans(model, ~ Sample | Sex)

# Step 2: Apply treatment vs control contrasts within each Sex
# Replace "Control" with the actual control Sample level
contrasts_by_sex <- contrast(emm, method = "trt.vs.ctrl", ref = "SDA-500", adjust = "none")

# Step 3: View the results
summary(contrasts_by_sex)
