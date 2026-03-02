
path <- "clean_data/SurvivalData.xlsx"

data <- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path, sheet = .x), .id = "line")

data <- data %>% 
  mutate(Genotype = case_when(Genotype == "SDA-500" ~ "WT",
            Genotype == "CdKO Het" ~ "HET",
            Genotype == "CdKO Hom" ~ "HOM",
            .default = as.character(Genotype))) %>% 
  mutate(Sex = str_to_lower(Sex))


model <- glm(cbind(`Adults`,200-`Adults`) ~ line * Genotype * Sex, family = binomial, data = data)
summary(model)