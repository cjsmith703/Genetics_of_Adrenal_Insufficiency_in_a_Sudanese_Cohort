#title "Genetic analysis of a Sudanese PAI cohort"
#author: Chris J Smith

#load packages
library(tidyverse)
library(stringr)
library(packcircles)
library(gtools)
library(openxlsx)

#Upload data
raw_data <- read.csv("redcap_data_anon_trials.csv", na.strings = c(""))
patient_codes <- read.csv("patient_codes.csv")
patient_codes <- patient_codes %>%
rename(study_id = patient)

#create functions for cleaning data
ethnicity_assigner <- function(value) {
  if (is.na(value)) {
    value <- "Unknown"
  } else if (value == "1") {
    return("Non-Finnish European")
  } else if (value == "2") {
    return("Finnish")
  } else if (value == "3") {
    return("African")
  } else if (value == "4") {
    return("Latino")
  } else if (value == "5") {
    return("South Asian")
  } else if (value == "6") {
    return("Ashkenazi Jewish")
  } else if (value == "7") {
    return("East Asian")
  } else if (value == "8") {
    return("Other")
  } else {
    return("Unknown")
  }
}

consanguinity_assigner <- function(value) {
  if (is.na(value)) {
    value <- "Unknown"
  } else if (value == "0") {
    return("Non-Consanguineous")
  } else if (value == "1") {
    return("Consanguineous")
  } else {
    return("Unknown")
  }
}

hyperpigmentation_assigner <- function(value) {
  if (is.na(value)) {
    value <- "Unknown"
  } else if (value == "0") {
    return("No Hyperpigmentation")
  } else if (value == "1") {
    return("Hyperpigmentation")
  } else {
    return("Unknown")
  }
}

empty_as_na <- function(x) {
  if ("factor" %in% class(x)) x <- as.character(x)
  ifelse(as.character(x) != "", x, NA)
}

solved_assigner <- function(value) {
  if (is.na(value)) {
    value <- "Unknown"
  } else if (value == "0") {
    return("Unsolved")
  } else if (value == "1") {
    return("Solved")
  } else {
    return("Unknown")
  }
}

#--------------------------------------
#clean data from raw redcap output

#data is for family data
data <- raw_data %>%
  mutate(gene1 = coalesce(mutation_diagnosis, mutation_diagnosis2)) %>% #mutation_diagnosis is post CGS, diagnosis2 is post WES
  mutate(gene = str_extract(gene1, "(MC2R|MRAP|STAR|NROB1|AIRE|CYP11A1|NNT|AAAS| 
                                    |ABCD1|SGPL1|MCM4|TXNRD2|POMC|CYP21A2|
                                    |POR|GPX1|CYP11B1|HSD3B2|ARSA)")) %>%
  mutate(ethnicity1 = sapply(ethnicity, ethnicity_assigner)) %>% #apply cleaning functions
  mutate(consanguinity = sapply(consaguinity, consanguinity_assigner)) %>%
  mutate(hyperpigmentation1 = sapply(hyperpigmentation,
                                    hyperpigmentation_assigner)) %>%
  filter(country == "Sudan") %>% #filter database for Sudanese patients
  filter(patient_or_family_member == 1) %>% #filter for affected proband
  mutate(mutation_solved = ifelse(is.na(mutation_solved) | #mutation_solved is solved at CGS
                                  !is.numeric(mutation_solved),
                                  0, mutation_solved)) %>%
  mutate(mutation_solved2 = ifelse(is.na(mutation_solved2) | #mutation_solved2 is solved after WES
                                  !is.numeric(mutation_solved2),
                                  0, mutation_solved2))

data$gender[data$gender___1 > 0] <- "Male"
data$gender[data$gender___2 > 0] <- "Female"

#patient_data is all patients including affected siblings
patient_data <- raw_data %>%
  mutate(gene1 = coalesce(mutation_diagnosis, mutation_diagnosis2)) %>% #mutation_diagnosis is post CGS, diagnosis2 is post WES
  mutate(gene = str_extract(gene1, "(MC2R|MRAP|STAR|NROB1|AIRE|CYP11A1|NNT|AAAS|
                                    |ABCD1|SGPL1|MCM4|TXNRD2|POMC|CYP21A2|
                                    |POR|GPX1|CYP11B1|HSD3B2|ARSA)")) %>%
  mutate(ethnicity1 = sapply(ethnicity, ethnicity_assigner)) %>% #apply cleaning functions
  mutate(consanguinity = sapply(consaguinity, consanguinity_assigner)) %>%
  mutate(hyperpigmentation1 = sapply(hyperpigmentation,
                                     hyperpigmentation_assigner)) %>%
  filter(country == "Sudan") %>% #filter database for Sudanese patients
  filter(patient_or_family_member == 1) %>% #filter for affected proband
  mutate(mutation_solved = ifelse(is.na(mutation_solved) | #mutation_solved is solved at CGS
                                  !is.numeric(mutation_solved),
                                  0, mutation_solved)) %>%
  mutate(mutation_solved2 = ifelse(is.na(mutation_solved2) | #mutation_solved2 is solved after WES
                                  !is.numeric(mutation_solved2),
                                  0, mutation_solved2))

patient_data$gender[patient_data$gender___1 > 0] <- "Male"
patient_data$gender[patient_data$gender___2 > 0] <- "Female"


#set colours for each gene
colours <- c("#C0392B", "#E74C3C", "#9B59B6",
             "#23c73e", "#8E44AD", "#2980B9",
             "#3498DB", "#1ABC8C", "#0f1de4",
             "#16A085", "#27AE60", "#2ECC71",
             "#F1C40F", "#F39C12", "#E37E22",
             "#D35400", "#BDC37C", "#7F8C8D",
             "#34495E")

names(colours) <- c("AAAS", "ABCD1", "AIRE",
                    "ARSA", "CYP11A1", "CYP11B1",
                    "CYP21A2", "GPX1", "HSD3B2",
                    "MC2R", "MRAP", "NNT",
                    "NROB1", "POMC", "POR",
                    "SGPL1", "STAR",
                    "TXNRD2")

#assign disease categories
#patient numbers removed for sharing purposes

autoimmune_patient_list <- c()
other_patient_list <- c()
exclude <- c()

pai_patients <- patient_data %>%
select(study_id) %>%
filter(!study_id %in% autoimmune_patient_list &
      !study_id %in% other_patient_list) %>%
filter(!study_id %in% exclude) %>%
unlist()

#function to add disease categorisation
patient_assigner <- function(value) {
  if (value %in% autoimmune_patient_list) {
    return("A-AI") #Autoimmune-Adrenal Insufficiency
    } else if (value %in% other_patient_list) {
    return("Other")
    } else if (value %in% pai_patients) {
    return("NA-AI") #Non-autoimmune Adrenal Insufficiency
    } else {
      return("Unknown")
    }
}

#Family solved data
family_solved <- data %>%
mutate(solved = ifelse(mutation_solved | mutation_solved2 == 1, #state whether solved/unsolved
                      "Solved", "Unsolved")) %>%
select(study_id, gender, solved, gene, gene1) %>%
inner_join(patient_codes, by = "study_id") %>% #assign number for published article
select(family, study_id, gender, solved, gene, gene1) %>%
arrange(family) %>%
mutate(family = round(family, digits = 0)) %>% #round number to produce family number
mutate(disease = sapply(study_id, patient_assigner))

family_solved_count <- family_solved %>% #count solved/unsolved
count(solved, disease)

#Patient solved data
patient_solved <- patient_data %>%
mutate(solved = ifelse(mutation_solved | mutation_solved2 == 1, #state whether solved/unsolved
                       "Solved", "Unsolved")) %>%
select(study_id, gender, solved, gene, gene1) %>%
inner_join(patient_codes, by = "study_id") %>% #assign number for published article
select(family, study_id, gender, solved, gene, gene1) %>%
arrange(family) %>%
mutate(disease = sapply(study_id, patient_assigner))

patient_solved_count <- patient_solved %>% #count solved/unsolved
count(solved)

#-----------------------------------------------
#Data Visualisation

#Bar chart of family solved data
ggplot(family_solved_count, aes(x = "", y = n, fill = solved)) +
  geom_bar(width = 1, stat = "identity") +
  theme_minimal() +
  coord_polar("y", start = 0) +
  ggtitle("Solved Rate") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#Solved rate by gender
solved_gender <- family_solved %>%
count(solved, gender)

ggplot(solved_gender, aes(x = gender, y = n, fill = solved)) +
  geom_col() +
  theme_minimal() +
   theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  ylab("Number of Families") +
  ggtitle("Solved Rate")

#Genes in our cohort pie chart
genes_count <- data %>%
  select(gene) %>%
  drop_na(gene) %>%
  count(gene)

genes_count <- genes_count %>%
arrange(desc(n))

ggplot(genes_count, aes(x = "", y = n, fill = gene)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = colours) +
  theme_minimal() +
  coord_polar("y", start = 0) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))+
  ggtitle("Genes of Cohort") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


#Cicle plot of genes in the cohort
packing <- circleProgressiveLayout(genes_count$n, sizetype = "area")

circle_plot <- cbind(genes_count, packing)
plot(circle_plot$radius, circle_plot$n)
dat_gg <- circleLayoutVertices(packing, npoints = 50)

ggplot() +
geom_polygon(data = dat_gg,
            mapping = aes(x, y, fill = as.factor(id), group = id, alpha = 0.6),
            colour = "white") +
geom_text(data = circle_plot, aes(x, y, size = n, label = gene)) +
scale_size_continuous(range = c(1, 13)) +
theme_void() +
theme(legend.position = "none") +
coord_equal()

#Bar chart of genes in our cohort
ggplot(genes_count, aes(x = reorder(gene, -n), n, fill = gene)) +
  geom_col() +
  scale_fill_manual(values = colours) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 30, face = "bold"),
        axis.title = element_text(size = 30, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ggtitle("Sudan Cohort PAI Genes") +
  ylab("Number of Diagnosed Families")

##Genes with genders
genes_genders <- data %>%
  select(gene, gender) %>%
  drop_na(gene) %>%
  drop_na(gender) %>%
  count(gene, gender)

ggplot(genes_genders, aes(x = reorder(gene, -n), n, fill = gender)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cohort PAI Genes by Gender") +
  ylab("Number of Diagnosed Families")


##Cortisol boxplot
cortisolug <- patient_data %>%
  select(study_id, gene, cortisol_level2) %>%
  drop_na(cortisol_level2)%>%
  arrange(desc(cortisol_level2)) %>%
  mutate(gene = ifelse(is.na(gene), "Unknown", gene))

cortisolug$cortisol_level2 <- as.numeric(cortisolug$cortisol_level2) #set values to numeric and check
class(cortisolug$cortisol_level2)

ggplot(cortisolug, aes(gene, y = cortisol_level2, fill = gene)) +
  scale_fill_manual(values = colours)+
  ylab("Cortisol Concentration (ug/dl)")+
  xlab(element_blank())+
  #ggtitle("Cortisol Concentration")+
  ylim(0,20)+
  theme_minimal()+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(alpha = 0.5)+
  geom_point(alpha = 0.5)+
  #geom_text(data = head(cortisol, 15), aes(label = record_id),colour = "black", nudge_x = 1, show.legend = FALSE ) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text=element_text(size=15, face = 'bold'),
        axis.title=element_text(size=15, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size=15, face = 'bold'), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none")

#Age of onset per gene
age_assigner <- function(value) { #convert all age of onset values to years
  if (grepl("birth|Birth|neonatal|Neonatal", value) == TRUE) {
    return(0.01)
  } else if (grepl("years|year|YEARS", value) == TRUE) {
    return(as.numeric(str_extract(value, "\\d+")))
  } else if (grepl("month|months|m", value) == TRUE) {
    return(as.numeric(str_extract(value, "\\d+")) / 12)
  } else if (grepl("weeks", value) == TRUE) {
    return(as.numeric(str_extract(value, "\\d+")) / 52)
  } else if (grepl("days|day", value) == TRUE) {
    return(as.numeric(str_extract(value, "\\d+")) / 365)
  } else if (grepl("hours", value) == TRUE) {
    return(as.numeric(str_extract(value, "\\d+")) / (365 * 24))
  } else if (grepl("?", value) == TRUE) {
    return(as.numeric(str_extract(value, "\\d+")))
  } else {
    return(value)
  }
}

age <- patient_data %>%
  select(study_id, gene, age_of_onset) %>%
  drop_na(age_of_onset) %>%
  mutate(age_years = sapply(age_of_onset, age_assigner)) %>%
  drop_na(age_years) %>%
  arrange(gene)

ggplot(age, aes(gene, age_years, fill = gene)) +
  scale_fill_manual(values = colours) +
  ylab("Age (years)") +
  xlab(element_blank()) +
  ggtitle("Age of Onset") +
  theme_minimal() +
  ylim(0, 24) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

#----------------------------------------------------------------
#Statistical analysis of disease categorisation

#set colours for solved/unsolved and disease category
solved_colours <- c("#e2c4c0", "#a7c5ec")
names(solved_colours) <- c("Unsolved", "Solved")

disease_colours <- c("#297992", "#70d1c1")
names(disease_colours) <- c("NA-AI", "A-AI")

my_comparisons <- list( c("NA-AI", "A-AI"))

#Analyse solved rates between disease categories
solved_2 <- family_data %>%
mutate(solved = ifelse(mutation_solved |
                      mutation_solved2 == 1,
                      "Solved", "Unsolved")) %>%
select(study_id, disease, solved, gene, gene1) %>%
arrange(solved, gene)

solved_disease <- solved_2 %>%
count(solved, disease)

#Set categories as factors
solved_disease$solved <- factor(solved_disease$solved,
                                levels = c("Unsolved", "Solved"),
                                ordered = TRUE)

solved_disease$disease <- factor(solved_disease$disease,
                                 levels = c("NA-AI", "A-AI"),
                                 ordered = TRUE)

#Bar chart of solved categories
ggplot(solved_disease, aes(x = disease, y = n, fill = solved)) +
  scale_fill_manual(values = solved_colours) +
  geom_col() +
  theme_minimal() +
   theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  ylab("Number of Families") +
  ggtitle("Solved Rate")

#Cortisol plot of solved categories 
cortisolug <- patient_data %>% 
  select(record_id, gene, disease, cortisol_level2, study_id) %>%
  drop_na(cortisol_level2)%>%
  arrange(desc(cortisol_level2)) 

cortisolug$disease <- factor(cortisolug$disease, #set disease categories as factors
                                 levels = c('NA-AI', 'A-AI'), ordered=TRUE)

ggplot(cortisolug, aes(disease, y = cortisol_level2, fill = disease)) + #plot boxplot with anova analysis
  scale_fill_manual(values = disease_colours)+
  ylab("Cortisol Concentration (ug/dl)")+
  xlab(element_blank())+
  ggtitle("Cortisol Concentration")+
  ylim(0,20)+
  theme_minimal()+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(alpha = 0.5)+
  geom_point(alpha = 0.5)+
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text=element_text(size=15, face = 'bold'),
        axis.title=element_text(size=15, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size=15, face = 'bold'), 
        legend.position="none")+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(method = "t.test")

#t test for na-ai vs a-ai
t.test(cortisol_level2 ~ disease, data = cortisolug)

#age of onset between categories
age$disease <- factor(age$disease, 
                      levels = c("PAI", "Autoimmune", "Transient", "Unknown"),
                      ordered = TRUE)

#t test for na-ai vs a-ai
t.test(age_years ~ disease, data = age)

#plot age of onset with all data points and statistical comparisons
ggplot(age, aes(disease, age_years, fill = disease)) +
  scale_fill_manual(values = disease_colours) +
  ylab("Age (years)") +
  xlab(element_blank()) +
  ggtitle("Age of Onset") +
  theme_minimal() +
  ylim(0, 24) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 20)+
  stat_compare_means(method = "t.test")
