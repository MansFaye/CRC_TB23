library(tidyverse)

## filtering segments and genes above limits of quantitation (LOQ)

# read annot data
data_annot <- read.csv("Master files/TB_TME_AL_target8.1_annot.csv", check.names=FALSE)

# filter to get df with epithelium TB with % alligned reads >=60%
data_annot_tb_tme <- data_annot %>% filter(Compartment=="TME" & TB_AL=="TB" & `% aligned`>=60)

# read raw data 
data_raw <- read.csv("Master files/TB_TME_AL_target8.1_count.csv", check.names = FALSE)

# filter data_raw using segmentdisplayname in data_annot_tb_e
data_raw_tb_tme <- data_raw %>% select(TargetName, data_annot_tb_e$SegmentDisplayName)

write.csv(data_annot_tb_tme, "data_annot_tb_tme.csv")
write.csv(data_raw_tb_tme, "data_raw_tb_tme.csv")


#---------------------------------------------------------
## CREATE LOOP 
# FOR EACH SAMPLE - CHECK WHICH GENES ARE >LOQ
# LOQ LIST

loq_data <- data_annot_tb_e %>% select(SegmentDisplayName, `LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0)`)
names(loq_data)[2] <- "loq"


# list of genes
genes <- data_raw_tb_e %>% select(TargetName)
samples <- data_annot_tb_e$SegmentDisplayName
results <- genes

for (i in 1:length(samples)) {
  
  sample <- samples[[i]]
  
  sample_data <- data_raw_tb_e %>% select(TargetName, samples[[i]])
  
  #if value below threshold replace with na
  sample_data.1 <- sample_data %>% 
    mutate("temp"=ifelse(sample_data[,2]>loq_data$loq, 1, 0))
  
  #rename thresholded column
  names(sample_data.1)[names(sample_data.1) == "temp"] <- sample
  
  #delete middle column of sample_data.1
  sample_data.2 <- sample_data.1[,-2]
  
  #left join to genes df
  thresholded_data <- left_join(sample_data.2, genes, by="TargetName")
  
  thresholded_data.1 <- thresholded_data %>% select(samples[[i]])
  
  results <- cbind(thresholded_data.1, results)
  
}

results1 <- results

## get segment level summary of genes over LOQ
segment_genes_loq <- as.data.frame(colSums(results[,1:29], na.rm = FALSE, dims = 1))
names(segment_genes_loq)[names(segment_genes_loq) == "colSums(results[, 1:29], na.rm = FALSE, dims = 1)"] <- "N genes above LOQ"
segment_genes_loq <- segment_genes_loq %>% mutate("% of panel genes (n=18677)"=`N genes above LOQ`/18677*100)

## filter out segments with <5% of genes over LOQ
segment_genes_loq_5 <- segment_genes_loq %>% filter(`% of panel genes (n=18677)`>=5)
segment_genes_loq_5 <- tibble::rownames_to_column(segment_genes_loq_5, "SegmentDisplayName")
results_f5 <- results %>% select(TargetName, segment_genes_loq_5$SegmentDisplayName)

library(janitor)
# transpose df 
results_t <- as.data.frame(t(results_f5))
colnames(results_t) <- results_t[1,]
results_t <- results_t[-1, ] 
results_t <- tibble::rownames_to_column(results_t, "SegmentDisplayName")
temp <- data_annot_tb_e %>% select(subgroup.a, SegmentDisplayName) 
results_t <- left_join(temp, results_t, by="SegmentDisplayName")
results_t <- drop_na(results_t)
results_t[ ,3:18679] <- sapply(results_t[ ,3:18679],as.numeric)

table(results_t$subgroup.a)

# bio group agnostic summary
#--------------------------------
gene_segments <- results_t %>%
  summarise_if(is.numeric, sum)

gene_segments.1 <- gene_segments/26*100
gene_segments.2 <- as.data.frame(t(gene_segments.1))
table(gene_segments.2)
names(gene_segments.2)[names(gene_segments.2) == "V1"] <- "Percent segments over LOQ"

## assess number of bio groups with genes above 10%
gene_segments.2 <- gene_segments.2 %>%
  mutate("genes over LOQ in n groups"=ifelse(`Percent segments over LOQ`>=20, 1, 0))
table(gene_segments.2$`genes over LOQ in n groups`)


## TMP 
TMP_loq_freq <- results_t %>%
  filter(subgroup.a=="TMP") %>%
  summarise_if(is.numeric, sum)

TMP_loq_freq.1 <- TMP_loq_freq/9*100
TMP_loq_freq.2 <- as.data.frame(t(TMP_loq_freq.1))
table(TMP_loq_freq.2)
names(TMP_loq_freq.2)[names(TMP_loq_freq.2) == "V1"] <- "TMP"


## TSM 
TSM_loq_freq <- results_t %>%
  filter(subgroup.a=="TSM") %>%
  summarise_if(is.numeric, sum)

TSM_loq_freq.1 <- TSM_loq_freq/7*100
TSM_loq_freq.2 <- as.data.frame(t(TSM_loq_freq.1))
table(TSM_loq_freq.2)
names(TSM_loq_freq.2)[names(TSM_loq_freq.2) == "V1"] <- "TSM"

## TSS 
TSS_loq_freq <- results_t %>%
  filter(subgroup.a=="TSS") %>%
  summarise_if(is.numeric, sum)

TSS_loq_freq.1 <- TSS_loq_freq/10*100
TSS_loq_freq.2 <- as.data.frame(t(TSS_loq_freq.1))
table(TSS_loq_freq.2)
names(TSS_loq_freq.2)[names(TSS_loq_freq.2) == "V1"] <- "TSS"

# bind all three columns
gene_loq_freq_biogroup <- cbind(TMP_loq_freq.2, TSM_loq_freq.2, TSS_loq_freq.2)

## assess number of bio groups with genes above 10%
gene_loq_freq_biogroup <- gene_loq_freq_biogroup %>%
  mutate("genes over LOQ in n groups"=ifelse(TSM>=10 & TMP>=10 & TSS>=10, 3,
                                             ifelse(TSM<10 & TMP<10 & TSS<10, 0,
                                                    ifelse(TSM>=10 & TMP>=10 & TSS<10|
                                                             TSM>=10 & TMP<10 & TSS>=10|
                                                             TSM<10 & TMP>=10 & TSS>=10, 2, 1))))

table(gene_loq_freq_biogroup$`genes over LOQ in n groups`)








