#-------------------------------------------------------------------------------
# Whole-Exome sequencing analysis - visualization
#-------------------------------------------------------------------------------

library(maftools)


maf_files = list.files(path = "D:/JH_Lee/Processed_vcf/MAF/",pattern = "*.\\.maf", full.name=TRUE)
mymaf = maftools::merge_mafs(mafs = maf_files, verbose = TRUE)

#Filtration of mutation

mymaf_filt <- subsetMaf(maf = mymaf, query = "BIOTYPE=='protein_coding'")
mymaf_filt <- subsetMaf(maf = mymaf_filt, query = "IMPACT == 'HIGH'| IMPACT == 'MODERATE'")
mymaf_filt <- subsetMaf(maf = mymaf_filt, query = "FILTER == 'PASS'")

getSampleSummary(mymaf_filt)
getGeneSummary(mymaf_filt)
getClinicalData(mymaf_filt)
getFields(mymaf_filt)
write.mafSummary(maf=mymaf_filt,basename = 'maf_filt')


#Load reference mutation list
genelist_file <- read.csv("D:/JH_Lee/Processed_vcf/MAF/WES GENE LIST V3.csv")
genes <- genelist_file$gene
mymaf_filt_gene <- subsetMaf(maf = mymaf_filt, genes = genes)

getSampleSummary(mymaf_filt_gene)
getGeneSummary(mymaf_filt_gene)
getClinicalData(mymaf_filt_gene)
getFields(mymaf_filt_gene)
write.mafSummary(maf=mymaf_filt_gene,basename = 'maf_filt')

#visualization (Figure 2C, 2D)
order = c("_01-organoid","_01-parental","_02-organoid","_02-parental","_03-organoid","_03-parental","_04-organoid","_04-parental","_05-organoid","_05-parental","_06-organoid","_06-parental","_07-organoid","_07-parental","_08-organoid","_08-parental","_09-organoid","_09-parental")
plotmafSummary(maf = mymaf_filt_gene, rmOutlier = TRUE, addStat = 'median',dashboard = TRUE,titvRaw = FALSE)
oncoplot(maf=mymaf_filt_gene, top=47,sampleOrder = order, draw_titv = TRUE, gene_mar =6, fontSize = 0.5,removeNonMutated = FALSE, showTumorSampleBarcodes = TRUE,barcodeSrt=45)
oncoplot(maf=mymaf_filt, top=47,sampleOrder = order,logColBar = TRUE, draw_titv = TRUE, gene_mar =6, fontSize = 0.5,removeNonMutated = FALSE, showTumorSampleBarcodes = TRUE,barcodeSrt=45)
oncoplot(maf=mymaf_filt, top=47,sampleOrder = order,logColBar = FALSE, draw_titv = TRUE, gene_mar =6, fontSize = 0.5,removeNonMutated = FALSE, showTumorSampleBarcodes = TRUE,barcodeSrt=45)

#organoid
org_files = list.files(path = "D:/JH_Lee/Processed_vcf/MAF/Organoid/",pattern = "*.\\.maf", full.name=TRUE)
org = maftools::merge_mafs(mafs = org_files, verbose = TRUE)
org

org_filt <- subsetMaf(maf = org, query = "BIOTYPE=='protein_coding'")
org_filt <- subsetMaf(maf = org_filt, query = "IMPACT == 'HIGH'| IMPACT == 'MODERATE'")
org_filt <- subsetMaf(maf = org_filt, query = "FILTER == 'PASS'")

plotmafSummary(maf = org_filt, rmOutlier = TRUE, addStat = 'median',dashboard = TRUE,titvRaw = FALSE)

#parental
parental_files = list.files(path = "D:/JH_Lee/Processed_vcf/MAF/parental/",pattern = "*.\\.maf", full.name=TRUE)
parental = maftools::merge_mafs(mafs = parental_files, verbose = TRUE)
parental

parental_filt <- subsetMaf(maf = parental, query = "BIOTYPE=='protein_coding'")
parental_filt <- subsetMaf(maf = parental_filt, query = "IMPACT == 'HIGH'| IMPACT == 'MODERATE'")
parental_filt <- subsetMaf(maf = parental_filt, query = "FILTER == 'PASS'")

plotmafSummary(maf = parental_filt, rmOutlier = TRUE, addStat = 'median',dashboard = TRUE,titvRaw = FALSE)


#Mutation Concordance (Figure 2A)

library(dplyr)
data <- mymaf_filt@data

data <- data %>% select(Tumor_Sample_Barcode, Hugo_Symbol, Start_Position, Variant_Classification, Tumor_Seq_Allele1, Tumor_Seq_Allele2)

primary_data <- data %>% filter(grepl("-parental", Tumor_Sample_Barcode))
organoid_data <- data %>% filter(grepl("-organoid", Tumor_Sample_Barcode))

primary_data$Tumor_Sample_Barcode <- gsub("-.+", "", primary_data$Tumor_Sample_Barcode)
organoid_data$Tumor_Sample_Barcode <- gsub("-.+", "", organoid_data$Tumor_Sample_Barcode)

common_rows <- inner_join(primary_data, organoid_data, by = c("Tumor_Sample_Barcode", "Hugo_Symbol","Start_Position", "Variant_Classification", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"))

get_sample_number <- function(sample) {
  gsub("^.*-(\\d+).*$", "\\1", sample)
}

primary_data <- primary_data %>% mutate(Sample_Number = get_sample_number(Tumor_Sample_Barcode))
organoid_data <- organoid_data %>% mutate(Sample_Number = get_sample_number(Tumor_Sample_Barcode))


common_rows <- common_rows %>% mutate(Sample_Number = get_sample_number(Tumor_Sample_Barcode))

grouped_counts <- bind_rows(
  common_rows %>%
    group_by(Sample_Number) %>%
    summarize(concordant = n()),
  primary_data %>%
    anti_join(common_rows, by = c("Tumor_Sample_Barcode", "Hugo_Symbol","Start_Position", "Variant_Classification", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")) %>%
    group_by(Sample_Number) %>%
    summarize(primary_only = n()),
  organoid_data %>%
    anti_join(common_rows, by = c("Tumor_Sample_Barcode", "Hugo_Symbol","Start_Position", "Variant_Classification", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")) %>%
    group_by(Sample_Number) %>%
    summarize(organoid_only = n())
) %>%
  group_by(Sample_Number) %>%
  summarize(
    total = sum(concordant, primary_only, organoid_only, na.rm = TRUE),
    concordant = sum(concordant, na.rm = TRUE),
    primary_only = sum(primary_only, na.rm = TRUE),
    organoid_only = sum(organoid_only, na.rm = TRUE)
  ) %>%
  mutate(
    concordant_ratio = concordant / total,
    primary_only_ratio = primary_only / total,
    organoid_only_ratio = organoid_only / total
  ) %>%
  arrange(Sample_Number)

common_rows %>%
  group_by(Sample_Number) %>%
  summarize(concordant = n())

print(grouped_counts)
write.table(grouped_counts, file = "D:/JH_Lee/2023-09-27 WES/concordance.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#sample swap identification (Figure 2B)

bam_files = list.files(path = "D:/JH_Lee/2023-09-27 WES/1.bam/",pattern = "*.\\.bam", full.name=TRUE)
res = maftools::sampleSwaps(bams = bam_files, build = "hg19", ncores = 1)


res$pairwise_comparison
cor_table = cor(res$AF_table)
pheatmap::pheatmap(cor_table,breaks = seq(0,1,0.01),color = colorRampPalette(c("yellow","white","red"))(100))
