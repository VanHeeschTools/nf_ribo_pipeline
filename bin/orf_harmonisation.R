library(dplyr)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(ggplot2)
library(BSgenome.Homo.sapiens.RMScontainer)


basedir = "/hpc/pmc_vanheesch/projects/jvandiner/rms_analysis"

rna_dir = paste(basedir,
                "01_rnaeq",sep="/")
ribo_dir = paste(basedir,
                 "02_ribseq", sep = "/")

txdb_loc = paste0(rna_dir,"/analysis/rnaseqpipeline/cutomannotation/",
                   "RMS_conainer/",
                   "RMS_full_novel_fitered_corrected.gtf_TxDb")



## Data load
#RMS transcriptome

rms_gtf <- as.data.frame(rtracklayer::import.gff(paste(rna_dir,"analysis","rnaseq_pipeline","customannotation","RMS_full_novel_filtered_corrected.sorted.gtf", sep = "/")))

# Ensembl CDS
txdb <- AnnotationDbi::loadDb(txdb_loc)
tx2gene <- AnnotationDbi::select(txdb, 
                  AnnotationDbi::keys(txdb, keytype = "TXNAME"), 
                  "GENEID", 
                  "TXNAME") %>%
  dplyr::rename(transcript_id = TXNAME, gene_id = GENEID)
main_cds <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = T)
main_cds_sequence <- GenomicFeatures::extractTranscriptSeqs(x = BSgenome.Homo.sapiens.RMScontainer,
                                                   transcripts = main_cds)
main_cds_sequence <- main_cds_sequence[which((width(main_cds_sequence) %% 3) == 0),]
cds_proteins <- Biostrings::translate(main_cds_sequence) %>%
  as.data.frame() %>%
  dplyr::rename(protein = x) %>%
  tibble::rownames_to_column(var = "transcript_id") %>%
  dplyr::left_join(tx2gene,
                   by = "transcript_id") %>%
  dplyr::filter(grepl("^M",protein))



#First run PRICE_annotation.Rmd and ORFQUANT_annotation.Rmd with their attached
#scripts to calculate PPM, CDS similarity, protein similarity, and annotate 
#ORF-specific things.


price_orf_table <- read.delim(paste(ribo_dir,
                         "results/orf_reannotation",
                         "price_orf_table.csv",
                         sep = "/"), sep = ",") %>%
  dplyr::mutate(ORF_ranges = paste0(seqnames,":",start,"-",end),
                translated = ifelse(number_patient_samples > 3, T,F)) %>%
  dplyr::select(-c("seqnames","start","end","width","strand","RMS_enriched","uniprot_gn_ids")) %>%
  dplyr::rename(orf_id = name) %>%
  dplyr::filter(!(is.na(ppm_sum)))

orfquant_table <- read.delim(paste(ribo_dir,
                         "results/orf_reannotation",
                         "RMS_ORFquant_table.csv",
                         sep = "/"), sep = ",") %>%
  dplyr::select(-c("expr_threshold","start_dif",
                   "cdsRangeSimilarity","stop_same","orf_cds",
                   "start_check","same_as_cds",
                   "P_sites_raw","P_sites_raw_uniq","transcript_biotype",
                   "ORF_category_Tx","ORF_category_Gen")) %>%
  dplyr::rename(orf_id = ORF_id_tr) %>%
  dplyr::mutate(translated = ifelse(number_patient_samples > 3, T,F),
                start_codon = "ATG")

rms_orf_table <- rbind(orfquant_table,
                       price_orf_table)


# RMS enrichment

rms_gtf <- as.data.frame(rtracklayer::import.gff(paste(rna_dir,"analysis","rnaseq_pipeline","customannotation","RMS_full_novel_filtered_corrected.sorted.gtf", sep = "/")))

fp_table <- read.delim(paste(rna_dir,
                                "results/quantification/",
                                "RMS-FP_enrichment.csv", sep="/"), 
                          sep = ",")
fp_enriched <- fp_table %>% 
                      dplyr::filter(enriched == T &
                                    gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
                      dplyr::pull(gene_id)

fn_table <- read.delim(paste(rna_dir, 
                                "results/quantification/",
                                "RMS-FN_enrichment.csv", sep="/"), 
                          sep = ",")
fn_enriched <- fn_table %>% 
                      dplyr::filter(enriched == T &
                                    gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
                      dplyr::pull(gene_id)

rms_table <- read.delim(paste(rna_dir, "results/quantification/",
                                 "RMS_enrichment.csv", sep="/"), 
                           sep = ",")
rms_enriched <- rms_table %>% 
                      dplyr::filter(enriched == T &
                                    gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
                      dplyr::pull(gene_id)

enriched_genes <- c(fp_enriched, 
                    rms_enriched,
                    fn_enriched) %>% unique()

fp_specific <- fp_table %>% 
                      dplyr::filter(specific == T &
                                    gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
                      dplyr::pull(gene_id)

fn_specific <- fn_table %>% 
                      dplyr::filter(specific == T &
                                    gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
                      dplyr::pull(gene_id)

rms_specific <- rms_table %>% 
                      dplyr::filter(specific == T &
                                    gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
                      dplyr::pull(gene_id)

specific_genes <- c(fp_specific, 
                    rms_specific,
                    fn_specific) %>% unique()

rms_orf_table <- rms_orf_table %>%
  dplyr::mutate(RMS_enriched = ifelse(gene_id %in% enriched_genes,T,F),
                RMS_specific = ifelse(gene_id %in% specific_genes,T,F))

## intORFs reannotation

#We want to only select out of frame intORFs. Some nested ORFs are in the same frame
#as the CDS region


intorfs = read.delim(paste(
    ribo_dir,
    "results/orf_reannotation",
    "jorge_psites_intorfs.csv",
    sep = "/"),
  sep = ",")

rms_orf_table <- rms_orf_table %>%
  dplyr::mutate(orf_category_new = 
                  dplyr::case_when(orf_category_new == "nested_ORF" & 
                                     !(orf_id %in% intorfs$orf_id) ~ "not_nested",
                                   T ~ orf_category_new))

table(rms_orf_table_ints$orf_category_new)

# Load P site counts
ppm_mat <- read.table(paste(ribo_dir,"analysis/orf_quantification",
                            "ppm_orfquant.txt", sep = "/"),
            sep = ",", check.names = F) %>%
  dplyr::select(1:50) %>%
  dplyr::mutate(ppm_sum = rowSums(.),
                ppm_mean = rowMeans(.),
                ppm_median = Biobase::rowMedians(as.matrix(.))) %>%
  tibble::rownames_to_column(var = "orf_id") %>%
  dplyr::filter(orf_id %in% c(nested_orfs$orf_id, main_nested_orfs$orf_id)) %>%
  dplyr::select(orf_id,ppm_mean,ppm_median,ppm_sum)

ppm_mat_price <- read.table(paste(ribo_dir,"analysis/orf_quantification",
                            "ppm_price.txt", sep = "/"),
            sep = ",", check.names = F) %>%
  dplyr::select(1:50) %>%
  dplyr::mutate(ppm_sum = rowSums(.),
                ppm_mean = rowMeans(.),
                ppm_median = Biobase::rowMedians(as.matrix(.))) %>%
  tibble::rownames_to_column(var = "orf_id") %>%
  dplyr::filter(orf_id %in% c(nested_orfs$orf_id, main_nested_orfs$orf_id)) %>%
  dplyr::select(orf_id,ppm_mean,ppm_median,ppm_sum)

ppm_mat <- rbind(ppm_mat,ppm_mat_price)

data.frame(orf_id_i = nested_orfs$orf_id,
           transcript_id_i = nested_orfs$transcript_id,
           gene_id = nested_orfs$gene_id,
           sim_score = nested_orfs$similarity_score) %>%
  dplyr::left_join(ppm_mat, by = c("orf_id_i" = "orf_id")) %>%
  dplyr::rename(nested_mean = ppm_mean,
                nested_median = ppm_median,
                nested_sum = ppm_sum) %>%
  dplyr::left_join(main_nested_orfs,
                   by = "gene_id") %>%
  dplyr::select(-ppm_mean) %>%
  dplyr::left_join(ppm_mat, by = "orf_id") %>%
  dplyr::mutate(nested_mean = log2(nested_mean + 0.01),
                nested_median = log2(nested_median + 0.01),
                nested_sum = log2(nested_sum + 0.01),
                ppm_mean = log2(ppm_mean + 0.01),
                ppm_median = log2(ppm_median + 0.01),
                ppm_sum = log2(ppm_sum + 0.01)) %>%
  dplyr::filter(!(is.na(nested_mean) |
                  is.na(ppm_mean))) %>%
  ggplot() +
  geom_density(aes(x = ppm_median), 
               col = "darkred") +
  geom_density(data = . %>% 
                 dplyr::filter(sim_score > 66),
               aes(x = nested_median), 
               col = "#FCD667",
               fill = "#FCD667",
               alpha=.33) +
  geom_density(data = . %>% 
                 dplyr::filter(sim_score < 40),
               aes(x = nested_median), 
               col = "#BF8FCC",
               fill = "#BF8FCC",
               alpha=.33) +
  theme_minimal()



## CDS reannotation

#Load CDS overlap 


orfquant_cds <- read.delim(paste(ribo_dir,
                       "results/orf_reannotation",
                       "orfquant_cds_similarity.csv",
                       sep= "/"), 
                       sep = ",", 
                       row.names = 1) %>%
  dplyr::mutate(start_codon = "ATG")

price_cds <- read.delim(paste(ribo_dir,
                       "results/orf_reannotation",
                       "price_cds_similarity.csv",
                       sep= "/"), 
                       sep = ",",
                       row.names = 1)
# Don't use this one for now
price_orfquant <- read.delim(paste(ribo_dir,
                       "results/orf_reannotation",
                       "price_orfquant_similarity.csv",
                       sep= "/"), 
                       sep = ",",
                       row.names = 1)

cds_check <- rbind(orfquant_cds,price_cds)

#Check CDS regions that are very similar to called ORFs based on CDS loci and 
#annotated proteins.
#Check whether PRICE ORF calls overlap with ORFquant calls

#```{r CDS annotation}
rms_orf_table <- rms_orf_table %>%
  dplyr::left_join(cds_check[,c("orf_id","same_as_cds")]) %>%
  dplyr::mutate(same_as_cds = tidyr::replace_na(same_as_cds,
                                                replace = 0),
                similarity_score = tidyr::replace_na(similarity_score,
                                                     replace = 0),
                orf_category_new = ifelse(same_as_cds == T,"ORF_annotated",
                                          ifelse(similarity_score > 80, "ORF_annotated",orf_category_new))) %>%
  dplyr::select(-same_as_cds)

rms_orf_table <- rms_orf_table %>%
  dplyr::mutate(orf_caller = ifelse(grepl("__",orf_id),"PRICE","ORFquant"),
                orf_category_new = ifelse(gene_biotype == "lncRNA","lncORF",orf_category_new))

 write.table(rms_orf_table, file = paste(ribo_dir,
                                             "results/orf_reannotation",
                                             "RMS_harmonised_ORF_table.csv",
                                             sep = "/"),
             sep = ",", 
             quote = F, 
             row.names = F)

#```


## Plot ORF calling comparisons

#```{r }

table(rms_orf_table[which(rms_orf_table$translated == T & rms_orf_table$RMS_enriched == T),]$orf_category_new,
      rms_orf_table[which(rms_orf_table$translated == T & rms_orf_table$RMS_enriched == T),]$orf_caller)

table(rms_orf_table[which(rms_orf_table$translated == T),]$orf_category_new,
      rms_orf_table[which(rms_orf_table$translated == T),]$start_codon)

rms_orf_table %>%
  dplyr::filter(orf_category_new %in% c("novel","ORF_annotated","nested_ORF","uORF","overl_dORF","overl_uORF","dORF")) %>%
  ggplot(aes(x = orf_category_new, fill = orf_caller)) +
  geom_bar(position = "dodge") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("All called ORFs")

rms_orf_table %>%
  dplyr::filter(orf_category_new %in% c("novel","ORF_annotated","nested_ORF",
                                        "uORF","overl_dORF","overl_uORF","dORF") &
                translated == T) %>%
  ggplot(aes(x = orf_category_new, fill = orf_caller)) +
  geom_bar(position = "dodge") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("All translated ORFs")

rms_orf_table %>%
  dplyr::filter(orf_category_new %in% c("novel","ORF_annotated","nested_ORF",
                                        "uORF","overl_dORF","overl_uORF","dORF") &
                translated == T &
                RMS_enriched == T) %>%
  ggplot(aes(x = orf_category_new, fill = orf_caller)) +
  geom_bar(position = "dodge") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("All translated ORFs in enriched RMS genes")


#```

## Check

#Some checks
#```{r}

# Check genes with translated ORFs
table(rms_orf_table %>%
  dplyr::filter(translated == T) %>%
    dplyr::distinct(gene_id,gene_biotype) %>%
    dplyr::select(gene_biotype))

# Check translated ORF categories
table(rms_orf_table %>%
  dplyr::filter(translated == T) %>%
    dplyr::select(orf_category_new))

# Check novel ORFs in novel genes
nrow(rms_orf_table %>% dplyr::filter(translated == T & 
                                       gene_biotype == "stringtie"))

# Check lncORFs
nrow(rms_orf_table %>% dplyr::filter(translated == T & 
                                       gene_biotype == "lncRNA"))

# Check translated ORF in RMS enriched genes categories
rms_spec_trans_orf_table <- rms_orf_table %>%
  dplyr::filter(translated == T & RMS_enriched == T)

table(rms_spec_trans_orf_table %>%
  dplyr::filter(translated == T) %>%
    dplyr::select(orf_category_new))

# Check novel ORFs in novel genes
nrow(rms_spec_trans_orf_table %>% dplyr::filter(gene_biotype == "stringtie"))

# Check lncORFs
nrow(rms_spec_trans_orf_table %>% dplyr::filter(gene_biotype == "lncRNA"))

# Check enriched RMS genes with ORFs
table(rms_spec_trans_orf_table %>%
    dplyr::distinct(gene_id,gene_biotype) %>%
    dplyr::select(gene_biotype))


## PPM plots


ppm_price <- read.delim(paste(ribo_dir,"analysis/orf_quantification",
                            "ppm_price.txt", sep = "/"),sep =",")

ppm_orfquant <- read.delim(paste(ribo_dir,"analysis/orf_quantification",
                            "ppm_orfquant.txt", sep = "/"),sep =",") %>%
  dplyr::select(-RMS_merged)

ppms <- rbind(ppm_price,ppm_orfquant)

colors <- c("#0073C2FF", "#EFC000FF", "#DC0000FF", "#7AA6DCFF", "#FF7F0EFF", "#17BCEFFF", "#009E73")

## Translated
data.frame(means = rowMeans(rbind(ppm_price,ppm_orfquant)),
                        orf_id = rownames(rbind(ppm_price,ppm_orfquant))) %>%
  dplyr::filter(orf_id %in% rms_orf_table[which(rms_orf_table$translated == T),]$orf_id) %>%
  dplyr::left_join(rms_orf_table[,c("orf_id","orf_category_new")]) %>%
  dplyr::filter(!(grepl("extension|truncation",orf_category_new))) %>%
  ggplot(aes(x = orf_category_new, y = means, fill = orf_category_new, color = orf_category_new)) +
  ggbeeswarm::geom_quasirandom(size = 0.1, color = "grey20", fill = "grey20", shape = 16) +
  scale_y_continuous(trans = "log2", labels = scales::comma) +
  coord_flip() +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", width = 0.6) +
  scale_fill_manual(values = colors) +
  labs(x = "Translated ORFs", y = "Normalized P-site Counts") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor.y = element_blank())

ggsave(file = paste(ribo_dir,
                    "results/quantification/figures",
                    "PPM_per_ORF_type.pdf",
                    sep = "/"),device = "pdf",
       height = 4,
       width = 5)

## Translated Enriched
data.frame(means = rowMeans(rbind(ppm_price,ppm_orfquant)),
                        orf_id = rownames(rbind(ppm_price,ppm_orfquant))) %>%
  dplyr::filter(orf_id %in% rms_orf_table[which(rms_orf_table$translated == T),]$orf_id) %>%
  dplyr::left_join(rms_orf_table[,c("orf_id","orf_category_new","RMS_enriched")]) %>%
  dplyr::filter(!(grepl("extension|truncation",orf_category_new)) &
                  RMS_enriched == T) %>%
  ggplot(aes(x = orf_category_new, y = means, fill = orf_category_new, color = orf_category_new)) +
  ggbeeswarm::geom_quasirandom(size = 0.1, color = "grey20", fill = "grey20", shape = 16) +
  scale_y_continuous(trans = "log2", labels = scales::comma) +
  coord_flip() +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", width = 0.6) +
  scale_fill_manual(values = colors) +
  labs(x = "Translated enriched ORFs", y = "Normalized P-site Counts") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor.y = element_blank())

ggsave(file = paste(ribo_dir,
                    "results/quantification/figures",
                    "PPM_per_ORF_type_RMS_enriched.pdf",
                    sep = "/"),device = "pdf",
       height = 4,
       width = 5)


#### Trans no-Mean


ppms %>%
  tibble::rownames_to_column(var = "orf_id") %>%
  tidyr::pivot_longer(cols = -orf_id, 
                      names_to = "sample_id") %>%
  dplyr::filter(orf_id %in% rms_orf_table[which(rms_orf_table$translated == T &
                                                rms_orf_table$RMS_enriched == T),]$orf_id) %>%
  dplyr:: mutate(type = ifelse(grepl("ORG.*FP",sample_id),"ORG_FP",
                        ifelse(grepl("ORG.*FN",sample_id),"ORG_FN",
                        ifelse(grepl("TIS.*aRMS",sample_id),"PAT_FP",
                        ifelse(grepl("TIS.*eRMS",sample_id),"PAT_FN",
                        "other"))))) %>%
dplyr::left_join(rms_orf_table[,c("orf_id","orf_category_new","RMS_enriched")]) %>%
  dplyr::filter(!(grepl("extension|truncation",orf_category_new))) %>%
  dplyr::mutate(orf_category_new = ifelse(orf_category_new == "lncORF",
                                          "novel",orf_category_new)) %>%
  ggplot(aes(x = orf_category_new, y = value, fill = orf_category_new, color = orf_category_new)) +
  geom_hline(yintercept = 4) +
  ggbeeswarm::geom_quasirandom(size = 0.1, color = "grey20", fill = "grey20", shape = 16) +
  scale_y_continuous(trans = "log2", labels = scales::comma) +
  coord_flip() +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", width = 0.6) +
  scale_fill_manual(values = colors) +
  labs(x = "Translated enriched ORFs", y = "Normalized P-site \n Unique sample counts") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~ type, ncol = 2)
                 

ggsave(file = paste(ribo_dir,
                    "results/quantification/figures",
                    "PPM_per_ORF_per_type_per_sample_RMS_enriched.pdf",
                    sep = "/"),device = "pdf",
       height = 4,
       width = 5)


#### Trans Enr Demarc

org_fp_means <- rowMeans(ppms[grepl("ORG.*FP",colnames(ppms)),])
org_fn_means <- rowMeans(ppms[grepl("ORG.*FN",colnames(ppms)),])
tis_fn_means <- rowMeans(ppms[grepl("TIS.*eRMS",colnames(ppms)),])
tis_fp_means <- rowMeans(ppms[grepl("TIS.*aRMS",colnames(ppms)),])

## Translated -> demarcated
data.frame(means = c(org_fp_means,org_fn_means,
                     tis_fp_means,tis_fn_means),
           orf_id = rownames(ppms),
           type = c(rep("ORG-FP",length(org_fp_means)),
                    rep("ORG-FN",length(org_fn_means)),
                    rep("PAT-FP",length(tis_fp_means)),
                    rep("PAT-FN",length(tis_fn_means))),
           batch = c(rep("ORG",length(org_fp_means)+length(org_fn_means)),
                    rep("PAT",length(tis_fp_means) + length(tis_fn_means))
                    )) %>%
  dplyr::filter(orf_id %in% rms_orf_table[which(rms_orf_table$translated == T),]$orf_id) %>%
  dplyr::left_join(rms_orf_table[,c("orf_id","orf_category_new","RMS_enriched")]) %>%
  dplyr::mutate(orf_category_new = ifelse(orf_category_new == "lncORF",
                                          "novel",orf_category_new)) %>%
  dplyr::filter(!(grepl("extension|truncation",orf_category_new))) %>%
  ggplot(aes(x = orf_category_new, y = means, fill = orf_category_new, color = orf_category_new)) +
  geom_hline(yintercept = 4) +
  ggbeeswarm::geom_quasirandom(size = 0.1, color = "grey20", fill = "grey20", shape = 16) +
  scale_y_continuous(trans = "log2", labels = scales::comma) +
  coord_flip() +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", width = 0.6) +
  scale_fill_manual(values = colors) +
  labs(x = "Translated enriched ORFs", y = "Normalized P-site Counts") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~ batch, ncol = 1)

ggsave(file = paste(ribo_dir,
                    "results/quantification/figures",
                    "PPM_per_batch_per_ORF_type.pdf",
                    sep = "/"),device = "pdf",
       height = 4,
       width = 5)

## Translated Enriched -> demarcated
data.frame(means = c(org_fp_means,org_fn_means,
                     tis_fp_means,tis_fn_means),
           orf_id = rownames(ppms),
           type = c(rep("ORG-FP",length(org_fp_means)),
                    rep("ORG-FN",length(org_fn_means)),
                    rep("PAT-FP",length(tis_fp_means)),
                    rep("PAT-FN",length(tis_fn_means))),
           batch = c(rep("ORG",length(org_fp_means)+length(org_fn_means)),
                    rep("PAT",length(tis_fp_means) + length(tis_fn_means))
                    )) %>%
  dplyr::filter(orf_id %in% rms_orf_table[which(rms_orf_table$translated == T),]$orf_id) %>%
  dplyr::left_join(rms_orf_table[,c("orf_id","orf_category_new","RMS_enriched")]) %>%
  dplyr::filter(!(grepl("extension|truncation",orf_category_new)) &
                  RMS_enriched == T) %>%
  ggplot(aes(x = orf_category_new, y = means, fill = orf_category_new, color = orf_category_new)) +
  ggbeeswarm::geom_quasirandom(size = 0.1, color = "grey20", fill = "grey20", shape = 16) +
  scale_y_continuous(trans = "log2", labels = scales::comma) +
  coord_flip() +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", width = 0.6) +
  scale_fill_manual(values = colors) +
  labs(x = "Translated enriched ORFs", y = "Normalized P-site Counts") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~ batch, ncol = 1)

ggsave(file = paste(ribo_dir,
                    "results/quantification/figures",
                    "PPM_per_batch_per_ORF_type_RMS_enriched.pdf",
                    sep = "/"),device = "pdf",
       height = 4,
       width = 5)

# colors <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#DC0000FF", "#7AA6DCFF", "#1B1B1BFF", "#FF7F0EFF", "#17BCEFFF", "#009E73", "#CC79A7")

counts_summary <- data.frame(mean = rowMeans(counts(dds_NC, normalized = TRUE))) %>%
  rownames_to_column("orf_id") %>%
  left_join(orfs_tx_df[, c("orf_id", "orf_category_new")], by = c("orf_id" = "orf_id")) %>%
  column_to_rownames("orf_id")

ORF_df <- data.frame(orf_id = ORFs_translated) %>%
  left_join(orfs_tx_df, by = c("orf_id" = "orf_id"))

ORF_df_summary <- ORF_df %>%
  group_by(orf_category_new) %>%
  summarize(count = n())

df_row <- data.frame(orf_id = rownames(mat_vsd_NC)) %>%
  left_join(orfs_tx_df[, c("orf_id", "orf_category_new")], by = c("orf_id" = "orf_id")) %>%
  column_to_rownames("orf_id") %>%
  filter(!grepl("extension|truncation", .$orf_category_new))

# Exclude ORF types with 'extension' or 'truncation' in their names
counts_summary <- counts_summary[!grepl("extension|truncation", counts_summary$orf_category_new), ]
ORF_df_summary <- ORF_df_summary[!grepl("extension|truncation", ORF_df_summary$orf_category_new), ]
# df_row <- df_row[!grep("extension|truncation", df_row$orf_category_new), ]

counts_summary$orf_category_new <- factor(counts_summary$orf_category_new, levels = c(ORF_df_summary[order(ORF_df_summary$count, decreasing = TRUE), ]$orf_category_new))
counts_summary$orf_category_new <- droplevels(counts_summary$orf_category_new)

ggplot(data = counts_summary, aes(x = orf_category_new, y = mean, fill = orf_category_new, color = orf_category_new)) +
  geom_quasirandom(size = 0.1, color = "grey20", fill = "grey20", shape = 16) +
  scale_y_continuous(trans = "log2", labels = scales::comma) +
  coord_flip() +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", width = 0.6) +
  scale_x_discrete(limits = rev(levels(counts_summary$orf_category_new))) +
  # scale_fill_manual(values = setNames(colors[1:length(unique(df_row$orf_category_new))],
                                      # unique(df_row$orf_category_new))) +
  scale_fill_manual(values = colors) +
  labs(x = "Translated ORFs", y = "Normalized P-site Counts") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor.y = element_blank())

ggsave(
filename = "/hpc/pmc_vanheesch/projects/Da/Neuroblastoma_neoantigens/02_riboseq_analysis/ribo_nbl_merged/analysis/orf_quantification/translated_orfs_normalized_psitecounts.pdf", device = "pdf", path = , width = unit(4, "cm"), height = unit(3, "cm"))
