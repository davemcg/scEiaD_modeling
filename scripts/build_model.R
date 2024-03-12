# stage 1
## pull data from one species 
## filter to age group (dev or mature)
## select random (up to 50k) cells per study and output those barcodes
## output barcodes for the non-selected cells
## run scvi on biowulf
library(tidyverse)
sample_meta <- data.table::fread('~/git/scEiaD_quant/sample_meta.scEiaD_v1.2024_02_28.01.tsv.gz')
cell_meta <- data.table::fread('~/data/scEiaD_2024_02_28/hs111.adata.solo.2024_03_07.obs.csv.gz')[,-1] %>% 
  relocate(barcode) %>% 
  filter(solo_doublet == "FALSE")

hs111_eye <- cell_meta %>% 
  filter(study_accession != 'SRP362101') %>% 
  mutate(stage = case_when(as.numeric(age) <= 10 ~ 'Developing', 
                           TRUE ~ 'Mature'), 
         side = case_when(tissue == 'Brain Choroid Plexus' ~ 'Brain Choroid Plexus',
                          grepl("Choroid|RPE", tissue) ~ 'BackEye',
                          grepl("Retina", tissue) ~ 'BackEye',
                          grepl("Outf", tissue) ~ 'FrontEye',
                          grepl("Iris", tissue) ~ 'FrontEye',
                          grepl("Sclera", tissue) ~ 'FrontEye',
                          grepl("Cornea", tissue) ~ 'FrontEye',
                          grepl("Macula", tissue) ~ 'BackEye',
                          grepl("Trabecul", tissue) ~ 'FrontEye',
                          grepl("Optic", tissue) ~ 'BackEye',
                          TRUE ~ tissue)) %>% 
  filter(organism == 'Homo sapiens',
         !grepl("^#", sample_accession),
         source == 'Tissue',
         #tissue %in% c("Macula", "Retina"),
         side %in% c("FrontEye", "BackEye"),
         #capture_type == 'cell',
         kb_tech %in% c("10xv1","10xv2","10xv3"),
         stage == 'Mature') 

hs111_study <- hs111_eye %>% 
  group_by(study_accession) %>% count() %>% filter(n>1000)
set.seed(101294)
hs111_ref_bcs <- hs111_eye %>% 
  filter(study_accession %in% hs111_study$study_accession) %>% 
  group_by(study_accession) %>% 
  sample_n(50000, replace = TRUE) %>% 
  unique() %>% 
  pull(barcode)

hs111_query_bcs <- hs111_eye %>% 
  filter(!barcode %in% hs111_ref_bcs) %>% 
  pull(barcode)

hs111_ref_bcs %>% write(gzfile('~/git/scEiaD_modeling/data/hs111_mature_alleye_ref_bcs.csv.gz'))
hs111_query_bcs %>% write(gzfile('~/git/scEiaD_modeling/data/hs111_mature_alleye_query_bcs.csv.gz'))

# now go to biowulf2:/data/OGVFB_BG/scEiaD/2024_02_28
# mamba deactivate; mamba activate scvi1.0.4
# python ~/git/scEiaD_modeling/workflow/scripts/run_scvi.py \
# hs111.adata.solo.2024_03_07.h5ad    \
# /home/mcgaugheyd/git/scEiaD_modeling/data/hs111_mature_alleye_ref_bcs.csv.gz   hs111_mature_alleye_4k   scviModel.hs111_mature_alleye_4k   hs111_mature_alleye_cells4k.h5ad   hs111_mature_alleye_cells4k.obs.csv.gz   --n_top_genes 4000 --hvg_span 1 --mode barcode

obs <- data.table::fread('~/data/scEiaD_2024_02_28/hs111_mature_alleye_cells4k_50e.obs.csv.gz')[,-1] %>% relocate(barcode) %>% 
  mutate(side = case_when(tissue == 'Brain Choroid Plexus' ~ 'Brain Choroid Plexus',
                          grepl("Choroid|RPE", tissue) ~ 'BackEye',
                          grepl("Retina", tissue) ~ 'BackEye',
                          grepl("Outf", tissue) ~ 'FrontEye',
                          grepl("Iris", tissue) ~ 'FrontEye',
                          grepl("Sclera", tissue) ~ 'FrontEye',
                          grepl("Cornea", tissue) ~ 'FrontEye',
                          grepl("Macula", tissue) ~ 'BackEye',
                          grepl("Trabecul", tissue) ~ 'FrontEye',
                          grepl("Optic", tissue) ~ 'BackEye',
                          TRUE ~ tissue)) 

obs_labels <- obs %>% group_by(leiden, MajorCellType) %>% 
  filter(MajorCellType != '') %>% 
  summarise(Count = n()) %>% mutate(Ratio = Count/sum(Count)) %>% 
  filter(Ratio > 0.10) %>% 
  arrange(-Ratio) %>% 
  summarise(MCT = paste0(MajorCellType, collapse = ', '), 
            CT = gsub(", .*","",MCT)) 
obs_leiden_ct <- obs %>% group_by(leiden) %>% summarise(lCount = n())


obs %>% 
  left_join(obs_labels, by = 'leiden') %>% 
  ggplot(aes(x=umap1,y=umap2)) +
  scattermore::geom_scattermore(aes(color = MCT), pointsize = 0.8, alpha = 0.5) +
  ggrepel::geom_label_repel(data = . %>% group_by(CT) %>% 
                              summarise(umap1 = mean(umap1),
                                        umap2 = mean(umap2)),
                            aes(label = CT)) +
  scale_color_manual(values = c(pals::alphabet2(), pals::glasbey()) %>% unname()) + 
  cowplot::theme_cowplot() + theme(legend.position = "none")


obs %>% 
  left_join(obs_labels, by = 'leiden') %>% 
  ggplot(aes(x=umap1,y=umap2)) +
  scattermore::geom_scattermore(aes(color = MCT)) +
  facet_wrap(~CT) +
  scale_color_manual(values = c(pals::alphabet2(), pals::glasbey()) %>% unname()) + 
  cowplot::theme_cowplot()

set.seed(234)
hs111_ref_bcs_r2 <- obs %>% 
  left_join(obs_labels, by = 'leiden') %>% 
  select(barcode, MCT, study_accession) %>% 
  group_by(MCT, study_accession) %>% 
  sample_n(1000, replace = TRUE) %>% 
  unique() %>% 
  ungroup() %>% 
  group_by(MCT) %>% 
  sample_n(5000, replace = TRUE) %>% 
  unique()
study_count <- hs111_ref_bcs_r2 %>% group_by(study_accession) %>% 
  summarise(Count = n()) %>% filter(Count > 1000)
hs111_ref_bcs_r2 <- hs111_ref_bcs_r2 %>% filter(study_accession %in% study_count$study_accession) %>% 
  pull(barcode)
hs111_query_bcs_r2 <- obs %>% 
  filter(!barcode %in% hs111_ref_bcs_r2) %>% 
  pull(barcode)


hs111_ref_bcs_r2 %>% write(gzfile('~/git/scEiaD_modeling/data/hs111_mature_alleyeR2_ref_bcs.csv.gz'))
hs111_query_bcs_r2 %>% write(gzfile('~/git/scEiaD_modeling/data/hs111_mature_alleyeR2_query_bcs.csv.gz'))

obs_r2 <- data.table::fread('~/data/scEiaD_2024_02_28/hs111_mature_alleye_cells2k_50e_30l_R2.obs.csv.gz')[,-1] %>% relocate(barcode) %>% 
  mutate(side = case_when(tissue == 'Brain Choroid Plexus' ~ 'Brain Choroid Plexus',
                          grepl("Choroid|RPE", tissue) ~ 'BackEye',
                          grepl("Retina", tissue) ~ 'BackEye',
                          grepl("Outf", tissue) ~ 'FrontEye',
                          grepl("Iris", tissue) ~ 'FrontEye',
                          grepl("Sclera", tissue) ~ 'FrontEye',
                          grepl("Cornea", tissue) ~ 'FrontEye',
                          grepl("Macula", tissue) ~ 'BackEye',
                          grepl("Trabecul", tissue) ~ 'FrontEye',
                          grepl("Optic", tissue) ~ 'BackEye',
                          TRUE ~ tissue)) 

obs_labels_r2 <- obs_r2 %>% group_by(leiden, MajorCellType) %>% 
  filter(MajorCellType != '') %>% 
  summarise(Count = n()) %>% mutate(Ratio = Count/sum(Count)) %>% 
  filter(Ratio > 0.10) %>% 
  arrange(-Ratio) %>% 
  summarise(MCT = paste0(MajorCellType, collapse = ', '), 
            CT = gsub(", .*","",MCT)) 
obs_leiden_ct_r2 <- obs_r2 %>% group_by(leiden) %>% summarise(lCount = n())


p_r2 <- obs_r2 %>% 
  left_join(obs_labels_r2, by = 'leiden') %>% 
  filter(barcode %in% hs111_ref_bcs_r2) %>% 
  ggplot(aes(x=umap1,y=umap2)) +
  scattermore::geom_scattermore(aes(color = CellType), pointsize = 4.8, alpha = 0.5) +
  scale_color_manual(values = c(pals::alphabet2(), pals::glasbey(), pals::brewer.set1(n=8)) %>% unname()) + 
  cowplot::theme_cowplot() + theme(legend.position = "none")

p_r2 +  ggrepel::geom_label_repel(data = . %>% group_by(CT) %>% 
                                     summarise(umap1 = mean(umap1),
                                               umap2 = mean(umap2)),
                                   aes(label = CT)) 

p_r2 + facet_wrap(~MCT)
