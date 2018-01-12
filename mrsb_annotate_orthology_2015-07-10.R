# 2015-07-10
# Michael C. Saul
# msaul [at] illinois.edu

# MRSB orthology annotation

load("~/Desktop/mrsb/orthology/OrthoDB/three_species_OrthoDB_parse_2015-04-07.Rdata")

stopifnot(require("biomaRt"))

# Setting biomaRt to generate annotations from Ensembl gene IDs (the row names)
mm_maRt = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                  host = "may2012.archive.ensembl.org",
                  dataset = "mmusculus_gene_ensembl")
mm_maRt_filter = "ensembl_peptide_id"
mm_maRt_attributes = c("ensembl_peptide_id","ensembl_gene_id")

mm_mrsb_maRt_query = unlist(strsplit(threespeciesTriplets$mouse_protein_ids, ";"))

# Grabbing biomart data to annotate the analyzed files
mm_mrsb_biomaRt = getBM(mm_maRt_attributes, mm_maRt_filter, mm_mrsb_maRt_query, mm_maRt)
mm_peptides_not_in_old_annotation = mm_mrsb_maRt_query[which(!mm_mrsb_maRt_query %in% mm_mrsb_biomaRt$ensembl_peptide_id)]

# Setting biomaRt to generate annotations from new Ensembl peptide IDs (the row names)
mm_new_maRt = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
               host = "feb2014.archive.ensembl.org",
               dataset = "mmusculus_gene_ensembl")
mm_mrsb_new_biomaRt = getBM(mm_maRt_attributes, mm_maRt_filter, mm_peptides_not_in_old_annotation, mm_new_maRt)
mm_gene_annotations = rbind(mm_mrsb_biomaRt, mm_mrsb_new_biomaRt)

sb_maRt = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  host = "feb2014.archive.ensembl.org",
                  dataset = "gaculeatus_gene_ensembl")
sb_maRt_filter = "ensembl_peptide_id"
sb_maRt_attributes = c("ensembl_peptide_id","ensembl_gene_id")

sb_mrsb_maRt_query = unlist(strsplit(threespeciesTriplets$stickleback_protein_ids, ";"))
sb_mrsb_new_biomaRt = getBM(sb_maRt_attributes, sb_maRt_filter, sb_mrsb_maRt_query, sb_maRt)
sb_gene_annotations = sb_mrsb_new_biomaRt

threespeciesTriplets$bee_gene_ids = as.character(gsub("-PA","",threespeciesTriplets$bee_protein_ids))
threespeciesTriplets$mouse_gene_ids = as.character(threespeciesTriplets$mouse_protein_ids)
threespeciesTriplets$stickleback_gene_ids = as.character(threespeciesTriplets$stickleback_protein_ids)

for (i in 1:nrow(mm_gene_annotations)) {
  current_gene = mm_gene_annotations[i,"ensembl_gene_id"]
  current_peptide = mm_gene_annotations[i,"ensembl_peptide_id"]
  threespeciesTriplets$mouse_gene_ids = gsub(current_peptide, current_gene, threespeciesTriplets$mouse_gene_ids)
}

for (i in 1:nrow(sb_gene_annotations)) {
  current_gene = sb_gene_annotations[i,"ensembl_gene_id"]
  current_peptide = sb_gene_annotations[i,"ensembl_peptide_id"]
  threespeciesTriplets$stickleback_gene_ids = gsub(current_peptide, current_gene, threespeciesTriplets$stickleback_gene_ids)
}

orthology_edges = as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(orthology_edges) = c("orthology_group","species_1","species_2","gene_species_1","gene_species_2")
for (i in 1:nrow(threespeciesTriplets)) {
  current_orthology_group = threespeciesTriplets[i,"odb8_og_id"]
  current_bee_genes = unlist(strsplit(threespeciesTriplets[i,"bee_gene_ids"],";"))
  current_mouse_genes = unlist(strsplit(threespeciesTriplets[i,"mouse_gene_ids"],";"))
  current_stickleback_genes = unlist(strsplit(threespeciesTriplets[i,"stickleback_gene_ids"],";"))
  l_bee = length(current_bee_genes)
  l_mouse = length(current_mouse_genes)
  l_stickleback = length(current_stickleback_genes)
  for (j in 1:l_bee) {
    for (k in 1:l_mouse) {
      current_df = data.frame(orthology_group = current_orthology_group,
                              species_1 = "honeybee",
                              species_2 = "mouse",
                              gene_species_1 = current_bee_genes[j],
                              gene_species_2 = current_mouse_genes[k])
      orthology_edges = rbind(orthology_edges, current_df)
    }
  }
  for (j in 1:l_bee) {
    for (k in 1:l_stickleback) {
      current_df = data.frame(orthology_group = current_orthology_group,
                              species_1 = "honeybee",
                              species_2 = "stickleback",
                              gene_species_1 = current_bee_genes[j],
                              gene_species_2 = current_stickleback_genes[k])
      orthology_edges = rbind(orthology_edges, current_df)
    }
  }
  for (j in 1:l_mouse) {
    for (k in 1:l_stickleback) {
      current_df = data.frame(orthology_group = current_orthology_group,
                              species_1 = "mouse",
                              species_2 = "stickleback",
                              gene_species_1 = current_mouse_genes[j],
                              gene_species_2 = current_stickleback_genes[k])
      orthology_edges = rbind(orthology_edges, current_df)
    }
  }
}