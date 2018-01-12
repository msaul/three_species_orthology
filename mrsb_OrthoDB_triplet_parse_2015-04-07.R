# 2015-04-07
# Michael C. Saul
# msaul [at] illinois.edu

# MRSB orthology analysis: OrthoDB
# Finding OrthoDB triplets using the big OrthoDB table

# Setting working directory
setwd("~/Desktop/mrsb/orthology/OrthoDB/")

# Reading in the table
# Note: OrthoDB v8 table can be downloaded from:
# ftp://cegg.unige.ch/OrthoDB8/Eukaryotes/Genes_to_OGs/ODB8_EukOGs_genes_ALL_levels.txt.gz
ODB.file = gzfile("ODB8_EukOGs_genes_ALL_levels.txt.gz","rt")
ODB = read.table(ODB.file,
                 header=T,
                 sep="\t",
                 stringsAsFactors=F,
                 quote="",
                 fill = T)
close(ODB.file)
rm(list = c("ODB.file"))

# Finding lines for each species of interest
beeLines = grep("Apis mellifera",ODB$organism)
mouseLines = grep("Mus musculus",ODB$organism)
sticklebackLines = grep("Gasterosteus aculeatus",ODB$organism)

# Getting tables for each species of interest
beeODB = ODB[beeLines,]
mouseODB = ODB[mouseLines,]
sticklebackODB = ODB[sticklebackLines,]

# Removing the rest of the ODB for size and lines objects for clarity
rm(list = c("ODB"))
rm(list = ls(pattern = "Lines"))

# Finding overlap between each species (gene level)
beeGenes = unique(beeODB$odb8_og_id)
mouseGenes = unique(mouseODB$odb8_og_id)
sticklebackGenes = unique(sticklebackODB$odb8_og_id)
threespeciesGenes = mouseGenes[which(mouseGenes %in% sticklebackGenes)]
threespeciesGenes = threespeciesGenes[which(threespeciesGenes %in% beeGenes)]

# Building a data frame for the OrthoDB triplets
threespeciesTriplets = data.frame(row.names = threespeciesGenes,
                                  odb8_og_id = threespeciesGenes,
                                  bee_protein_ids = rep(NA, times = length(threespeciesGenes)),
                                  mouse_protein_ids = rep(NA, times = length(threespeciesGenes)),
                                  stickleback_protein_ids = rep(NA, times = length(threespeciesGenes)))

# Looping through the three species and making a list of triplets
for (i in 1:length(threespeciesGenes)) {
  currentGene = threespeciesTriplets[i,"odb8_og_id"]
  threespeciesTriplets[i,"bee_protein_ids"] = paste(beeODB[which(beeODB$odb8_og_id == currentGene),"protein_id"], collapse = ";")
  threespeciesTriplets[i,"mouse_protein_ids"] = paste(mouseODB[which(mouseODB$odb8_og_id == currentGene),"protein_id"], collapse = ";")
  threespeciesTriplets[i,"stickleback_protein_ids"] = paste(sticklebackODB[which(sticklebackODB$odb8_og_id == currentGene),"protein_id"], collapse = ";")
}

# Defining proteins that have single definitions for each gene
beeSingle = grep(";", threespeciesTriplets$bee_protein_ids, invert = T)
mouseSingle = grep(";", threespeciesTriplets$mouse_protein_ids, invert = T)
sticklebackSingle = grep(";", threespeciesTriplets$stickleback_protein_ids, invert = T)
threespeciesSingle = sticklebackSingle[sticklebackSingle %in% mouseSingle]
threespeciesSingle = threespeciesSingle[threespeciesSingle %in% beeSingle]
threespeciesNotsingle = threespeciesTriplets[-threespeciesSingle,]
threespeciesSingle = threespeciesTriplets[threespeciesSingle,]
row.names(threespeciesSingle) = threespeciesSingle$mouse_protein_ids

# Using biomaRt to find mouse IDs for each single protein
stopifnot(require("biomaRt"))

# Setting biomaRt to generate annotations from Ensembl gene IDs (the row names)
maRt = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
               host = "may2012.archive.ensembl.org",
               dataset = "mmusculus_gene_ensembl")
maRt_filter = "ensembl_peptide_id"
maRt_attributes = c("mgi_symbol","mgi_description","chromosome_name","start_position",
                    "end_position","strand","band","ensembl_gene_id","ensembl_peptide_id")

# Grabbing biomart data to annotate the analyzed files
mm_mrsb_maRt_query = threespeciesSingle$mouse_protein_ids
mm_mrsb_biomaRt = getBM(maRt_attributes, maRt_filter, mm_mrsb_maRt_query, maRt)

# Resolving duplicates in biomaRt output
mm_mrsb_biomaRt_duplicates = mm_mrsb_biomaRt[duplicated(mm_mrsb_biomaRt$ensembl_gene_id),]
mm_mrsb_biomaRt = mm_mrsb_biomaRt[!duplicated(mm_mrsb_biomaRt$ensembl_gene_id),]
row.names(mm_mrsb_biomaRt) = mm_mrsb_biomaRt$ensembl_gene_id
mm_mrsb_biomaRt$duplicate = rep(FALSE,times=nrow(mm_mrsb_biomaRt))

if(nrow(mm_mrsb_biomaRt_duplicates) != 0) {
  for (i in 1:nrow(mm_mrsb_biomaRt_duplicates)) {
    current_ensembl_id = mm_mrsb_biomaRt_duplicates$ensembl_gene_id[i]
    for (k in 1:ncol(mm_mrsb_biomaRt_duplicates)) {
      current_attribute = colnames(mm_mrsb_biomaRt_duplicates)[k]
      current_value = mm_mrsb_biomaRt_duplicates[i,k]
      current_nonduplicate = mm_mrsb_biomaRt[current_ensembl_id,current_attribute]
      if (current_nonduplicate == current_value) {
        mm_mrsb_biomaRt[current_ensembl_id,current_attribute] = current_nonduplicate
      } else {
        mm_mrsb_biomaRt[current_ensembl_id,current_attribute] = paste(
          current_nonduplicate," / ",current_value,sep="")
        mm_mrsb_biomaRt[current_ensembl_id,"duplicate"] = TRUE
      }
    }
  }
} else {
  mm_mrsb_biomaRt = mm_mrsb_biomaRt
}
mm_mrsb_biomaRt = mm_mrsb_biomaRt[order(row.names(mm_mrsb_biomaRt)),]

# Adding the three species OrthoDB annotation into this annotation
mm_mrsb_biomaRt$orthodb_id = threespeciesSingle[mm_mrsb_biomaRt$ensembl_peptide_id,"odb8_og_id"]
mm_mrsb_biomaRt$bee_protein_id = threespeciesSingle[mm_mrsb_biomaRt$ensembl_peptide_id,"bee_protein_ids"]
mm_mrsb_biomaRt$stickleback_protein_id = threespeciesSingle[mm_mrsb_biomaRt$ensembl_peptide_id,"stickleback_protein_ids"]

# Running stickleback maRt
sb_maRt = maRt = useMart(biomart = "ensembl",
                         dataset = "gaculeatus_gene_ensembl")
sb_maRt_attributes = c("ensembl_gene_id","ensembl_peptide_id")
sb_mrsb_biomaRt = getBM(sb_maRt_attributes, maRt_filter, mm_mrsb_biomaRt$stickleback_protein_id, sb_maRt)
row.names(sb_mrsb_biomaRt) = sb_mrsb_biomaRt$ensembl_peptide_id

# Adding bee and stickleback gene ID
mm_mrsb_biomaRt$bee_gene_id = gsub("-PA", "", mm_mrsb_biomaRt$bee_protein_id)
mm_mrsb_biomaRt$stickleback_gene_id = sb_mrsb_biomaRt[mm_mrsb_biomaRt$stickleback_protein_id,"ensembl_gene_id"]

write.table(mm_mrsb_biomaRt,
            "three_species_OrthoDB_triplets.tsv",
            row.names=F,
            quote=F,
            sep="\t")

# Saving data into a .Rdata file and grabbing an md5 checksum
stopifnot(require("tools"))
save(list = ls(), file = "three_species_OrthoDB_parse_2015-04-07.Rdata")
md5sum("three_species_OrthoDB_parse_2015-04-07.Rdata")

# Writing session info for long-term reproducibility
sessionInfo()
