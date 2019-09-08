# Activating Qiime2
source activate qiime2-2018.11

#16S
#################################################################################################

# Directories
Database_fasta_16S='/Metagenetics_Projects/Databases/SILVA/SILVA_132_QIIME_release/silva_132_99_18S.fna'
Database_txt_16S='/home/olivierl/apps/Metagenetics_Projects/Databases/SILVA/SILVA_132_QIIME_release/18S_taxonomy_7_levels.txt'
Outputs_16S='16S'
Feature_ASVs_table_16S='/16S/dada_table.txt' 
Feature_ASVs_repseq_16S='/16S/rep_ASVs.fasta' 

# Training a Naive Bayes classifier - assigning taxonomy from different databases
# Importing the fasta and taxonomy files
printf "\nImporting reference database\n" 
qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path $Database_fasta_16S \
   --output-path $Outputs_16S/database_seq_16S.qza
 
 qiime tools import \
   --type 'FeatureData[Taxonomy]' \
   --input-format HeaderlessTSVTaxonomyFormat \
   --input-path $Database_txt_16S \
   --output-path $Outputs_16S/database_taxonomy_16S.qza
 
# Extracting the sequence portion we have targeted
 printf "\nExtracting the sequence portion we have targeted\n" 
 qiime feature-classifier extract-reads \
   --i-sequences $Outputs_16S/database_seq_16S.qza \
   --p-f-primer CCTACGGGNGGCWGCAG \
   --p-r-primer GACTACHVGGGTATCTAATCC \
   --o-reads $Outputs_16S/database_seq_16S.qza
   
# Trainning the classifier
printf "\nTrainning the classifier\n" 
qiime feature-classifier fit-classifier-naive-bayes \
   --i-reference-reads $Outputs_16S/database_seq_16S.qza \
   --i-reference-Outputs_16S $Outputs_16S/database_taxonomy_16S.qza \
   --o-classifier $Outputs_16S/taxo_classifier_16S.qza \
   --verbose

# Importing features table and representative sequences
qiime tools import \
  --input-path $Feature_ASVs_repseq_16S \
  --output-path $Outputs_16S/rep_seqs_dada2_16S.qza \
  --type 'FeatureData[Sequence]'

biom convert -i $Feature_ASVs_table_16S  -o $Outputs_16S/feature_ASVs_16S.biom --to-hdf5 --table-type="OTU table"
qiime tools import \
  --input-path $Outputs_16S/feature_ASVs_16S.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path $Outputs_16S/table_dada2_16S.qza

# Assigning taxonomy with the trainned classifier and looking at the results
printf "\nAssigning taxonomy with the trainned classifier\n" 

qiime feature-classifier classify-sklearn \
  --i-classifier $Outputs_16S/taxo_classifier_16S.qza \
  --i-reads $Outputs_16S/rep_seqs_dada2_16S.qza \
  --o-classification $Outputs_16S/taxonomy_16S.qza

# Exporting table
 printf "\nExporting tables\n" 
 qiime tools export \
   --input-path $Outputs_16S/table_dada2_16S.qza \
   --output-path $Outputs_16S/Data_analysis

 biom convert -i $Outputs_16S/Data_analysis/feature-table_16S.biom  -o $Outputs_16S/Data_analysis/feature-table_16S.txt --to-tsv

qiime tools export \
  --input-path $Outputs_16S/taxonomy_16S.qza \
  --output-path $Outputs_16S/Data_analysis

# Sequence clustering at 97% similarity and exporting the data
printf "\nClustering at 97% similarity with vsearch\n" 
qiime vsearch cluster-features-de-novo \
  --i-table $Outputs_16S/table_dada2_16S.qza \
  --i-sequences $Outputs_16S/rep_seqs_dada2_16S.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table $Outputs_16S/OTU_table_dada2_16S.qza \
  --o-clustered-sequences $Outputs_16S/rep_seqs_dada2_16S_97.qza \
  --verbose

printf "\nExporting data\n" 
qiime tools export \
  --input-path $Outputs_16S/OTU_table_dada2_16S.qza \
  --output-path $Outputs_16S/Data_analysis/OTU_table_dada2_16S.biom

biom convert -i $Outputs_16S/Data_analysis/OTU_table_dada2_16S.biom  -o $Outputs_16S/Data_analysis/OTU_table_dada2_16S.txt --to-tsv

qiime tools export \
  --input-path $Outputs_16S/rep_seqs_dada2_16S_97.qza \
  --output-path $Outputs_16S/rep_seqs_dada2_16S_97.fasta



#18S
#################################################################################################

# Directories
Database_fasta_18S='/Metagenetics_Projects/Databases/SILVA/SILVA_132_QIIME_release/silva_132_99_18S.fna'
Database_txt_18S='/Metagenetics_Projects/Databases/SILVA/SILVA_132_QIIME_release/18S_taxonomy_7_levels.txt'
Outputs_18S='18S'
Feature_ASVs_18S='/18S/dada_table.txt' 
Feature_ASVs_table_18S='/18S/dada_table.txt' 
Feature_ASVs_repseq_18S='/18S/rep_ASVs.fasta' 

# Training a Naive Bayes classifier - assigning taxonomy from different databases
# Importing the fasta and taxonomy files
printf "\nImporting reference database\n" 
qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path $Database_fasta_18S \
   --output-path $Outputs_18S/database_seq_18S.qza
 
 qiime tools import \
   --type 'FeatureData[Taxonomy]' \
   --input-format HeaderlessTSVTaxonomyFormat \
   --input-path $Database_txt_18S \
   --output-path $Outputs_18S/database_taxonomy_18S.qza
 
# Extracting the sequence portion we have targeted
 printf "\nExtracting the sequence portion we have targeted\n" 
 qiime feature-classifier extract-reads \
   --i-sequences $Outputs_18S/database_seq_18S.qza \
   --p-f-primer AGGGCAAKYCTGGTGCCAGC \
   --p-r-primer GRCGGTATCTRATCGYCTT \
   --o-reads $Outputs_18S/database_seq_18S.qza
   
# Trainning the classifier
printf "\nTrainning the classifier\n" 
qiime feature-classifier fit-classifier-naive-bayes \
   --i-reference-reads $Outputs_18S/database_seq_18S.qza \
   --i-reference-Outputs_18S $Outputs_18S/database_taxonomy_18S.qza \
   --o-classifier $Outputs_18S/taxo_classifier_18S.qza \
   --verbose

# Importing features table and representative sequences
qiime tools import \
  --input-path $Feature_ASVs_repseq_18S \
  --output-path $Outputs_18S/rep_seqs_dada2_18S.qza \
  --type 'FeatureData[Sequence]'

biom convert -i $Feature_ASVs_table_18S  -o $Outputs_18S/feature_ASVs_18S.biom --to-hdf5 --table-type="OTU table"
qiime tools import \
  --input-path $Outputs_18S/feature_ASVs_18S.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path $Outputs_18S/table_dada2_18S.qza

# Assigning taxonomy with the trainned classifier and looking at the results
printf "\nAssigning taxonomy with the trainned classifier\n" 

qiime feature-classifier classify-sklearn \
  --i-classifier $Outputs_18S/taxo_classifier_18S.qza \
  --i-reads $Outputs_18S/rep_seqs_dada2_18S.qza \
  --o-classification $Outputs_18S/taxonomy_18S.qza

# Exporting data
printf "\nExporting tables\n" 
qiime tools export \
  --input-path $Outputs_18S/table_dada2_18S.qza \
  --output-path $Outputs_18S/Data_analysis

biom convert -i $Outputs_18S/Data_analysis/feature-table_18S.biom  -o $Outputs_18S/Data_analysis/feature-table_18S.txt --to-tsv

qiime tools export \
  --input-path $Outputs_18S/taxonomy_18S.qza \
  --output-path $Outputs_18S/Data_analysis

# Sequence clustering at 97% similarity and exporting the data
printf "\nClustering at 97% similarity with vsearch\n" 
qiime vsearch cluster-features-de-novo \
  --i-table $Outputs_18S/table_dada2_18S.qza \
  --i-sequences $Outputs_18S/rep_seqs_dada2_18S.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table $Outputs_18S/OTU_table_dada2_18S.qza \
  --o-clustered-sequences $Outputs_18S/rep_seqs_dada2_18S_97.qza \
  --verbose

printf "\nExporting data\n" 
qiime tools export \
  --input-path $Outputs_18S/OTU_table_dada2_18S.qza \
  --output-path $Outputs_18S/Data_analysis/OTU_table_dada2_18S.biom

biom convert -i $Outputs_18S/Data_analysis/OTU_table_dada2_18S.biom  -o $Outputs_18S/Data_analysis/OTU_table_dada2_18S.txt --to-tsv

qiime tools export \
  --input-path $Outputs_18S/rep_seqs_dada2_18S_97.qza \
  --output-path $Outputs_18S/rep_seqs_dada2_18S_97.fasta
