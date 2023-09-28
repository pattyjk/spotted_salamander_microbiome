## Sequence processing in QIIME2
```
source activate qiime2-2023.5
cd  /hpcstor6/scratch01/p/patrick.kearns/Spot_sal_fungi/

#load raw FASTQ reads into QIIME
qiime tools import --type EMPSingleEndSequences --input-path ./data --output-path spot_sal_seqs.qza

#demultiplex reads
qiime demux emp-single \
  --i-seqs spot_sal_seqs.qza \
 --m-barcodes-file ITS_FINAL_RUN_Mapping_File.txt \
 --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences spot_sal_demux.qza \
  --o-error-correction-details  spot_sal_demux-details.qza \
  --p-no-golay-error-correction \
  
#quality filer
qiime quality-filter q-score \
--i-demux  spot_sal_demux.qza \
--o-filtered-sequences  spot_sal_demux-filtered.qza \
--o-filter-stats  spot_sal_demux-filter-stats.qza \

 #export filter stats
  qiime tools export --input-path   spot_sal_demux-filter-stats.qza --output-path filt_stats
 
  #call ASVs with deblur
  qiime deblur denoise-16S \
  --i-demultiplexed-seqs  spot_sal_demux-filtered.qza \
  --p-trim-length 120 \
  --o-representative-sequences  spot_sal_rep-seqs-deblur.qza \
  --o-table  spot_sal_table-deblur.qza \
   --o-stats  spot_sal_deblur-stats.qza
 
 #make phylogenetic tree with fasttree
 qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences spot_sal_rep-seqs-deblur.qza \
  --output-dir phylogeny-align-to-tree-mafft-fasttree
  
  #export deblur stats
  qiime tools export --input-path spot_sal_deblur-stats.qza --output-path deblur_stats
  
  #export rep seqs
  qiime tools export --input-path spot_sal_rep-seqs-deblur.qza --output-path rep_seqs
  
#export tree as NWK format
qiime tools export --input-path phylogeny-align-to-tree-mafft-fasttree/tree.qza --output-path tree
 
#pull trained taxonomy dataset from UNITE (https://unite.ut.ee/repository.php)
#https://john-quensen.com/tutorials/training-the-qiime2-classifier-with-unite-its-reference-sequences/
#wget https://files.plutof.ut.ee/public/orig/98/AE/98AE96C6593FC9C52D1C46B96C2D9064291F4DBA625EF189FEC1CCAFCF4A1691.gz
#tar xzf 98AE96C6593FC9C52D1C46B96C2D9064291F4DBA625EF189FEC1CCAFCF4A1691.gz
#rm 98AE96C6593FC9C52D1C46B96C2D9064291F4DBA625EF189FEC1CCAFCF4A1691.gz
cd cd sh_qiime_release_04.02.2020
#awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' sh_refs_qiime_ver8_99_04.02.2020.fasta | tr -d ' ' > sh_refs_qiime_ver8_99_04.02.2020_dev_uppercase.fasta

#qiime tools import \
#--type FeatureData[Sequence] \
#--input-path sh_refs_qiime_ver8_99_04.02.2020_dev_uppercase.fasta \
#--output-path unite-ver8-seqs_99_04.02.2020.qza

#qiime tools import \
#--type FeatureData[Taxonomy] \
--input-path sh_taxonomy_qiime_ver8_99_04.02.2020.txt \
--output-path unite-ver8-taxonomy_99_04.02.2020.qza \
--input-format HeaderlessTSVTaxonomyFormat

#train it (if needed)
#qiime feature-classifier fit-classifier-naive-bayes \
#--i-reference-reads  unite-ver8-seqs_99_04.02.2020.qza \
#--i-reference-taxonomy  unite-ver8-taxonomy_99_04.02.2020.qza \
#--o-classifier unite-ver8-99-classifier-04.02.2020.qza

cd ..
#assign taxonomy to  SILVA with sklearn
qiime feature-classifier classify-sklearn   --i-classifier sh_qiime_release_04.02.2020/unite-ver8-99-classifier-04.02.2020.qza   --i-reads spot_sal_rep-seqs-deblur.qza   --o-classification spot_sal_taxonomy.qza

#summarize taxonomy
qiime taxa barplot \
  --i-table spot_sal_table-deblur.qza \
  --i-taxonomy spot_sal_taxonomy.qza \
  --m-metadata-file ITS_FINAL_RUN_Mapping_File.txt \
  --o-visualization taxa-bar-plots.qzv
  ```
 #export summary of taxonomy
 qiime tools export --input-path taxa-bar-plots.qzv --output-path taxa_sum
 
 #export ASV table
 qiime tools export --input-path spot_sal_table-deblur.qza --output-path asv_table
 biom convert -i asv_table/table.biom -o asv_table.txt --to-tsv
