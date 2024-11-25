# Exploring the temporal dynamics of transcriptional responses to high day and night temperatures in maize
## SNK for upstream transcriptome analsis

```
# Upstream transcriptome analysis using snakemake. Output: 06_transcripts_quant; 05_stringtie_merged/stringtie_merged.gtf
snakemake -j 00_Snakemake_rna-seq.smk

# 00_extract_exp.py was used to extract FPKM values for all genes of different samples
python3 00_extract_exp.py <input_path> <output_path>

# 00_ID_transition.py was used for the gene_id modification
python3 00_ID_transition.py <merge_gtf> <input_path> <output_path>

# 00_rm_MGSRT.py was used for remove unpair gene
python3 00_rm_MGSRT.py <input_paht> <output_path>
```

## DESeq2 for analyzing differential genes
Scripts in deseq2 folder
## Co-expression network building
Scripts in Co-expression folder
## GO and TF family enrichment analysis
Scripts in enrichment folder
