import pandas as pd

samples_path = "/farmshare/user_data/gday/mayo/data/CoMMpass/CoMMpass_clin_info.csv"
sra_table_path = "/farmshare/user_data/gday/mayo/data/CoMMpass/20181105_SraRunTable.txt"
out_samples_path = "/farmshare/user_data/gday/mayo/projects/igFinder_CoMM_RNAseq/CoMMpass_clin_RNAseq_avail.csv"


samples = pd.read_table(samples_path, sep = ",").set_index("sample", drop=False)
sra_table = pd.read_table(sra_table_path, sep = "\t").set_index("Sample_Name", drop = False)

sra_table = sra_table.query('Assay_Type == "RNA-Seq"')

RNAseq_samples = [sample for sample in samples["sample"] if (len(sra_table[sra_table['Sample_Name'].str.contains(str(sample))]) >= 1)]

RNAseq_table = samples[[sample in RNAseq_samples for sample in samples["sample"]]]
# sample_list = samples_of_int["sample"]
RNAseq_table.drop(["sample", "Unnamed: 0"], axis = 1, inplace=True)
RNAseq_table.to_csv(out_samples_path)


