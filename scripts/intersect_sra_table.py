import pandas as pd

samples_path = "/farmshare/user_data/gday/mayo/data/CoMMpass/CoMMpass_clin_info.csv"
sra_table_path = "/farmshare/user_data/gday/mayo/data/CoMMpass/20181105_SraRunTable.txt"
out_samples_path = "/farmshare/user_data/gday/mayo/projects/igFinder_CoMMpass/CoMMpass_clin_WXS_avail.csv"

samples = pd.read_table(samples_path, sep = ",").set_index("sample", drop=False)
sra_table = pd.read_table(sra_table_path, sep = "\t").set_index("Sample_Name", drop = False)

sra_table = sra_table.query('Assay_Type == "WXS"')

WXS_samples = []
srr_list = []
for sample in samples["sample"]:
	sra_match = sra_table[sra_table['Sample_Name'].str.contains(str(sample))]
	if (len(sra_match) == 1):
		# print(len(sra_match))
		WXS_samples += [sample]
		srr_list += list(sra_match["Run"])[:1]
	elif (len(sra_match) > 1):
		print("More than one match!")

WXS_table = samples[[sample in WXS_samples for sample in samples["sample"]]]
# WXS_table.loc[:, "SRR"] = srr_list
WXS_table = WXS_table.assign(SRR=srr_list)
WXS_table.drop(["sample", "Unnamed: 0"], axis = 1).to_csv(out_samples_path)

# WXS_table[""]

# WXS_samples = [sample for sample in samples["sample"] if (len(sra_table[sra_table['Sample_Name'].str.contains(str(sample))]) >= 1)]

# WXS_table = samples[[sample in WXS_samples for sample in samples["sample"]]]
# sample_list = samples_of_int["sample"]
# WXS_table.drop(["sample", "Unnamed: 0"], axis = 1).to_csv(out_samples_path)
# WXS_table.to_csv(out_samples_path)


