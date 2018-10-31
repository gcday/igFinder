import pandas as pd
import os
import sys
import shutil
import subprocess
from Bio import SeqIO

sys.stdout = open(str(snakemake.log), "w")

clones = pd.read_table(snakemake.input["clones"], sep="\t").set_index("cloneId", drop=False)
sig_clones = clones.query("cloneFraction >= 0.05 | cloneCount >= 10")["cloneId"]
print(sig_clones, flush = True)
if os.path.exists(snakemake.output.tempdir):
	shutil.rmtree(snakemake.output.tempdir)
os.makedirs(snakemake.output.tempdir, exist_ok = True)
vOpt_dir = os.path.join(snakemake.output.tempdir, "vOpt_temp")
os.makedirs(vOpt_dir, exist_ok = True)
reads_file = os.path.join(snakemake.output.tempdir, "reads.fastq")
vdj_seq_recs = []
master_dir = os.getcwd()
for clone in sig_clones:
	os.chdir(master_dir)
	clone_dir = os.path.join(vOpt_dir, str(clone))
	os.makedirs(clone_dir, exist_ok = True)
	subprocess.run(["mixcr", "exportReadsForClones", "-f", "-s", 
									str(snakemake.input.clna), "--id", str(clone), reads_file], stdout=sys.stdout, stderr=sys.stdout)
	r1_file = os.path.abspath(os.path.join(snakemake.output.tempdir, "reads_cln{}_R1.fastq".format(str(clone))))
	r2_file = os.path.abspath(os.path.join(snakemake.output.tempdir, "reads_cln{}_R2.fastq".format(str(clone))))
	vOpt_clone = os.path.abspath(os.path.join(snakemake.output.tempdir, "cln_{}_vOptOut".format(str(clone))))
	os.chdir(clone_dir)

	vOpt_cmd = ["VelvetOptimiser.pl", "-v", "-s 21", "-e 91", "-x 4", "-t {}".format(str(snakemake.threads)),
							"-c max", "-k max", "-a", "-d {}".format(vOpt_clone),  
							'-f "-fastq -shortPaired -separate {} {}"'.format(r1_file, r2_file)]
	print(" ".join(vOpt_cmd), flush = True)
	subprocess.run(" ".join(vOpt_cmd), shell=True, stdout=sys.stdout, stderr=sys.stdout)
	contigs_path = os.path.join(vOpt_clone, "contigs.fa")
	if os.path.exists(contigs_path):
		print("cln {} assembly successful".format(str(clone)), flush = True)
		contigs = list(SeqIO.parse(contigs_path, "fasta"))
		print("Found {} contigs".format(len(contigs)), flush = True)
		longest = contigs[0]
		cdr3_seq = clones.loc[(clone), ["nSeqCDR3"]].dropna()[0]
		mixcr_bestVGene = clones.loc[(clone), ["bestVGene"]].dropna()[0]
		cloneCount = clones.loc[(clone), ["cloneCount"]].dropna()[0]
		print(cdr3_seq, flush = True)
		if cdr3_seq in longest.seq or cdr3_seq in longest.seq.reverse_complement():
			print("Found CDR3 in expected contig.", flush = True)
			print(longest.id, flush = True)
			longest.id = "|".join([str(snakemake.wildcards.sample), "clone_{}".format(str(clone)), 
													 mixcr_bestVGene, "{} reads".format(cloneCount)])
			vdj_seq_recs.append(longest)
		else:
			print("CDR3 not found in expected contig.", flush = True)
			# "{}|clone_{}|{}".format(str(snakemake.wildcards.sample), str(clone), longest.id)	
	else:
		print("cln {} assembly failed".format(str(clone)), flush = True)
os.chdir(master_dir)
SeqIO.write(vdj_seq_recs, snakemake.output.contigs, "fasta")

		# shell("mixcr exportReadsForClones -f -s {snakemake.input.clna} --id {clone} {reads_file}")