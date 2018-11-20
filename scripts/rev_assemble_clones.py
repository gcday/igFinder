import pandas as pd
import os
import sys
import shutil
import subprocess
from Bio import SeqIO

class assembler(object):
	def __init__(snakemake):
		self.tmpdir = snakemake.output.tempdir
		if os.path.exists(self.tmpdir):
			shutil.rmtree(self.tmpdir)
		self.vOpt_dir = os.path.join(self.tmpdir, "vOpt_temp")
		os.makedirs(vOpt_dir, exist_ok = True)
		self.reads_file = os.path.join(self.tmpdir, "reads.fastq")
		self.master_dir = os.getcwd()
		self.vdj_seq_recs = []
		self.clones = pd.read_table(snakemake.input["clones"], sep="\t").set_index("cloneId", drop=False)
		(self.min_frac, 
		 self.min_reads, 
		 self.max_clones) = (snakemake.config["assembly"]["min_frac_to_assemble"],
												 snakemake.config["assembly"]["min_reads"], 
												 snakemake.config["assembly"]["max_clones_to_assemble"])


	def get_clones_to_assemble():
		query_str = "cloneFraction >= {} & cloneCount >= {}".format(self.min_frac, self.min_reads)
		self.sig_clones = clones.query(query_str)["cloneId"][:self.max_clones]
		print(self.sig_clones, flush = True)

	def assemble_clone(clone):
		os.chdir(self.master_dir)
		clone_dir = os.path.join(self.vOpt_dir, str(clone))
		os.makedirs(clone_dir, exist_ok = True)
		subprocess.run(["mixcr", "exportReadsForClones", "-f", "-s", 
										str(snakemake.input.clna), "--id", str(clone), self.reads_file], 
										stdout=sys.stdout, stderr=sys.stdout)
		r1_file = os.path.join(self.tmpdir, "reads_cln{}_R1.fastq".format(str(clone)))
		r2_file = os.path.join(self.tmpdir, "reads_cln{}_R2.fastq".format(str(clone)))
		vOpt_clone = os.path.join(self.tmpdir, "cln_{}_vOptOut".format(str(clone)))
		os.chdir(clone_dir)
		vOpt_cmd = ["VelvetOptimiser.pl", "-v", "-s 21", "-e 91", "-x 4", "-t {}".format(str(snakemake.threads)),
								"-c max", "-k max", "-a", "-d {}".format(vOpt_clone),
								'-f "-fastq -shortPaired -separate {} {}"'.format(r1_file, r2_file)]
		print(" ".join(vOpt_cmd), flush = True)
		subprocess.run(" ".join(vOpt_cmd), shell=True, stdout=sys.stdout, stderr=sys.stdout)
		self.parse_vOpt_output(clone, os.path.join(vOpt_clone, "contigs.fa"))
			
	def parse_vOpt_output(clone, contigs_path):
		if not os.path.exists(contigs_path):
			print("cln {} assembly failed".format(str(clone)), flush = True)
			return
		print("cln {} assembly successful".format(str(clone)), flush = True)
		cdr3_seq = self.clones.loc[(clone), ["nSeqCDR3"]].dropna()[0]
		mixcr_bestVGene = self.clones.loc[(clone), ["bestVGene"]].dropna()[0]
		cloneCount = self.clones.loc[(clone), ["cloneCount"]].dropna()[0]
		contigs = list(SeqIO.parse(contigs_path, "fasta"))
		print("Found {} contigs".format(len(contigs)), flush = True)
		len_longest = 0
		longest = -1
		for contig in contigs:
			if cdr3_seq in contig.seq or cdr3_seq in contig.seq.reverse_complement():
				if len(contig.seq) > len_longest:
					longest = contig
					len_longest = len(longest.seq)
		if longest == -1:
			print("CDR3 not found in contigs.", flush = True)
		else:
			print("Found CDR3 in contig {}.".format(longest.id), flush = True)
			longest.id = "|".join([str(snakemake.wildcards.sample), 
														"clone_{}".format(str(clone)), 
													 	mixcr_bestVGene, "{} reads".format(cloneCount)])
			self.vdj_seq_recs.append(longest)
	def assmb_clones():
		self.get_clones_to_assemble()
		for clone in self.sig_clones:
			self.assemble_clone(clone)
		os.chdir(self.master_dir)
		SeqIO.write(self.vdj_seq_recs, snakemake.output.contigs, "fasta")


sys.stdout = open(str(snakemake.log), "w")
assmb = assembler(snakemake)
assmb.assmb_clones()