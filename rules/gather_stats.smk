def sample_to_clones():
  return {os.path.abspath("data/mixcr/clone_summary/{0}_clone_summary.txt".format(s)) : s for s in sample_list}
def sample_to_read_counts():
  return {os.path.abspath("data/read_counts/{0}.txt".format(s)) : s for s in sample_list}

def sample_to_igblast():
  return {os.path.abspath("results/igblast/{0}_igblast_output.txt".format(s)) : s for s in sample_list}

def sample_to_igblast():
  return {s : os.path.abspath("results/igblast/{0}_igblast_output.txt".format(s)) for s in sample_list}

# rule gather_igblast:
#   input:
#     expand(os.path.abspath("results/igblast/{sample}_igblast_output.txt"), sample = sample_list)
#   output:
#     config["igblast_file"]
#   log:
#     os.path.abspath("logs/gather_stats/gather_igblast.log")
#   resources:
#     mem_mb=24000
#   run:
#     files = sample_to_igblast()
#     igblast_samples = samples.copy()
#     max_V_ident = []
#     for sample in samples["sample"]:
#       with open(files[sample]) as igblast:
#         max_ident = ""
#         for line in igblast.readlines():
#           if "Total\tN/A" in line:
#             split = line.strip().split()
#             curr_ident = float(split[-1])
#             max_ident = str(curr_ident)
#             if (curr_ident < 99.0):
#               break
#         max_V_ident += [max_ident]
#     igblast_samples["max_V_ident"] = max_V_ident
#     # igblast_samples = igblast_samples.drop(["sample", "mutect2.vcf"], axis = 1)
#     igblast_samples.to_csv(output[0], sep = "\t")
rule gather_igblast:
  input:
    expand(os.path.abspath("results/igblast/{sample}_igblast_output.txt"), sample = sample_list)
  output:
    config["igblast_file"]
  log:
    os.path.abspath("logs/gather_stats/gather_igblast.log")
  resources:
    mem_mb=24000
  run:
    import re
    files = sample_to_igblast()
    igblast_samples = samples.copy()
    max_light_ident = []
    max_heavy_ident = []
    light_read_counts = []
    heavy_read_counts = []
    for sample in samples["sample"]:
      with open(files[sample]) as igblast:
        max_ident, light_ident, heavy_ident , heavy_reads, light_reads = ("", "", "", "", "")
        for line in igblast.readlines():
          if "# Query: " in line:
            line = line.replace("# Query: ", "")
            split = line.strip().split("|")
            chain_type = "Heavy" if "IGH" in split[2] else "Light"
            # if (chain_type == "Heavy"):
            #   break
            # max_reads = split[-1]
            num_reads = re.search(r'(\d*\.?\d*) reads', split[-1])
            if num_reads:
              num_reads = num_reads.group(1)
            print(num_reads)

          if "Total\tN/A" in line:
            split = line.strip().split()
            curr_ident = float(split[-1])
            if (chain_type == "Light" and light_ident == ""):
              light_ident = str(curr_ident)
              light_reads = num_reads
              # if (light_ident < 99.0):
              #   break
            elif (chain_type == "Heavy" and heavy_ident == ""):
              heavy_ident = str(curr_ident)
              heavy_reads = num_reads
            # if (curr_ident < 99.0):
            #     break
        max_light_ident += [light_ident]
        light_read_counts += [light_reads]
        max_heavy_ident += [heavy_ident]
        heavy_read_counts += [heavy_reads]

    igblast_samples["max_light_ident"] = max_light_ident
    igblast_samples["max_heavy_ident"] = max_heavy_ident
    igblast_samples["light_read_counts"] = light_read_counts
    igblast_samples["heavy_read_counts"] = heavy_read_counts


    # igblast_samples = igblast_samples.drop(["sample", "mutect2.vcf"], axis = 1)
    igblast_samples.to_csv(output[0], sep = "\t")
    
rule gather_read_counts:
  input:
    sample_to_read_counts()
  output:
    "read_counts.tsv"
  resources:
    mem_mb=24000
  log:
    os.path.abspath("logs/gather_stats/gather_read_counts.log")
  run:
    files = sample_to_read_counts()
    with open(output[0], "w+") as summary:
      summary.write("sample\tread_count\n")
      for sample in files:
        with open(sample) as file:
          count = file.readline()
          summary.write(files[sample] + "\t" + count)
          # if len(clones) < 2:
            # summary.write(files[sample] + "\n")
            # continue
          # for line in clones[1:]:
            # summary.write(files[sample] + "\t" + line)

rule gather_output:
  input:
    sample_to_clones()
  output:
    config["clones_file"]
  resources:
    mem_mb=24000
  log:
    os.path.abspath("logs/gather_stats/gather_output.log")
  run:
    files = sample_to_clones()
    with open(output[0], "w+") as summary:
      first = True
      for sample in files:
        with open(sample) as file:
          clones = file.readlines()
          if (first):
            summary.write("sample\t" + clones[0])
            first = False
          if len(clones) < 2:
            summary.write(files[sample] + "\n")
            continue
          for line in clones[1:]:
            summary.write(files[sample] + "\t" + line)
            
rule calc_clonality:
  conda: "../envs/R.yaml"
  input:
    config["clones_file"],
    "read_counts.tsv"
  resources:
    mem_mb=24000
  log:
    os.path.abspath("logs/gather_stats/calc_clonality.log")
  output:
    config["stats_file"]    
  script:
    "../scripts/calc_clonality.R"

localrules: calc_clonality, gather_output, gather_read_counts, gather_igblast
