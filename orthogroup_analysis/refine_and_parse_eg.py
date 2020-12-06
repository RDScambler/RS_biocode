			# Script for improving the readability of eggnog output files.
			# In principle should be applicable to all outputs, providing sp labels are suffixed with "_OGxxxxxxx".
			# Important not to discard raw .emapper.annotations files since GOs, kegg etc. are not transferred.
			# Script parses out all OGs that have no hits in eggnog to a separate file.
			# Script now also parses bacterial and non-bacterial OGs using a customisable cutoff point.
import re
import group

# Ask user for query taxa.
# Append to a list all the uniquely shared OGs between the query taxa using count_ogs.
while True:
	try:
		query = input("Enter query taxa (e.g. Metamonads_Discoba): \n")
		preanno_full_list = group.count_ogs("/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/OG_arb-fal/new_outputs/" + query + "_output.txt")
		if preanno_full_list is not None:
			break
	except FileNotFoundError as e:
		print(f'{e} not found. Please enter a valid query.')
		
# Define code map.
ortho = group.OrthogroupSearch('codes_alt_18')
output = open(query + "_refined_eg.txt", "w")
og_list = []

# regex used to split up the first column of eggnog output file.
sp_label_regex = r"(\w{,8})_\w*[tn|gn]\d+_(\w+\-\w+)_(OG\d+)"					

with open(query + ".emapper.annotations") as f:	
	for line in f:
		res = re.search(sp_label_regex, line)
		if res:
			sp_code = res.group(1)
			sp_name = res.group(2)					
			og = res.group(3)

			# Eukaryote group is derived from the codes dict.
			eugroup = ortho.code_map()[sp_code]

			# .ljust formats for appearance.						
			fnames = eugroup.ljust(20) + sp_name.ljust(30)				
			tabline = re.split("\t", line.strip())
			e = tabline[2].ljust(15)

			# Remove whitespace from between tax_hit strings (for future parsing).
			ftax_hit = tabline[4].replace(" ", "_").ljust(40)			
			gene_id = tabline[5].ljust(10)
			tax_scope = tabline[17].ljust(20)
			func_cat = tabline[-2]
			text_description = tabline[-1]

			# Append ogs to og_list for the purposes of identifying new ogs.
			if og not in og_list:							
				og_list.append(og)
				outputWrite = output.write("\n")
			outputWrite = output.write(f"{og}\t{fnames}{e}{ftax_hit}{tax_scope}{gene_id}{func_cat}\t{text_description}\n")

output.close()

# Here a cutoff point is requested, this can be customised as required.
# The proportion defines how many bacterial sequences can make up the eukaryote/archaea output file.
# I.e what level is being considered as homologous rather than a contaminant.
while True:
	cutoff = float(input("What is the bacterial cutoff point for orthogroups (as a proportion of the group's sequences)?\n"))
	if cutoff < 1 and cutoff > 0:
		break
	else:
		print("Enter a proportion for cutoff (between 0 and 1).")

# ean = Eukaryotes, Archaea, no hits.
ean_output = open(query + "_ean_" + str(cutoff * 100) + "percent.txt", "w")
bac_output = open(query + "_bac_" + str(cutoff * 100) + "percent.txt", "w")

euk_arc_list = []
bac_list = []

for og in og_list:
	j = 0
	tax_list = []

	# While loop necessary to iterate over the file twice.
	# The first loop appends all taxa in og to tax_list.
	# The second compares the frequency of Bacterial seqs to the cutoff.
	while j < 2:												
		i = 0												
		with open(query + "_refined_eg.txt") as f:							
			if j < 1:
				for line in f:

					# Match each og in turn.
					if og in line:								
						spline = re.split(r"\s+", line.strip(), maxsplit=6)
						tax = spline[5]
						tax_list.append(tax)
			else:
				for taxon in tax_list:
					if taxon == "Bacteria":
						i += 1

				# If bacteria count is below cutoff, ogs are written to the 'ean' output file.		
				if i < cutoff * len(tax_list):							
					for line in f:
						if og in line:
							outputWrite = ean_output.write(line)

							# For the purposes of collecting summary stats.
							if og not in euk_arc_list:				
								euk_arc_list.append(og)

				# Otherwise they are written to the bacterial output file.				
				else:										
					for line in f:
						if og in line:
							outputWrite = bac_output.write(line)
							if og not in bac_list:
								bac_list.append(og)
		
		# j ensures while loop occurs twice.
		j += 1												

bac_output.close()

# Identify ogs not analysed by eggnog.
# Access and parse out relevant sp and og data from pre-annotation file.
# non-hits are appended to the ean_output file.
# It is necessary to parse data from "xxxxxx_toannotate.fal' file since sp. name and code data are here.
for og in preanno_full_list:
	if og not in og_list:
		with open("/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/eggnog-mapper/eg_" + query + "_toannotate.fal") as f:
			for line in f:
				if og in line:
					res = re.search(sp_label_regex, line)
					if res:
						sp_code = res.group(1)
						sp_name = res.group(2)
						eugroup = ortho.code_map()[sp_code]
						fnames = eugroup.ljust(20) + sp_name
						no_hits_outputWrite = ean_output.write(f"{og}\t{fnames}\n")

ean_output.close()

output = open(query + "_refined_eg.txt", "a")

# Summary stats.
# Will be written out to "xxxxxx_refined_eg.txt'.
bac_prop = round(len(bac_list) / len(og_list) * 100, 2)
overall_prop = round(len(og_list) / len(preanno_full_list) * 100, 2)
stats = "__STATS__"

per_bac = "\nPercentage of bacterial-dominated OGs at " + str(cutoff) + " cutoff: " + str(bac_prop)
per_euk = "Percentage of Eukaryote/Archaeal OGs: " + str(100 - bac_prop)
per_overall = "Overall percentage of OGs with hits in eggNOG: " + str(overall_prop) + "\n"

total = "Total number of OGs in query: " + str(len(preanno_full_list))
tot_minus = "Total number of OGs after bacteria removed: " + str(len(preanno_full_list) - len(bac_list))
tot_euk = "Total number of OGs with significant non-bacterial hits: " + str(len(euk_arc_list))

print(per_bac)
print(per_euk)
print(per_overall)
print(total)
print(tot_minus)
print(tot_euk)

output.write(f"\n{stats}{per_bac}\n{per_euk}\n{per_overall}\n{total}\n{tot_minus}\n{tot_euk}")

output.close()
