				# The group module contains all necessary functions involved in the OG search pipeline.
				# Former code map and group list functions are now consolidated within OrthogroupSearch.
import re
import sys
import glob

class OrthogroupSearch:

	""" OrthogroupSearch requires a code map argument upon creation. This can be either "codes", "codes_alt" or "codes_alt_18",
	and determines how eukaryote groups in the dataset are delimited. Subsequent functions are derived from the selected code map,
	including find_group, the primary function used to search for uniquely shared OG between different eukaryote groups. """

	def __init__(self, codes):
		self.codes = codes

	def code_map(self):

		""" code_map returns a dictionary containing sp.codes as keys and eukaryote groups as values.
		This can be used to search through genomic data for presence/absence of particular groups. """

		cm = {}
		codes = self.codes

		try:
			with open("/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/OG_arb-fal/Eukaryote_" + codes + ".txt") as f:
				for line in f:
					fields = re.split("\t", line.strip())
					cm[fields[0]] = fields[1]

		except FileNotFoundError as e:
			print(f'{e} not found. Please use "codes", "codes_alt" or "codes_alt_18" to select the necessary code map.')

		return cm

	def group_list(self):

		""" group_list uses the dictionary defined in code_map() to create a list containing all eukaryote groups in the dataset.
		Note the group list returned always corresponds to the equivalent code_map since it is derived from the same source file. """

		group_names = []
		spcode = self.code_map()
		for group in spcode.values():
			if group not in group_names:
				group_names.append(group)

		return group_names


	def long_name_codes(self):

		""" long_name_codes is an alternative to the other code maps, returning full sp. names instead of just abbreviations.
		The code map in this function is derived from code_map. """

		code_map = self.code_map()
		long_name_code_map = {}

		for sp_code in code_map:
			group = code_map[sp_code]

			# Note: this file is currently used as a source for the names, but in principle any file containing all the long names could be used.
			with open("/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/iqtree/Ancyromonads_Collodictyonids/Ancyromonads_Collodictyonids.constr") as f:
				for line in f:
					res = re.search(sp_code + r"(_\w+-\w+)", line.strip())
					if res:
						long_name = res.group(1)
						full_name = sp_code + long_name
						long_name_code_map[full_name] = group

		return long_name_code_map


	def find_group(self, query):

		""" find_group takes a list of eukaryote groups as an argument, to be compared with the groups present in each .fal file in directory (/OG_arb-fal).
		Matching files are written to an output file with the list elements (i.e. group names) as the file title.
		The code_map is selected when an OrthogroupSearch object is created. This will affect how eukaryote groups are defined.

		It is important here that the data type of query is a list - this will be converted as appropriate within the function to match groups_present.
		The initial if/else loop ensures the data type is correct (using x). """

		code_map = self.code_map()
		to_parse = glob.glob("/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/OG_arb-fal/*.fal")

		# Join set into a sensible filename.
		# x necessary as a conditional in place of isinstance(query, list) - see loop further down.
		if isinstance(query, list):
			filename = "_".join(query)
			x = 0

		# [] necessary to match the str query to groups_present (which is a single element list).
		else:
			filename = query
			query = [query]
			x = 1
		genome = open(filename + "_output.txt", "w")
		genomewrite = genome.write(f'{query}\n')
		k = 0
		for file in to_parse:
			groups_present = []
			with open(file) as f:
				for line in f:
					if line.startswith('>'):
						fields = re.split('_', line)
						species_code = fields[0][1:]

                        # Map sp. code to group.
						for i in code_map:
							eugroup = code_map[species_code]

                            # Add to group list if new group.
							if eugroup not in groups_present:
								groups_present.append(eugroup)

            # Turns multiple element lists into sets.
			if x == 0:
				set_present = set(groups_present)
				query = set(query)

            # if isinstance(query, list) not used since this would convert it to a set after first iteration.
            # Single elements (i.e own group comparisons) must remain as they are.
			else:
				set_present = groups_present
			if query == set_present:
				file = re.search(r'(OG\d{7}.fal$)', file)

            	# Adds file to output if it matches query.
                # k is a counter for total number of OG matches.
				genomewrite = genome.write(f'{file.group(1)}\n')
				k += 1
		genomewrite = genome.write(f'Shared gene families: {k}')
		genome.close()


	def sp_total_dic(self):

		""" sp_total_dic iterates over every file in /sp_individual_total, and stores the total number of OGs that each sp. has in a dictionary. """

		codes = self.code_map()
		to_parse = glob.glob("/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/OG_arb-fal/sp_individual_total/*.txt")
		sp_list = []
		sp_total_dic = {}

		# Extract sp. codes from code map.
		for sp in codes:
			sp_list.append(sp)
		for sp in sp_list:
			for file in to_parse:

				# Match sp to the relevant file.
				# Append sp. code keys and total values to the sp_total_dic.
				sp_file = re.search(r"(\w{8})_total.txt$", file, re.IGNORECASE)
				if sp_file.group(1) == sp:
					total = parse_OG(file)
					sp_total_dic[sp] = total

		return sp_total_dic


def parse_OG(file):

	""" parse_OG is designed to parse out the 'Shared gene families: ' data from find_group() output files. """

	with open(file) as f:
		for line in f:

			# Match summary of OGs shared i.e. "Shared gene families: NNNN".
			og = re.search(r"\w*:\s(\w*)", line)
			if og:
				return og.group(1)


def count_ogs(file):

	""" count_ogs searches through any document and appends to a list all the different OG strings.
	This regex is not specified with any preceding "_" or trailing ".fal" due to the variety of contexts OGs may occur in. """

	og_list = []
	with open(file) as f:
		for line in f:
			res = re.search(r"(OG\d{7})", line.strip())
			if res:
				og = res.group(1)
				if og not in og_list:
					og_list.append(og)
	return og_list

def parse_total():

	""" parse_total returns a dictionary with the total number of OGs that each group has. This extracts totals data from /total_genome (not via group_total.txt as it did formerly). """

	# All group divisions are present as *_allOGs.txt files.
	to_parse = glob.glob("/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/OG_arb-fal/total_genome/*.txt")
	groups = []
	totals = []

	for file in to_parse:
		name = re.search(r"([A-Z]\w+)_allOGs.txt$", file)
		og = parse_OG(file)
		groups.append(name.group(1))
		totals.append(og)
	total_dict = dict(zip(groups, totals))

	return total_dict

