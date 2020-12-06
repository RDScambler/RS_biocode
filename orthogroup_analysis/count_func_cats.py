			# This file contains the count_func_cats function, used to count all the functional categories present in an eggnog output file.
import group
import re
import sys 


def count_func_cats(file):
	
	""" count_func_cats takes an eggnog output file as an argument, and counts different functional categories of each OG 
	(note only one category is recorded per OG). These frequncies are stored in a functional category frequency dictionary. 
	Note the function executes succesfully regardless of whether the eggnog output has been refined (i.e. parsed out) or is raw. """
	
	og_list = group.count_ogs(file)
	overall_list = []
	for og in og_list:
		with open(file) as f:
			
			# The cat_list must be appended for each OG in turn.
			cat_list = []
			for line in f:
				if og in line:
					
					# Each seq is assigned up to 3 cats, so re must account for this.
					res = re.search(r"\s([A-Z]{1,3})\s", line.strip())
					if res:
						func_cat = res.group(1)
					
						# Iterates over each letter (for cases where multiple cats are assigned).
						for cat in range(len(func_cat)):
							if func_cat[cat] not in cat_list:
								cat_list.append(func_cat[cat])
			
			# Ignore eggnog nonhits.
			if len(cat_list) == 0:
				pass
			else:
				for cat in range(len(cat_list)):
					overall_list.append(cat_list[cat])

	# List comprehension is used to generate a dict of counts for each cat in overall_list.
	func_cat_freq_dict = {i:overall_list.count(i) for i in overall_list}

	return func_cat_freq_dict

# Configure file at command line.
file = sys.argv[1]
print(count_func_cats(file))
