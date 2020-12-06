			# setquery.py can be used to search for OGs shared uniquely between any combination of groups, input by the user.
			# OrthogroupSearch requires an argument specifying the code map to be used. The corresponding group_list is then generated.
			# The particular group_list must be borne in mind when inputting groups.

import group
import sys

# sys.argv[1] must equal codes, codes_alt, or codes_alt_18.
og = group.OrthogroupSearch(sys.argv[1])
group_list = og.group_list()

query = []

# Obtain a list of eukaryote groups from the user.
print("Enter eukaryote group names (press enter when done):")
while True:
	answer = input()
	if answer == '':
		break
	elif answer not in group_list:
		print("Incorrect name. Please check spelling.")
	else:
		query.append(answer)


# Must not convert query into a set - this disrupts find_group function (data type should be list).
sorted_group = sorted(query)
og.find_group(sorted_group)

