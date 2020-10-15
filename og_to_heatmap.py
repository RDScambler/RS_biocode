			# Unique OG data is parsed from /new_outputs.
			# Add a second sys.argv (any) to convert data to proportional data.
			# Output data is converted into a heatmap, with own group values removed.

import re
import group
import glob
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors
from shift_colormap import shiftedColorMap


# Define code map and group_list with sys.argv[1]
og = group.OrthogroupSearch(sys.argv[1])
group_list = og.group_list()
group_list = sorted(group_list)
data = []

# Extract unique OG data from /new_outputs dir.
for name in group_list:
	to_parse = glob.glob("/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/OG_arb-fal/new_outputs/*.txt")
	for file in to_parse:
		filename = re.search(r"([A-Z]\w+)_([A-Z]\w+)_output.txt$", file)
		if filename:
			# Skips files in new_outputs (if certain groups are being excluded).
			# Must be done for group(1) and group(2) since name order in filename can vary.
			if filename.group(1) not in group_list:
				pass
			elif filename.group(2) not in group_list:
				pass
			elif filename.group(1) == name:
				shared_og = group.parse_OG(file)
				data.append(shared_og)
			elif filename.group(2) == name:
				shared_og = group.parse_OG(file)
				data.append(shared_og)


# Create an array from the data.
data_array = np.array(data, dtype = float)

# Adjust the array shape to match the length of group list.
data_array.shape = (len(group_list), len(group_list))


def convert_to_proportional(data_point, euk_group):

	""" convert_to_proportional takes OG data individually and uses eukaryote genome total data from parse_total to
	calculate the percentage of that genome the OGs constitute. """

	totals = group.parse_total()
	total = totals[euk_group]
	prop_data_point = round(int(data_point) / int(total) * 100, 2)

	return prop_data_point


# Convert unique OG data to Dataframe.
# DataFrame is necessary for indexing -> removal of own group values.
df = pd.DataFrame(data_array, index = group_list, columns = group_list)

# Remove own group values (these ought to be NaNs, as opposed to 0s.)
for groups in df:
	for index_group in df.index:
		if groups == index_group:
			df.loc[groups, index_group] = np.nan
		# Enter a 2nd sys.argv in order to convert raw data to proportional.
		elif len(sys.argv) > 2:
			df.loc[groups, index_group] = convert_to_proportional(df.loc[groups, index_group], groups)


# Mask NaNs by converting back to array.
df = df.to_numpy()
np.ma.masked_invalid(df)


### Create heatmap ###

fig, ax = plt.subplots(figsize=(20, 10))

# Choosing colour from list as below allows customisation of the colourmap.
# Use midpoint to adjust the centre of the colormap.
orig_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(" ", ["white", "gold", "red"])
newcmap = shiftedColorMap(orig_cmap, midpoint=0.2)
im = ax.imshow(df, cmap=newcmap)

# Show all ticks.
ax.set_xticks(np.arange(len(group_list)))
ax.set_yticks(np.arange(len(group_list)))
ax.set_xticklabels(group_list, fontsize = 13)
ax.set_yticklabels(group_list, fontsize = 13)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Create colorbar.
cbarlabel = "Number of orthogroups"
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.tick_params(labelsize=15)
cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom", labelpad = 30, fontsize = 15)

# Loop over data dimensions and create text annotations.
# This conditional loop seems to be the only way to remove text from the main diagonal and .0 from the other cells simultaneously.
# Converting df to int is the solution - not necessary when using prop_data, hence the conditional sys.argv length.
for i in range(len(group_list)):
	for j in range(len(group_list)):
		if i == j:
			text = ax.text(i, j, "")
		elif len(sys.argv) > 2:
			text = ax.text(j, i, df[i, j],ha="center", va="center", color="black", fontsize = 11)
		else:
			text = ax.text(j, i, int(df[i, j]),ha="center", va="center", color="black", fontsize = 11)

ax.set_title("Uniquely shared orthogroups", fontsize = 22, fontweight = "bold", pad = 30, loc = 'center')
fig.tight_layout()
plt.savefig("og_heatmap_minusown.pdf")

