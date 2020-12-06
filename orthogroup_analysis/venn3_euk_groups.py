			# This script is to be run in main /OG_arb-fal dir, drawing on total proteome data of groups in /gr_outputs.
			# 3 sys.args (different eukaryote groups) should be specified, which will produce a 3-way venn diagram, saved as a pdf.
			# The venn diagram will return the total proteome overlap of the 3 eukaryote groups.
			# NOTE: this does not return a venn diagram for uniquely shared OGs!
			# NOTE: /gr_outputs only currently contains data for the original 11-way split (need to update!).

from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib import pyplot as plt
import group
import sys

path = "/mnt/c/Users/scamb/Documents/uob_msc/Genome_data/OG_arb-fal/gr_outputs/gr_"

# Creates sets of each group's total proteome (identified by OG numbers).
# sys.args must be valid eukaryote group names in dataset (include caps).
group1 = set(group.count_ogs(path + sys.argv[1] + ".txt"))
group2 = set(group.count_ogs(path + sys.argv[2] + ".txt"))
group3 = set(group.count_ogs(path + sys.argv[3] + ".txt"))

# Add _unweighted to venn3 function to balance the resulting circles (regardless of data).
venn3([group1, group2, group3], (sys.argv[1], sys.argv[2], sys.argv[3]))

plt.savefig("venn3_euk_groups.pdf")

