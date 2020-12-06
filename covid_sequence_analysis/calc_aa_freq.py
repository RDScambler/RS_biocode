            # Code to determine significance of AA frequency differences along the alignment of two sets of sequences.
            # Looking for COVID spike protein seqs with variable sites - separated into two groups by collection date.
            # Differences in AA frequency at two different periods may signify immunologically relevant mutations.
            
            # All european spike data now sourced from two files - can be parsed as necessary. 
            # Total records in Europe: 637.
            # Optional: enter sys.argv[1] to filter sequences by region.
            # Regions: 'Finland', 'Italy', 'Sweden', 'Spain', 'Germany', 'France', 'Greece', 'Serbia', 'Czech Republic', 
            # 'Poland', 'Russia', 'Netherlands', 'United Kingdom', 'Malta', 'Romania'.

            # This analysis was inspired by:
            # Korber et al. (2020) Tracking Changes in SARS-CoV-2 Spike: Evidence that D614G Increases Infectivity of the COVID-19 Virus.
            

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import re
from collections import Counter
import scipy.stats as stats
import sys
import json
from check_date import check_date_format



# Obtain date boundaries for two spike seq sets.
# Use check_date_format to check for and return correct date format.
# Note: some dates do not include a specific day (YYYY-MM). In these instances they are considered the earliest point in the month.
date_boundaries = []
print('Enter two pairs of dates for spike sequence collection (YYYY-MM-DD). Earliest collection date = 2020-01-28.')
while True:
    answer = input()
    if check_date_format(answer)[0] == True:
        date_boundaries.append(check_date_format(answer)[1])
    else:
        print(check_date_format(answer)[1])
    if len(date_boundaries) == 4:
        break
        

# Access data file containing all SARS-COV-2 spike protein data isolated from human hosts in Europe.
# Also access annotations file (cannot access necessary data from seq object unfortunately).
data = '/mnt/c/Users/scamb/Documents/Programming/covid/spike_europe_data/ncbi_dataset/data/protein.gpff'
annotation_file = '/mnt/c/Users/scamb/Documents/Programming/covid/spike_europe_data/ncbi_dataset/data/data_report.jsonl'

# Generate list of annotation data.
# json.loads is ideal for data in this format.
annotation_list = []
for line in open(annotation_file, 'r'):
    annotation_list.append(json.loads(line))

# Create separate record lists to count number of seqs in each list.
date_1 = []
date_2 = []

# Loop over records to parse depending on desired collection date and location of seqs.
# Access the relevant location data for each accession in the .json file (mapped with the ID in the corresponding .gpff file).
# Access the relevant date data likewise. Note this must also be done via the .json file since the given date in 
# record.annotations['date'] is NOT the actual collection date.
file = 'spike_seqs_to_align.fasta'
with open(file, 'w') as f:
    records = SeqIO.parse(data, 'genbank')
    for record in records:
        for annotation in annotation_list:

            # Use continue to avoid KeyError.
            if 'annotation' not in annotation.keys():
                continue 
            elif record.id == annotation['annotation']['genes'][0]['cds'][0]['protein']['accessionVersion']:
                loc = annotation['location']['geographicLocation']
                date = annotation['isolate']['collectionDate']
     
                # Filter by region and date.
                # Write to file to be aligned.
                # Append IDs to date lists in order to check seq group sizes (i.e. length of date lists).
                if len(sys.argv) > 1:
                    if sys.argv[1] in loc and date_boundaries[0] < date < date_boundaries[1]:
                        SeqIO.write(record, f, 'fasta')
                        date_1.append(record.id)
                    elif sys.argv[1] in loc and date_boundaries[2] < date < date_boundaries[3]:
                        SeqIO.write(record, f, 'fasta')
                        date_2.append(record.id)

                # If no region specified in sys.argv[1], filter only by date.
                elif date_boundaries[0] < date < date_boundaries[1] or date_boundaries[2] < date < date_boundaries[3]:
                    if date_boundaries[0] < date < date_boundaries[1]:
                        date_1.append(record.id)
                    elif date_boundaries[2] < date < date_boundaries[3]:
                        date_2.append(record.id)
                            
                    SeqIO.write(record, f, 'fasta')


# Print length of each sequence group for user.
# Exit program if either group has no seqs.
group_lists = [date_1, date_2]
for i, ls in enumerate(group_lists):
    print('Sequences in group %i: %i' % (i + 1, len(ls)))
    if len(ls) == 0:
        print('Group %i dates did not return any sequences. Try expanding date boundaries (or omit location from sys.arg) for more results.' % (i + 1))
        sys.exit()

# Align sequences.
cline = MuscleCommandline(input = 'spike_seqs_to_align.fasta', out = 'spike_seqs.aln', clw=True)
cline()

# Read in alignment.
# Get alignment length for iterating (len(alignment) returns the number of seqs (i.e. rows not columns)).
alignment = AlignIO.read('spike_seqs.aln', 'clustal')
aln_length = alignment.get_alignment_length()

# Loop over seq positions.
for a in range(0, aln_length):
    pos = alignment[:, a]

    # Ignore Xs (unknowns) and indels - these typically result in empty group lists.
    # Is preferable to focus only on positions that are unambiguous.
    if '-' in pos or 'X' in pos:
        continue

    # Create empty lists to fill with AAs.
    pre_onset = []
    post_onset = []

    # Append AAs to onset list depending on collection date.
    for i, s in enumerate(alignment):
        if s.id in date_1:
            pre_onset.append(pos[i])
        elif s.id in date_2:
            post_onset.append(pos[i])

    # Identify the mode of the first list.
    mode = max(pre_onset, key = pre_onset.count)
    mode_2 = max(post_onset, key = post_onset.count)

    # Create a dict of AA frequencies for the list using Counter.
    dict_aa_freqs = Counter(pre_onset)
    
    # Check it is the dominant AA (e.g. >= 50%).
    if dict_aa_freqs[mode] >= 0.5 * len(pre_onset):
        percent_dominance = round(dict_aa_freqs[mode] / len(pre_onset) * 100, 3)
    
        # Specify relevant values for Fisher's exact test.
        # Create a new dict for the second list and index using the mode from the first.
        dict_aa_freqs_2 = Counter(post_onset)

        # Obtain the frequency of pre onset dominant AA in post onset seqs.
        # Obtain the non-dominant remainder (for FET).
        observed = dict_aa_freqs_2[mode]
        non_target_observed = len(post_onset) - observed

        # Obtain inverse values for FET.
        target_not_observed = dict_aa_freqs[mode]
        remainder_not_observed = len(pre_onset) - target_not_observed

        # Calculate the percent which dominant AA from first list makes up the second.
        percent_comprised_2 = round(observed / len(post_onset) * 100, 3)

        # Compute the probability using FET to determine whether the frequency difference is significant.
        # Note two-sided test is default.
        odds, p_value = stats.fisher_exact([[observed, non_target_observed], [target_not_observed, remainder_not_observed]])

        # Print results with AA percentages for context if FET returns significant result.
        if p_value <= 0.05:
            print('Pre onset dominant AA: %s (percent dominance: %f). Post onset is %f percent %s at this position (%i). FET: %f' % 
            (mode, percent_dominance, percent_comprised_2, mode, a, p_value))
            
            # Determine if the post onset group has a different dominant AA.
            if dict_aa_freqs_2[mode_2] >= 0.5 * len(post_onset):
                percent_dominance_2 = round(dict_aa_freqs_2[mode_2] / len(post_onset) * 100, 3)
                print('%s is the dominant post onset AA at this position (%i). Percent dominance: %f.' % (mode_2, a, percent_dominance_2))

# Spike is highly conserved. Can I expand search to entire COVID genome?
# Should be fairly trivial in theory - all record IDs and annotation map the same (I think).