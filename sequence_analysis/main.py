import re
import glob

all_species = {}

# Creates a list containing the directory of each sequence file
file_list = glob.glob("INSERT SEQUENCE DIRECTORY")

for file_name in file_list:
    # Reading in the study
    infile = open(file_name)
    alltext = infile.read()

    # Extracting the sequences
    seqs = re.split('>.+?\n', alltext)
    del seqs[0]
    # Removing whitespace from sequences
    for i in range(len(seqs)):
        seqs[i] = re.sub('\s', '', seqs[i])
    
    # Replacing all character that don't correspond to a base
    seqs = [sequence.replace('R','-').replace('Y','-').replace('S','-')\
            .replace('K','-').replace('M','-').replace('N','-')\
            .replace('W','-') for sequence in seqs]
    
    # Code used to skip the study if there is a deletion / addition (optional based on desired results)
    indel_test = iter(seqs)
    indel_next_len = len(next(indel_test))
    if not all(len(l) == indel_next_len for l in indel_test):
        continue
    # Excluding studies with one individual, can be altered to filter out studies with a low sample size
    if len(seqs) < 2:
        continue
    
    # Extracting the species name
    species_name = re.split('\\\\', file_name)
    species_name = re.split('_', species_name[-1])[:2]
    species_name = '_'.join(species_name)

    # Adds species to the dictionary
    if species_name not in all_species:
        all_species[species_name] = {'total_sequences': 0, 'total_positions': 0, 'nucleotide_diversity_sum': 0,\
                                     'r2_sum_within_50': 0, 'total_within_50': 0,\
                                     'r2_sum_within_500': 0, 'total_within_500': 0}
    
    # Finding and storing the length of the longest sequence of the study
    longest_sequence = len(max(seqs, key = len))
    all_species[species_name]['total_positions'] += longest_sequence
        

    # Creating a nested dictionary to store the abundances of each base
    bases = {}
    for i in range(longest_sequence):
        bases[i] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    
    # Iterating through each position in each sequence, and tallying the abundance of each base
    for sequence in seqs:
        for position in range(len(sequence)):

            if sequence[position] == 'A':
                bases[position]['A'] += 1
            if sequence[position] == 'T':
                bases[position]['T'] += 1
            if sequence[position] == 'G':
                bases[position]['G'] += 1
            if sequence[position] == 'C':
                bases[position]['C'] += 1
        
      
    # Identifying which sites have a polymorphism
    polymorphism_sites = []
    for position in bases:
        # Calculating how many different bases are found at the site
        diversity = bool(bases[position]['A']) + bool(bases[position]['T']) \
            + bool(bases[position]['G']) + bool(bases[position]['C'])
        if diversity > 1:
            polymorphism_sites.append(position)
    #Excluding studies with no polymorphisms
    if len(polymorphism_sites) == 0:
        continue
    
    # Variable used to store the total nucleotide diversity before normalising for sequence length
    nd_total = 0
    # Calculating nucleotide diversity for the study
    for position in bases:
        # Creating new dictionary which doesn't include bases with a value of 0
        bases_present = (({k:v for k, v in bases[position].items() if v}))
        total_bases = (sum(bases_present.values()))
        # Calculating nucleotide diversity for the position
        if total_bases > 0:
            f = (min(bases_present.values())) / total_bases
            nd_total += 2*f*(1-f)
    
    # Calculating the study nucleotide diversity mean, and adding it to the dictionary
    all_species[species_name]['nucleotide_diversity_sum'] += nd_total / longest_sequence

    # Calculating the mean linkage disequilibrium for the study
    ld_stats = {}
    # Selecting the first position in the haplotype 
    for first_position in polymorphism_sites[:-1]:
        # Selecting the second position in the haplotype
        for second_position in polymorphism_sites[polymorphism_sites.index(first_position)+1:]:
           distance = second_position - first_position
           if distance not in ld_stats:
               ld_stats[distance] = {'pairs': 0, 'r2_sum': 0, 'r2_mean': 0}
           haplotypes = {}
           # Comparing the two sites across all sequences in the study
           for sequence in seqs:
                # Ensures that there isn't a deletion at the positions looked at in the sequence
                if sequence[first_position] != '-' and sequence[second_position] != '-':
                    # Constructing the current haplotype, and adding it to the dictionary
                    current_haplotype = sequence[first_position] + sequence[second_position]
                    if current_haplotype not in haplotypes:
                        haplotypes[current_haplotype] = 0
                    haplotypes[current_haplotype] += 1
           
           # Comparing observed haplotypes to expected haplotypes SHOULDNT ITERATE THROUGH EACH HAPLO
           observed_haplotype = (list(haplotypes.keys()))
           # Calculating the occurance frequency of the bases in the haplotype
           first_frequency = bases[first_position][observed_haplotype[0][0]] / (sum(bases[first_position].values()))
           second_frequency = bases[second_position][observed_haplotype[0][1]] / (sum(bases[second_position].values()))
           
           # Calculating the r2 value
           D = (haplotypes[observed_haplotype[0]] / (sum(haplotypes.values())) - (first_frequency * second_frequency))
           r2 = D**2 / (first_frequency * (1 - first_frequency) * second_frequency * (1 - second_frequency))
           # Recording the calculated variables
           ld_stats[distance]['pairs'] += 1   
           ld_stats[distance]['r2_sum'] += r2
           
    # Calculating the mean r2
    total_within_50 = 0
    r2_sum_within_50 = 0
    total_within_500 = 0
    r2_sum_within_500 = 0
    for distance in ld_stats:
        ld_stats[distance]['r2_mean'] = ld_stats[distance]['r2_sum'] / ld_stats[distance]['pairs']
        # Extracting data for pairs of positions within 50 bases
        if distance <= 50:
            all_species[species_name]['total_within_50'] += ld_stats[distance]['pairs']
            all_species[species_name]['r2_sum_within_50'] += ld_stats[distance]['r2_sum']
        # Extracting data for pairs of positions within 500 bases
        if distance <= 500:
            all_species[species_name]['total_within_500'] += ld_stats[distance]['pairs']
            all_species[species_name]['r2_sum_within_500'] += ld_stats[distance]['r2_sum']
    

# Recording the data
outfile = open('LD_and_ND_Results.txt','w')
outfile.write('species' + '\t' + 'mean_r2_within_50' + '\t' + 'mean_r2_within_500' + '\t' + 'mean_nd' + '\n')
for species in all_species:
    if (all_species[species]['total_within_50'] != 0) and (all_species[species]['total_within_500'] != 0):
        outfile.write(str(species) + '\t' +
              str(all_species[species]['r2_sum_within_50'] / all_species[species]['total_within_50']) + '\t' + 
              str(all_species[species]['r2_sum_within_500'] / all_species[species]['total_within_500']) + '\t' +
              str(all_species[species]['nucleotide_diversity_sum'] / all_species[species]['total_positions']) + '\n')
outfile.close()
