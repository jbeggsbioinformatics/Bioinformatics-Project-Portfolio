from Bio import SeqIO
import random
import os
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import *
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib
import matplotlib
# Use the Agg backend which is non-interactive and doesn't pop up a window
matplotlib.use('Agg')
import numpy as np
import datetime
from matplotlib.backends.backend_pdf import PdfPages
import pandas

#Function to return command line arguments and print help message
def get_arguments_a3():
    import argparse
    parser = argparse.ArgumentParser(description="Assignment 3")
    parser.add_argument("-protname", type=str, help="Protein name",required=True)
    parser.add_argument("-blastcom", type=str, help="Blast executable",required=True)
    parser.add_argument("-musclecom", type=str, help="muscle executable",required=True)
    parser.add_argument("-fastafile", type=str, help="Source fasta file", required=True)
    parser.add_argument("-databasefasta", type=str, help="Target database file", required=True)
    args = parser.parse_args()
    return [arg for arg in vars(args).values()]

if __name__ == "__main__":
    input_args = get_arguments_a3()
    protein_name = input_args[0]
    blast_exe = input_args[1]
    muscle_exe = input_args[2]
    human_fasta = input_args[3]
    database_fasta = input_args[4]

print("Pipeline running with protein name",protein_name)
print("Running with BLAST as",blast_exe)
print("Running with muscle as",muscle_exe)
print("Source fasta:",human_fasta)
print("Database to search against:",database_fasta)
### Put your code to execute below

#1. Import required packages

import random
import os
from Bio import SeqIO

#2. Define file names and paths

# Generate random number
RND = random.randint(1, 10000)

# Define file names
blast_input = f"temp_{RND}_blast_input.txt"
blast_output = f"temp_{RND}_blast_output.tsv"
msa_input = f"temp_{RND}_msa_input.fasta"
msa_output = f"temp_{RND}_msa_output.fasta"

# BLAST and Muscle executable paths
blast_exe = "/Users/alexbeggs/miniconda3/envs/life733/bin/blastp"  # Update to the path where BLAST+ is installed
muscle_exe = "/Users/alexbeggs/miniconda3/envs/life733/bin/muscle"  # Update to the path where MUSCLE is installed

#3. Find protein for blast search and write to blast input

# Parse human proteome FASTA file to find the protein of interest
fasta_file = "UP000005640_9606.fasta"
# Replace with protein name to use as query for the blast search

found_record = None
for record in SeqIO.parse(fasta_file, "fasta"):
    if protein_name in record.description.split('|')[-1]:
        found_record = record
        break

# Write the found protein to the BLAST input file
if found_record:
    with open(blast_input, "w") as output_file:
        SeqIO.write(found_record, output_file, "fasta")
else:
    raise ValueError(f"Protein {protein_name} not found in {fasta_file}.")

command = blast_exe + " -db UP_mammals -query "+ blast_input +" -out " + blast_output + " -outfmt 6"

print(f"Running BLAST search: {command}")
exit_code = os.system(command)

# Check if BLAST command ran successfully
if exit_code == 0:
    print(f"BLAST search completed successfully, output written to {blast_output}.")
else:
    print(f"Error running BLAST search. Exit code: {exit_code}")

# Manipulating blast output file

#Setting file paths

blast_result_file = blast_output #chose name for output file to be manipulated
database_fasta = '10_uniprot_mammals.fasta'


#Dictionary to store best hit of each species
species_best_hits = {}

#Processing blast output file

with open(blast_result_file, 'r') as f:
    for line in f:
        data_row = line.strip().split('\t')
        query = data_row[0]
        hit = data_row[1]
        e_value = float(data_row[10])
        blast_score = float(data_row [11])

        #extract species name using slicing
        species = hit.split('|')[-1].split('_')[-1]

        #Specifiy only hits with e-value < 1e-05
        if (e_value < 1e-05):
        #use conditional statement to check if species has already been stored,
        #if not or if this hit is more significant, update it

            if species not in species_best_hits or species_best_hits[species]['e_value'] > e_value:
                species_best_hits[species] = {
                    'query': query,
                    'hit': hit,
                    'e_value': e_value,
                    'blast_score': blast_score
                }

#Print best blast hit of each species
#for species, hit_info in species_best_hits.items():
    #print(f"Species: {species}, Best Hit: {hit_info}")

#Load uniprot mammal database fasta and store records in dictionary
id_to_record = {record.id: record for record in SeqIO.parse(database_fasta, "fasta")}

#create list to hold hits you want to take forward for the MSA
records_hit = []

#append sequences from best hits into list
for species, hit_info in species_best_hits.items():
    hit_id = hit_info['hit']
    if hit_id in id_to_record:
        records_hit.append(id_to_record[hit_id])

# Write the selected records to a new FASTA file
with open(msa_input, 'w') as msa_output_handle:
    SeqIO.write(records_hit, msa_output_handle, 'fasta')

print(f"Filtered hits written to {msa_input}")


#Part 2: Multiple sequence alignment

muscle_command = muscle_exe + " -align " + msa_input + " -output "+ msa_output

print(f"Running MUSCLE alignment: {muscle_command}") #debugging step to ensure its running

exit_code = os.system(muscle_command
                      ) #execute muscle command

# Check if the MSA ran successfully
if exit_code == 0:
    print(f"MUSCLE alignment completed successfully. Output written to {msa_output}.")
else:
    print(f"Error running MUSCLE alignment. Exit code: {exit_code}")

#Read alignment using AlignIO

from Bio import AlignIO
alignment = AlignIO.read(msa_output, "fasta")

#Create list of amino acid conservation scores

AA_conservation_scores = []

#Loop through each position in the alignment

for i in range(alignment.get_alignment_length()):
    column = [record.seq[i] for record in alignment] #column is amino acids at this position

    #Create dictionary to hold amino acid counts

    AA_counts = {}

    #Count frequency of each amino acid in the column

    for AA in column:
        if AA in AA_counts:
            AA_counts[AA] += 1
        else:
            AA_counts[AA] = 1

    #Find the most common amino acid and its frequency

    most_common_AA = max(AA_counts, key=AA_counts.get)
    most_common_count = AA_counts[most_common_AA]

    #Calculate conservation ratio and append to conservation score list

    conservation_ratio = most_common_count / len(column)
    AA_conservation_scores.append(conservation_ratio)

# Plot the conservation scores
plt.figure(figsize=(10, 6))
plt.plot(AA_conservation_scores, label="Conservation Score", color="orange")
plt.xlabel("Position in Alignment")
plt.ylabel("Conservation")
plt.title("Amino Acid Conservation Plot of Alignment")
plt.ylim(0, 1.2)  # Conservation ratio ranges from 0 to 1
plt.grid(True)
plt.legend()

# Save the plot to a file
plot_filename = f"{RND}_consensus_plot.png"
plt.savefig(plot_filename)
plt.close()

print(f"Conservation plot saved as {plot_filename}")

os.remove(msa_input)
os.remove(msa_output)

#Phylogenetics Part

#NJ tree construction

calculator = DistanceCalculator("blosum62") #use blosum62 matrix
dm = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

#Remove inner from clades
for clade in tree.find_clades(): #for loop removes inner
    if "Inner" in clade.name:
        clade.name = ""

#Set branch colors based on branch length
for clade in tree.find_clades():
    if clade.branch_length:  # Check if branch length exists
        if clade.branch_length < 0.01:
            clade.color = "brown"
        elif clade.branch_length < 0.02:
            clade.color = "orange"
        elif clade.branch_length < 0.03:
            clade.color = "maroon"
        else:
            clade.color = "red"
    else:
        clade.color = "brown" #make backbone brown


#Function to show just gene name

def extract_gene_name(clade):
    if clade.name:  # Ensure the clade has a name
        return clade.name.split('|')[-1]  # Keep only the part after '|'
    return None

tree_nj = f"temp_{RND}_phylo_nj.png" #define output file name

# Draw the tree
fig = plt.figure(figsize=(8, 8))
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes, label_func=extract_gene_name)  # Use custom label function

#Save figure as png
fig.savefig(tree_nj, format='png', bbox_inches='tight', dpi=300)
plt.close(fig)#close figure to stop it opening automatically

print(f"NJ tree saved to {tree_nj}") #debugging to ensure code runs



#UPGMA tree construction

calculator = DistanceCalculator("blosum62") #use blosum62 matrix
dm = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor()
upgma_tree = constructor.upgma(dm)

#Remove inner from clades
for clade in upgma_tree.find_clades(): #for loop removes inner
    if "Inner" in clade.name:
        clade.name = ""

#Set branch colors based on branch length
for clade in upgma_tree.find_clades():
    if clade.branch_length:  # Check if branch length exists
        if clade.branch_length < 0.01:
            clade.color = "brown"
        elif clade.branch_length < 0.02:
            clade.color = "orange"
        elif clade.branch_length < 0.03:
            clade.color = "maroon"
        else:
            clade.color = "red"
    else:
        clade.color = "brown" #make backbone brown


#Function to show just gene name

def extract_gene_name(clade):
    if clade.name:  # Ensure the clade has a name
        return clade.name.split('|')[-1]  # Keep only the part after '|'
    return None

tree_upgma = f"temp_{RND}_phylo_UPGMA.png"

# Draw the tree
fig = plt.figure(figsize=(8, 8))
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(upgma_tree, axes=axes, label_func=extract_gene_name)  # Use custom label function

#Save figure as png
fig.savefig(tree_upgma, format='png', bbox_inches='tight', dpi=300)
plt.close(fig)#close figure to stop it opening automatically

print(f"UPGMA tree saved to {tree_upgma}") #debugging to ensure code runs

#Extra plots: Barchart of hits per species

# Define a list of colors for the bars (can be customized)
colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'pink', 'brown', 'gray', 'yellow']


species_counts = {}

with open(blast_output, 'r') as blast_file:
    for line in blast_file:
        data_row = line.strip().split('\t')
        species = data_row[1].split('|')[-1].split('_')[-1]  # Extract species
        species_counts[species] = species_counts.get(species, 0) + 1

# Plot species hit counts
species_names = list(species_counts.keys())
hit_counts = list(species_counts.values())

plt.figure(figsize=(10, 6))
plt.bar(species_names, hit_counts, color=colors, edgecolor='black')
plt.xticks(rotation=90)  # Rotate x-axis labels
plt.title('Number of Hits per Species')
plt.xlabel('Species')
plt.ylabel('Hit Count')

species_bar_file = f"temp_{RND}_species_hits.png"
plt.savefig(species_bar_file, format='png', dpi=300, bbox_inches='tight')
plt.close()
print(f"Species hit count plot saved to {species_bar_file}")




os.remove(blast_result_file)
os.remove(blast_input)

#Add all to pdf

pipeline_pdf = f"Protein_analysis_report_{RND}.pdf"

import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Define your generated plot filenames (ensure these are correct)
conservation_plot_file = f"{RND}_consensus_plot.png"
nj_tree_file = f"temp_{RND}_phylo_nj.png"
upgma_tree_file = f"temp_{RND}_phylo_UPGMA.png"
species_bar_file = f"temp_{RND}_species_hits.png"

# Create the PdfPages object to which we will save the pages
with PdfPages(pipeline_pdf) as pdf:
    # Conservation Score Plot
    plt.figure(figsize=(10, 6))
    img = plt.imread(conservation_plot_file)
    plt.imshow(img)
    plt.axis('off')  # Turn off axis for better presentation
    plt.title("Amino Acid Conservation Plot")
    pdf.savefig()  # Saves the figure into a pdf page
    plt.close()

    # NJ Phylogenetic Tree
    plt.figure(figsize=(8, 8))
    img = plt.imread(nj_tree_file)
    plt.imshow(img)
    plt.axis('off')  # Turn off axis
    plt.title("Neighbor Joining Phylogenetic Tree")
    pdf.savefig()
    plt.close()

    # UPGMA Phylogenetic Tree
    plt.figure(figsize=(8, 8))
    img = plt.imread(upgma_tree_file)
    plt.imshow(img)
    plt.axis('off')  # Turn off axis
    plt.title("UPGMA Phylogenetic Tree")
    pdf.savefig()
    plt.close()

    # Species Hit Count Bar Plot
    plt.figure(figsize=(10, 6))
    img = plt.imread(species_bar_file)
    plt.imshow(img)
    plt.axis('off')  # Turn off axis
    plt.title("Number of Blast Hits per Species")
    pdf.savefig()
    plt.close()

    # Optional: Add metadata to the PDF file
    d = pdf.infodict()
    d['Title'] = 'Protein Analysis Results'
    d['Author'] = 'John Alexander Beggs'
    d['Subject'] = 'Analysis of BLAST, MSA, and Phylogenetic Results'
    d['Keywords'] = 'BLAST, MSA, Phylogenetic Tree, Protein Analysis'
    d['CreationDate'] = datetime.datetime(2024, 11, 28)
    d['ModDate'] = datetime.datetime.today()

print(f"All the plots have been successfully integrated into a comprehensive {pipeline_pdf}.")


#Os.remove other files, pipeline produces single pdf
os.remove(plot_filename)
os.remove(tree_nj)
os.remove(tree_upgma)
os.remove(species_bar_file)