#!/usr/bin/env python3

"""
-. .-.   .-. .-.   .-. .-.   .  
||\|||\ /|||\|||\ /|||\|||\ /|
|/ \|||\|||/ \|||\|||/ \|||\||
~   `-~ `-`   `-~ `-`   `-~ `-

This program is designed to look for known elements in assembled bacterial genomes. \
Elements can be any reasonable sized nucleotide sequence for a gene, MGE, ORF, etc. \

It will process three types of element lists: 
	1) A folder of different elements in single fasta files.
	2) A folder of different GROUPS of mutually exclusive elements in multi-fasta files (i.e. a file with all PLEs).
	3) A folder of single elements that need to be hit EXACTLY. Useful for stuff like serotype.

You must provide a path to the assemblies, the single, grouped and exact elements, and \
a location for intermediate blast output.


The program uses BLASTn to find likely hits and currently just uses default settings. \
You can adjust the minimum length of contig that will be considered and minumum coverage \
if available in fasta headers (will be ignored of not there). You can also adjust the \
length percent cutoff for single alignment calls. Finally, you can adjust \
the minumum hit length required for a multi-fasta group call. There is a low cutoff and \
a high cutoff. Below the low cutoff it will not report a hit, above the high it will report \
a known match. Between the two and it will make a guess, but inform you that if may be new. \
This is good for detecting new PLEs for example.


Oh, also the assembly files and element files are hashed and that info is stored in the cached blast results. \
This allows the program to easily look up previous results and saves a LOT of computational time. \
If you cange anything about an assembly or element file, including the name, it will generate a new hash \
and run the blast again.

-. .-.   .-. .-.   .-. .-.   .  
||\|||\ /|||\|||\ /|||\|||\ /|
|/ \|||\|||/ \|||\|||/ \|||\||
~   `-~ `-`   `-~ `-`   `-~ `-
"""


#Imported modules
import os
import time
import hashlib
import argparse

from Bio import SeqIO
from pathlib import Path
from subprocess import call
from argparse import RawTextHelpFormatter
from statistics import stdev, mean, median

#Authorship information
__author__ = "Angus Angermeyer"
__copyright__ = "Copyright 2023, Angus Angermeyer"
__license__ = "GPL"
__version__ = "1.0.0"
__email__ = "angus.angermeyer@gmail.com"


#Parser stuff
parser = argparse.ArgumentParser(
	usage='Batch Assembly BLAST Analysis',
	description = __doc__,
	formatter_class=RawTextHelpFormatter)

parser.add_argument('-A', '--assemblies', type=str, required=True, help='Location of the folder with assembly files (must be fasta).')
parser.add_argument('-S', '--singles', type=str, required=True, help='Location of the folder with individual elements.')
parser.add_argument('-G', '--groups', type=str, required=True, help='Location of the folder with grouped elements.')
parser.add_argument('-E', '--exact', type=str, required=True, help='Location of the folder with exact elements.')
parser.add_argument('-R', '--reports', type=str, required=True, help='Location to cache output blast reports.')

parser.add_argument('-L', '--len_cutoff', type=int, default=1000, help='Required contig length. (default 1000 bp)')
parser.add_argument('-C', '--cov_cutoff', type=int, default=5, help='Required contig coverage depth. (default 5x)')
parser.add_argument('-s', '--single_cutoff', type=float, default=0.60, help='Required percent of best hit to consider real. (default 0.6)')
parser.add_argument('-H', '--high_cutoff', type=float, default=0.80, help='Required percent of best hit to consider real. (default 0.8)')
parser.add_argument('-W', '--low_cutoff', type=float, default=0.50, help='Required percent of best hit to consider real. (default 0.5)')


args = parser.parse_args()


len_cutoff = args.len_cutoff
cov_cutoff = args.cov_cutoff
single_cutoff = args.single_cutoff
guess_hi_cutoff = args.high_cutoff
guess_lo_cutoff = args.low_cutoff

genome_directory = args.assemblies
additive_single_dir = args.singles
additive_guess_dir = args.groups
exact_single_dir = args.exact
existing_report_dir = args.reports

out_time = time.strftime('%d-%m-%y_%H-%M-%S')
out_file = f"genome_blast_results_{out_time}.txt"

try:
	os.makedirs(f"{existing_report_dir}/cache")
except FileExistsError:
	# directory already exists
	pass



##########################################################################################	
##########################################################################################
def ref_file_hash(fasta_file):
	with open(fasta_file, 'rb') as input_file:
		hash = hashlib.md5()
		hash.update(input_file.read())
		hash = hash.hexdigest()
	return(hash)


def recursive_merge(intervals):
	"""https://stackoverflow.com/questions/43600878/merging-overlapping-intervals"""
	
	intervals.sort(key=lambda interval: interval[0])
	merged = [intervals[0]]
	for current in intervals:
		previous = merged[-1]
		if current[0] <= previous[1]:
			previous[1] = max(previous[1], current[1])
		else:
			merged.append(current)

	summed_ranges = sum([j-i for i,j in merged])
	return summed_ranges



def additive_blast_single(genome_info, ref_info_entry):
	genome_file = genome_info[0]
	genome_stem = genome_info[1]
	genome_hash = genome_info[2]	
	
	ref_file = ref_info_entry[0]
	ref_stem = ref_info_entry[1]
	ref_hash = ref_info_entry[2]
	ref__len = ref_info_entry[3]
	
	report_file = f"{genome_stem}_{genome_hash}_{ref_stem}_{ref_hash}.txt"

	if Path(f"{existing_report_dir}/cache/{report_file}").is_file():
		pass

	else:

		blast_call = ("blastn "
		f"-subject {genome_file} "
		f"-query {ref_file} "
		f"-outfmt '6 qaccver saccver pident qlen length slen mismatch gapopen qstart qend sstart send sstrand evalue' "
		f"> {existing_report_dir}/cache/{report_file}"
		)
	
	
		call(blast_call, shell = True)



	with open(f"{existing_report_dir}/cache/{report_file}", "r") as blast_result:
		result_length = 0
		hit_ranges = []
		
		for line in blast_result:
			if line[0] != "#":
				line = line.split()
				hit_len = int(line[4])
				contig_len = int(line[5])
				
				query_start = int(line[8])
				query_end = int(line[9])
				
				try:				
					coverage = float(line[1][line[1].index("cov_")+4:])
				except:
					coverage = cov_cutoff #Always accept contig if no cov info in assembly
					
				if contig_len >= len_cutoff and coverage >= cov_cutoff:
					hit_ranges.append(sorted([query_start,query_end]))

	if len(hit_ranges) > 0:
		result_length = recursive_merge(hit_ranges)
	
	else:
		result_length = 0
		
	hit_pct = round(float(result_length/ref__len), 4)
	
	if hit_pct >= single_cutoff:
		return hit_pct
	else:
		return "-"
	


def exact_blast_single(genome_info, ref_info_entry):
	genome_file = genome_info[0]
	genome_stem = genome_info[1]
	genome_hash = genome_info[2]	
	
	ref_file = ref_info_entry[0]
	ref_stem = ref_info_entry[1]
	ref_hash = ref_info_entry[2]
	ref__len = ref_info_entry[3]
	
	report_file = f"{genome_stem}_{genome_hash}_{ref_stem}_{ref_hash}.txt"

	if Path(f"{existing_report_dir}/cache/{report_file}").is_file():
		pass

	else:

		blast_call = ("blastn "
		f"-subject {genome_file} "
		f"-query {ref_file} "
		f"-outfmt '6 qaccver saccver pident qlen length slen mismatch gapopen qstart qend sstart send sstrand evalue' "
		f"-max_hsps 1 " ###FIX THIS FOR FINAL VERSION (cuts some MGEs in half e.g. O1 sometimes)###
		f"> {existing_report_dir}/cache/{report_file}"
		)
	
		call(blast_call, shell = True)


	with open(f"{existing_report_dir}/cache/{report_file}", "r") as blast_result:
		result_length = 0
		
		result_list = []
		hit_list = []
		for line in blast_result:
			if line[0] != "#":
				line = line.split()
				query = line[0]
				qlen_len = int(line[3])
				hit_len = int(line[4])
				contig_len = int(line[5])
				mismatch = int(line[6])
				
				hit_list.append([mismatch,query])
				
				if hit_len == qlen_len and mismatch == 0:
					result_list.append(query)
				
				
	if len(result_list) == 0:
		return("-")

	else:
		return(",".join(result_list))





def additive_blast_guess(genome_info, ref_info_entry):
	genome_file = genome_info[0]
	genome_stem = genome_info[1]
	genome_hash = genome_info[2]	
	
	ref_file = ref_info_entry[0]
	ref_stem = ref_info_entry[1]
	ref_hash = ref_info_entry[2]

	guess_entry_dict = ref_info_entry[3]
	
	report_file = f"{genome_stem}_{genome_hash}_{ref_stem}_{ref_hash}.txt"
	
	blast_list = [[key,0,0] for key in guess_entry_dict]


	if Path(f"{existing_report_dir}/cache/{report_file}").is_file():
		pass

	else:
		blast_call = ("blastn "
		f"-subject {genome_file} "
		f"-query {ref_file} "
		f"-outfmt '6 qaccver saccver pident qlen length slen mismatch gapopen qstart qend sstart send sstrand evalue' "
		f"> {existing_report_dir}/cache/{report_file}"
		)
	
		call(blast_call, shell = True)


	with open(f"{existing_report_dir}/cache/{report_file}", "r") as blast_result:
		pct_ID_list = []
		
		for line in blast_result:
			line = line.split()
			sub_key = line[0]
			qry_key = line[1]
			prct_ID = float(line[2])
			hit_len = int(line[4])
			qry_len = int(line[5])
			mismtch = int(line[6])
			gap_num = int(line[7])
			
			try:
				coverage = float(line[1][line[1].index("cov_")+4:])
			except:
				coverage = cov_cutoff #Always accept contig if no cov info in assembly

			if qry_len >= len_cutoff and coverage >= cov_cutoff and hit_len>150:
				pct_ID_list.append(prct_ID)

				for entry in blast_list:
					if sub_key == entry[0]:
						entry[1] += hit_len
						entry[2] += 1
	
	
	
	guess_list = []
	
	for entry in blast_list:
		guess_list.append(entry[1]/guess_entry_dict[entry[0]])
		entry[1] = entry[1]/guess_entry_dict[entry[0]]


	#Need to find a way to use percent ID to refine best guess...
	if max(guess_list) > guess_lo_cutoff:
		len_ratio = round(max(guess_list), 2)
		
		if max(guess_list) > guess_hi_cutoff:
			best_guess = blast_list[max(range(len(guess_list)), key=guess_list.__getitem__)][0]
		else:
			best_guess = f"Unk_{blast_list[max(range(len(guess_list)), key=guess_list.__getitem__)][0]}-like-{len_ratio}"
	
	
		if len(pct_ID_list) > 1:
			avg_pct_ID = round(mean(pct_ID_list), 2)
		elif len(pct_ID_list) == 1:
			avg_pct_ID = pct_ID_list[0]
		
		
	else:
		best_guess = "-"
		avg_pct_ID = "-"
		len_ratio = "-"


	return([best_guess, avg_pct_ID, len_ratio, blast_list])




##########################################################################################	
##########################################################################################
#Initialize output document and read in reference file details
with open(out_file, "w") as output:

	output.write(f"Genome\t")
	output.write(f"Assembly_length\t")
	output.write(f"Contigs(>{len_cutoff})\t")
	output.write(f"Coverage_median\t")
	output.write(f"Coverage_mean\t")
	output.write(f"Coverage_stdev\t")
	output.write(f"N50\t")
	output.write(f"L50\t")

	
	
	######################################################################################
	#Get details for individual presence/absence references
	additive_single_list = [] 
	for ref_file in Path(additive_single_dir).glob("*.fasta"):
		ref_stem = ref_file.stem
		ref_hash = ref_file_hash(ref_file)
		ref__len = len(SeqIO.read(ref_file, "fasta").seq)
	
		ref_info = [ref_file, ref_stem, ref_hash, ref__len]
		additive_single_list.append(ref_info)

		output.write(f"{ref_stem}\t")

	######################################################################################
	#Get details for individual presence/absence references
	exact_single_list = [] 
	for ref_file in Path(exact_single_dir).glob("*.fasta"):
		ref_stem = ref_file.stem
		ref_hash = ref_file_hash(ref_file)
# 		ref__len = len(SeqIO.read(ref_file, "fasta").seq)
		ref__len = 0
		
		ref_info = [ref_file, ref_stem, ref_hash, ref__len]
		exact_single_list.append(ref_info)

		output.write(f"{ref_stem}\t")

	######################################################################################
	#Get details for best guess from group references
	additive_guess_list = []
	for ref_file in Path(additive_guess_dir).glob("*.fasta"):
		ref_stem = ref_file.stem
		ref_hash = ref_file_hash(ref_file)
		
		##################################################################################
		#Create a special report for each guess database file
		with open(f"{existing_report_dir}/{ref_stem}_report_{out_time}.txt", "w") as additive_guess_report:
			additive_guess_report.write(f"Genome\t")
			guess_entry_dict = {}
			for record in SeqIO.parse(ref_file, "fasta"):
				guess_entry_dict[record.description] = len(record.seq)
				additive_guess_report.write(f"{record.description}_len_ratio\t")
				additive_guess_report.write(f"{record.description}_#contigs\t")

# 			additive_guess_report.write("hit_index\t")
			additive_guess_report.write("best_guess\t")
			additive_guess_report.write("\n")


		ref_info = [ref_file, ref_stem, ref_hash, guess_entry_dict]
		additive_guess_list.append(ref_info)

		output.write(f"{ref_stem}\t")
		output.write(f"Avg_%ID\t")
		output.write(f"Len_ratio\t")

	output.write(f"\n")




##########################################################################################	
##########################################################################################
#Iterate over genomes
for genome_file in Path(genome_directory).glob("*.fasta"):
	genome_stem = genome_file.stem
	genome_hash = ref_file_hash(genome_file)
	print(genome_stem)
	
	######################################################################################
	#Compute coverage, contig number, total length, and stats
	cov_list = []
	len_list = []
	contig_count = 0
	
	for contig_record in SeqIO.parse(genome_file, "fasta"):
		contig_count += 1
		len_list.append(len(contig_record.seq))
		
		#Some contig headers have stuff after cov info, this should strip that out...
		try:
			raw_string = contig_record.description[contig_record.description.index("cov_")+4:]
			cov_string = ""
			for character in raw_string:
				if character in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '.']:
					cov_string += character
				else:
					break
			coverage = float(cov_string)

		except:
			coverage = 0
		
		cov_list.append(coverage)
	
	######################################################################################
	#If coverage data is not available, just write no data to field
	cov_mean = round(mean(cov_list), 2)
	cov_median = round(median(cov_list), 2)
	
	
	if cov_mean == 0:
		cov_mean = "ND"
		cov_median = "ND"
		cov_stdev = "ND"
	else:
		
		cov_stdev = round(stdev(cov_list), 2)
	
	######################################################################################
	#Compute N50 and L50
	len_list.sort(reverse=True)
	iteration = 0
	for item in len_list:
		if sum(len_list[:iteration+1]) > sum(len_list)/2:				
			N50 = item
			L50 = iteration + 1
			break
		else:
			iteration += 1
	
	
	######################################################################################
	#Create list of all genome data needed for blast analysis subroutines
	genome_info = [genome_file, genome_stem, genome_hash]

	
	######################################################################################
	#Write out genome info to results file
	with open(out_file, "a") as output:
		output.write(f"{genome_stem}\t")
		output.write(f"{sum(len_list)}\t")
		output.write(f"{contig_count}\t")
		output.write(f"{cov_median}\t")
		output.write(f"{cov_mean}\t")
		output.write(f"{cov_stdev}\t")
		output.write(f"{N50}\t")
		output.write(f"{L50}\t")

		
		for ref_info_entry in additive_single_list:
			output.write(str(additive_blast_single(genome_info, ref_info_entry))+"\t")

		for ref_info_entry in exact_single_list:
			output.write(str(exact_blast_single(genome_info, ref_info_entry))+"\t")

		for ref_info_entry in additive_guess_list:
			results = additive_blast_guess(genome_info, ref_info_entry)
			
			output.write(str(results[0])+"\t")
			output.write(str(results[1])+"\t")
			output.write(str(results[2])+"\t")
			
			with open(f"{existing_report_dir}/{ref_info_entry[1]}_report_{out_time}.txt", "a") as additive_guess_report:
				additive_guess_report.write(genome_stem+"\t")
				for item in results[3]:
					additive_guess_report.write(str(item[1])+"\t")
					additive_guess_report.write(str(item[2])+"\t")
			
				additive_guess_report.write(str(results[0])+"\n")

		output.write(f"\n")

