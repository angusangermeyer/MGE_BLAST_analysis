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