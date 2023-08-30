# extra_NGS_related
This repository contains simple codes and functions often used in NGS (next Gen seq) and TGS (Third gen seq) analyses

## 1. Bulk download sequences from NCBI using python

folder called [python_download_sequences](https://github.com/shivanshss/extra_NGS_related/tree/main/python_download_sequences_NCBI)

NCBI Sequence Downloader

This Python script is designed to download nucleotide sequences from the NCBI database using accession numbers. It uses the Biopython library to interact with the NCBI Entrez system and fetch sequences in FASTA format. Accession numbers and corresponding names are provided in a text file, and the script downloads the sequences and saves them as FASTA files.

Usage:

    Ensure you have Python 3.x installed on your system.
    Install the Biopython library if you haven't already. You can do this using the command: pip install biopython.
    Replace "your_email_address" in the script with your actual email address. This is required by NCBI Entrez.
    Prepare a text file with accession numbers and names in tab-separated format. Example: accession_list.txt.
    Set the input_file variable in the script to the path of your input file (e.g., input_file = "/path/to/accession_list.txt").
    Specify the output_directory where you want to save the downloaded sequences.
    Run the script using the command: python script_name.py.

Script Functionality:

    The script utilizes the provided Biopython library and NCBI Entrez API to fetch sequences by accession numbers.
    Accession numbers and corresponding names are read from the input file.
    The script attempts to download sequences for each accession number and saves them as individual FASTA files in the specified output directory.
    If any errors occur during the download process, appropriate messages are displayed.

Note:

    This script assumes that you have access to the NCBI Entrez system and that the provided accession numbers are valid.
    Make sure to replace "your_email_address" with your actual email address to comply with NCBI's usage policy.
    The script requires Biopython to be installed, which can be done using the command: pip install biopython.

Feel free to use this script to conveniently download nucleotide sequences from NCBI using accession numbers for your research purposes.

## 2. Split multi-sample multi-chromosome vcf file by chromosomes

[vcf_to_chr.sh](https://github.com/shivanshss/extra_NGS_related/blob/main/split_vcf/vcf_to_chr.sh) file

Split VCF by Chromosome

This script is designed to split an uncompressed VCF (Variant Call Format) file into separate files based on chromosomes. It employs tabix and parallel processing to efficiently extract entries for each chromosome present in the VCF.

Usage:

    Ensure that the script has execution permissions. If not, you can grant it using the command: chmod +x split_vcf_by_chr.sh.
    Execute the script by providing the uncompressed VCF file as an argument: ./split_vcf_by_chr.sh input.vcf.

Script Functionality:

    The script first checks if the provided file is an uncompressed VCF. If not, an error message is displayed.
    It extracts the list of chromosomes present in the VCF and stores them in a text file.
    The script then sorts the VCF entries based on chromosome and position.
    The sorted entries are compressed and indexed using bgzip and tabix commands, resulting in a compressed and indexed VCF file.
    The script then proceeds to extract separate VCF files for each chromosome using parallel processing.

Note:

    This script assumes that the input VCF file is uncompressed.
    The extracted VCF files for each chromosome are named as ${input.vcf}_${chromosome}.vcf.

Feel free to use this script to split your VCF files into chromosome-specific files for further analysis or downstream processing.

## 3. random_seq_generator.py

This script generates two random DNA sequences of specified lengths ('s_length' and 't_length') and displays them, which can be useful for simulating data in bioinformatics analysis. Adjust the lengths as needed for your analysis.


## 4. Counting Disease Carriers
[AFRQ.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/AFRQ/AFRQ.py)
The given script calculates the probabilities that a randomly selected individual carries at least one copy of the recessive allele for each Mendelian factor in a diploid population. It takes an array A as input, where A[k] represents the proportion of homozygous recessive individuals for the k-th factor. The function calculate_probabilities computes the probabilities using the complement rule and stores them in array B. The sample input A is provided as an example, and the script prints the resulting array B. Feel free to replace the sample input A with your own array of proportions and run the script to calculate and print the corresponding probabilities.


## 5. Assessing Assembly Quality with N50 and N75
[ASMQ.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/ASMQ/ASMQ.py)

This script calculates the N50 and N75 statistics for a collection of DNA strings. It takes a set of DNA strings as input and calculates the lengths of contigs that contribute to 50% and 75% of the total combined length. Adjust the input strings accordingly, and then run the script to obtain the N50 and N75 values for the provided collection.


## 6. Global Multiple Alignment
[CLUS.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/CLUS/CLUS.py) 
This script calculates the most different nucleotide sequence in a set of sequences provided in FASTA format. It uses the Biopython library to perform sequence alignment and identity calculation. Run the script by providing an "input.fasta" file containing the nucleotide sequences. The script will identify the sequence with the highest average identity difference from the others and output its ID.


## 7. Error Correction in Reads
[CORR.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/CORR/CORR.py)
This script identifies sequencing errors in a collection of DNA reads and provides error corrections in the form of "[old read]->[new read]", where each correction represents a single symbol substitution. The script utilizes a dictionary to count the occurrences of each read and then identifies reads that are potentially erroneous based on their frequency and Hamming distance from other reads. To use the script, provide the DNA reads in FASTA format as input. The script will output a list of error corrections. Each correction represents a potentially erroneous read and its corrected version, taking into account Hamming distance and possible reverse complements. Make sure you have the input data in a file named "input.fasta" in the same directory as the script before running it.


## 8. Counting Optimal Alignments 
[CTEA.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/CTEA/CTEA.py)
This script calculates the total number of optimal alignments between two protein strings, 's' and 't', provided in FASTA format. It employs dynamic programming to compute the counts of alignments while considering the edit alignment score. Run the script with an "input.fasta" file containing the protein strings. The script will output the number of optimal alignments modulo 134,217,727 (2^27 - 1).


## 9. Constructing a De Bruijn Graph
[DBRU.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/DBRU/DBRU.py)
This script constructs the adjacency list of a de Bruijn graph based on a collection of DNA strings of equal length (up to 50 base pairs) and a given parameter 'k' (k-mer length). It builds the graph by considering overlapping k-mers and connecting them based on their subsequences. Run the script, providing 'k' as input followed by the DNA strings. The script will output the adjacency list of the de Bruijn graph representing the set of k-mers.


## 10. Edit Distance 
[EDIT.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/EDIT/EDIT.py)
This script calculates the edit distance between two protein strings, 's' and 't', provided in FASTA format. It uses dynamic programming to compute the minimum number of operations (insertions, deletions, or substitutions) required to transform one string into the other. To run the script, provide an "input.fasta" file containing the protein strings. The script will output the edit distance between the two strings.


## 11. Edit Distance Alignment 
[EDTA.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/EDTA/EDTA.py)
This script calculates the edit distance between two protein strings, 's' and 't', in FASTA format. It then constructs the optimal alignment of the two strings by performing dynamic programming and backtracking. The output includes the edit distance, as well as the aligned strings 's' and 't'. To use the script, provide an "input.fasta" file containing the protein strings. The script will output the edit distance followed by the optimal alignment of the two strings.


## 12. The Founder Effect and Genetic Drift
[FOUN.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/FOUN/FOUN.py)
The given script calculates an m×k matrix B where each element Bi,j represents the common logarithm of the probability that after i generations, no copies of the recessive allele for the j-th factor will remain in the population. It uses the Wright-Fisher model, and the function calculate_probability computes the probabilities based on the binomial distribution and logarithms. The sample input values N, m, and A are provided as an example. You can replace them with your own values and run the script to calculate and print the resulting matrix B.


## 13. Global Alignment with Scoring Matrix and Affine Gap Penalty
[GAFF.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/GAFF/GAFF.py)
This script calculates the maximum alignment score between two protein strings, 's' and 't', in FASTA format. It uses the BLOSUM62 scoring matrix and performs global sequence alignment with an affine gap penalty. To use the script, provide an "input.fasta" file containing the protein strings. The script will output the maximum alignment score between the two strings.


## 14. Genome Assembly Using Reads 
[GASM.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/GASM/GASM.py)
This script constructs a cyclic superstring of minimal length containing a collection of error-free reads of equal length provided as input. It does so by utilizing the de Bruijn graph formed by the reads and their reverse complements. The script iteratively identifies and appends the next character to the cyclic superstring until both directed cycles in the graph are formed. To use the script, provide the error-free reads as input. The script will output the constructed cyclic superstring.


## 15. Global Alignment with Constant Gap Penalty 
[GCON.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/GCON/GCON.py)
This script calculates the maximum alignment score between two protein strings, 's' and 't', in FASTA format. It employs the BLOSUM62 scoring matrix and performs global sequence alignment with a constant gap penalty of -5. To use the script, provide an "input.fasta" file containing the protein strings. The script will output the maximum alignment score between the two strings.


## 16. Global Alignment with Scoring Matrix 
[GLOB.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/GLOB/GLOB.py)
This script calculates the maximum alignment score between two protein strings, 's' and 't', in FASTA format. It performs this by counting the total number of optimal alignments using dynamic programming with linear gap penalties. The BLOSUM62 scoring matrix is used for alignment, and the gap penalty is set to -5. To use the script, provide an "input.fasta" file containing the protein strings. The script will output the maximum alignment score between the two strings.


## 17. Genome Assembly with Perfect Coverage and Repeats 
[GREP.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/GREP/GREP.py)
This script constructs circular strings assembled from complete cycles in the de Bruijn graph Bk of error-free (k+1)-mers provided as input. It does so by iteratively connecting k-mers based on their suffixes and prefixes. To use the script, provide a list of (k+1)-mers as input. The script will output the circular strings, each beginning with the first (k+1)-mer given in the input.


## 18. Counting Point Mutations
[HAMM.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/HAMM/HAMM.py)
This script calculates and prints the Hamming distance between two DNA strings, 's' and 't', and their complementary version 't' and 's', respectively. The Hamming distance between two strings is the count of positions at which the symbols differ. It is computed for both directions. To use the script, provide the sequences 's' and 't' and run the script. The script will output the Hamming distances for both directions.

## 19. Longest Increasing Subsequence 
[LGIS.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/LGIS/LGIS.py)
This script finds the longest increasing subsequence and the longest decreasing subsequence of a given permutation. The longest increasing subsequence (LIS) is calculated using dynamic programming, and the longest decreasing subsequence (LDS) is derived by reversing the permutation and finding the LIS of the reversed sequence. To use the script, provide a positive integer n followed by a space-separated permutation of length n. The script will then output the longest increasing subsequence and the longest decreasing subsequence of the given permutation.


## 20. Genome Assembly as Shortest Superstring 
[LONG.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/LONG/LONG.py)
This script aims to reconstruct a linear chromosome from a set of DNA strings in FASTA format by iteratively merging overlapping pairs until a shortest superstring containing all input strings is obtained. Provide input in "input.fasta" and run the script to get the reconstructed superstring.
Make sure the DNA strings do not exceed 1 kbp and satisfy the given conditions.
The code uses a greedy approach to iteratively merge overlapping strings until a single superstring is obtained.


## 21. Multiple Alignment 
[MULT.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/MULT/MULT.py)
The script calculates the highest alignment score for a set of four DNA strings (up to 10 bp each) in FASTA format. It uses Biopython's pairwise2 module, assigning 0 for matches and -1 for mismatches, and prints the score. Prepare "input.fasta" with sequences and run the script to find the maximum alignment score.


## 22. Pairwise Global Alignment
[NEED.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/NEED/NEED.py)
The script calculates the maximum global alignment score between DNA sequences associated with two given GenBank IDs. It utilizes the EMBOSS Needle tool for pairwise sequence alignment, considering a gap open penalty of 10 and a gap extension penalty of 1. Replace the example sequences with the actual ones for accurate results. Run the script after providing the GenBank IDs to find and print the maximum alignment score.


## 23. Overlap Alignment 
[OAP.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/OAP/OAP.py)
This script computes the score of an optimal overlap alignment between two DNA strings (s and t) in FASTA format. Using Biopython's pairwise2 module, it calculates alignment scores with a custom scoring scheme and prints the score along with aligned sequences. Prepare "input.fasta" with sequences and run the script to find the optimal overlap alignment score and the corresponding aligned sequences.
Ensure that both DNA strings have lengths of at most 10 kbp.


## 24. Genome Assembly with Perfect Coverage
[PCOV.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/PCOV/PCOV.py)
This script constructs a cyclic superstring of minimal length containing error-free DNA k-mers (k≤50) from a circular chromosome strand. It forms the superstring using a de Bruijn graph and a dictionary to store k-1 mers and their suffixes. The script iteratively appends overlapping k-mers to create the cyclic superstring, which corresponds to a candidate cyclic chromosome. Run the script, providing the k-mers as input, to obtain the cyclic superstring. This script assumes k-mers are provided as input, so make sure to provide them accordingly.


## 25. Creating a Distance Matrix 
[PDST2.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/PDST2/PDST2.py)
This script calculates the p-distance matrix (D) corresponding to the p-distance (dp) between a collection of n DNA strings (n≤10) given in FASTA format. It uses pairwise sequence alignment with no penalty for gaps and counts the proportion of differing characters to compute the p-distance. The resulting matrix represents the pairwise p-distances between input sequences. Run the script, providing "input.fasta" with sequences, to obtain the distance matrix with an error margin of 0.001.


## 26. 
[PERM.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/PERM/PERM.py)
This script calculates the total number of permutations of length n (where n≤7) and generates a list of all such permutations. It utilizes the itertools library to generate permutations and then prints the count and the list of permutations. Run the script, providing n as input, to obtain the desired output. Make sure to provide a positive integer value for n when running the script.



## 27. Enumerating Gene Orders
[PPER.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/PPER/PPER.py)
This script calculates the total number of partial permutations P(n,k) where n and k are positive integers satisfying 100≥n>0 and 10≥k>0. It computes the partial permutation using the formula and then takes the modulo 1,000,000 of the result. Run the script, providing n and k as input, to get the total number of partial permutations modulo 1,000,000. Input positive integers for n and k when running the script.



## 28. Reversal Distance
[REAR.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/REAR/REAR.py)
This script calculates the reversal distance between pairs of permutations. For each pair of permutations provided (up to 5 pairs), it calculates the reversal distance using pancake sorting. The flip function performs the flip operation on an array, and the pancake_sort function applies pancake sorting to sort an array. The reversal_distance function computes the total reversal distance between two permutations. Run the script, providing input pairs of permutations, to obtain the reversal distances for each pair. The input should contain pairs of permutations, each on a separate line.



## 29. Enumerating Oriented Gene Orderings
[SIGN.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/SIGN/SIGN.py)
This script calculates the total number of signed permutations of length n (where n≤6) and generates a list of all such permutations. It utilizes the itertools library to generate permutations and signs for each number. The script then prints the count and the list of signed permutations. Run the script, providing n as input, to obtain the desired output. Make sure to provide a positive integer value for n when running the script.



## 30. Semiglobal Alignment
[SMGB.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/SMGB/SMGB.py)
This script calculates the maximum semiglobal alignment score between two DNA strings (s and t at most 10kbp) in FASTA format. It uses the Biopython library's pairwise2 module to perform semiglobal alignment with a scoring scheme where matches count +1, substitutions count -1, and there's a linear gap penalty of 1. The script then prints the maximum alignment score and the corresponding aligned sequences.



## 31. Sorting by Reversals 
[SORT.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/SORT/SORT.py)
This script calculates the reversal distance drev(π, γ) between two given permutations π and γ of length 10. It finds a collection of reversals that sort permutation π into γ using the pancake sorting algorithm. The flip function is used to reverse segments of an array, and the pancake_sort function performs pancake sorting to sort the permutation. The reversal_distance function calculates the reversal distance between two permutations. Run the script, providing the two permutations as input, to obtain the reversal distance and the collection of reversals. Input should consist of two lines, each containing a space-separated permutation of length 10.



 
## 32. Suboptimal Local Alignment 
[SUBO.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/SUBO/SUBO.py)
This script calculates the total number of occurrences of a given inexact repeat (with up to 3 changes/indels) within two DNA strings s and t provided in FASTA format. It uses the find_occurrences function to identify substrings that closely match the repeat within a certain length range (32-40 bp). The script then prints the count of occurrences for both s and t.



## 33. Transitions and Transversions 
[TRAN2.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/TRAN2/TRAN2.py)
This script calculates the transition/transversion ratio (R) between two given DNA strings s1 and s2 of unequal length (up to 1 kbp). The transition_transversion_ratio function computes the ratio by counting the number of transitions (A↔G or C↔T) and transversions (all other substitutions) between corresponding bases in the two sequences. Run the script, providing the sequences in FASTA format within "input.fasta," to calculate and print the transition/transversion ratio.



## 34. The Wright-Fisher Model of Genetic Drift
[WFMD.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/WFMD/WFMD.py)
The provided script calculates the probability that, in a population of N diploid individuals initially possessing m copies of a dominant allele, after g generations, at least k copies of a recessive allele will be observed. It uses the Wright-Fisher model, and the function calculate_probability computes the probability based on the binomial distribution. The sample input values N, m, g, and k are provided as an example. You can replace them with your own values and run the script to calculate and print the corresponding probability.


## 35.Counting DNA Nucleotides
[DNA.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/DNA/DNA.py)
This script reads a DNA string from a FASTA file and calculates the counts of 'A', 'C', 'G', and 'T' symbols in the string. The count_symbols function calculates the counts, and the process_fasta_file function processes the input FASTA file, extracts sequences, and prints the symbol counts for each sequence.

## 36. Transcribing DNA into RNA
[RNA.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/RNA/RNA.py)
This script reads DNA strings from a multi-FASTA file and transcribes them into RNA sequences. The transcribe_dna_to_rna function replaces 'T' with 'U' to perform the transcription. The process_multi_fasta_file function processes the multi-FASTA file, extracts sequences, and prints the corresponding transcribed RNA sequences.

## 37. reverse complement of a DNA string
[REVC.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/REVC/REVC.py)
This script reads DNA sequences from a multi-FASTA file, calculates the reverse complement for each sequence, and writes the reverse complements to an output FASTA file. The reverse_complement function calculates the reverse complement of a DNA string using the given complement dictionary.


## 38. Rabbits and Recurrence Relations
[FIB.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/FIB/FIB.py)
The script uses the Fibonacci-like recurrence relation to calculate the total number of rabbit pairs after n months, considering the provided reproduction rules. Here's how the script works:

    The function rabbit_pairs(n, k) takes two parameters: n (the number of months) and k (the number of offspring pairs produced by each reproduction-age pair).

    The function starts by handling the base cases. If n is 1 or 2, there will be only one pair of rabbits in the population for both cases.

    For months beyond the base cases (i.e., when n is greater than 2), the script uses an iterative approach to calculate the total number of rabbit pairs.

    The variables prev_prev, prev, and current are used to keep track of the number of pairs for the current month, the previous month, and two months ago, respectively.

    The loop runs from the 3rd month (i starts from 3) up to the desired n month.

    In each iteration, the number of pairs for the current month is calculated using the recurrence relation: current = prev + k * prev_prev.

    The values of prev_prev and prev are updated to move to the next iteration.

    Once the loop completes, the variable current holds the total number of rabbit pairs after n months.

    The script returns the final value of current.

    Example usage is provided at the end of the script. You input values for n and k using input(), and the script calculates and prints the result.


## 39. Computing GC Content
[GC.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/GC/GC.py)
The script calculates the GC-content of DNA sequences in FASTA format, identifies the sequence with the highest GC-content, and outputs its ID along with the calculated percentage.


## 40. Mendel's First Law
[IPRB.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/IPRB/IPRB.py)
The script calculates the probability of producing an individual with a dominant allele phenotype from a population containing homozygous dominant (k), heterozygous (m), and homozygous recessive (n) organisms. It reads input values from a file, computes the probability for each set of values, and outputs the results in tab-separated format.

## 41. Translating RNA into Protein
[PROT.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/PROT/PROT.py)
The script converts RNA sequences into protein sequences based on the codon table. It reads sequences from a multi-FASTA file, translates each RNA sequence to a protein sequence, and prints the results with labels.

## 42. Finding a Motif in DNA
[SUBS.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/SUBS/SUBS.py)
The script takes two input files, each containing strings separated by newlines. It then finds the locations of the second string (t) within each of the first strings (s) and prints the locations in space-separated format.


## 43. Consensus and Profile
[CONS.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/CONS/CONS.py)
The script aligns a collection of DNA strings using pairwise alignment, then builds a profile matrix from the aligned sequences, and generates a consensus string based on the profile. The consensus and profile are then printed.


## 44. Mortal Fibonacci Rabbits
[FIBD.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/FIBD/FIBD.py)
This script calculates the total number of rabbit pairs in a population after n months, considering a dynamic programming approach. It reads input data from a CSV file containing pairs of values n and m, representing the number of months and the lifespan of a rabbit pair, respectively. It then calculates the total number of pairs using the provided function rabbit_pairs_after_n_months, and writes the results to an output text file.


## 45. Overlap Graphs
[GRPH.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/GRPH/GRPH.py)
This script generates an overlap graph from a collection of DNA strings. It reads input data from a FASTA file containing DNA sequences, and then constructs an overlap graph by comparing the last k characters of one string with the first k characters of another string. If an overlap is found, an edge is added to the adjacency list. The output is a list of edges representing the adjacency relationships in the overlap graph.

## 46. Calculating Expected Offspring
[IEV.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/IEV/IEV.py)
This script calculates the expected number of offspring displaying the dominant phenotype given the number of individuals for each genotype pairing in a population. It reads input data from a text file containing the counts of homozygous dominant, heterozygous, and homozygous recessive genotypes, as well as the probabilities of producing offspring with dominant phenotypes for each genotype pairing. The output is the expected number of offspring displaying the dominant phenotype.


## 47. Finding a Shared Motif
[LCSM.py](https://github.com/shivanshss/extra_NGS_related/blob/main/learning_python_Rosalind/LCSM/LCSM.py)
This script finds the longest common substring that is present in all the given DNA strings. It reads input data from a FASTA file containing DNA sequences and extracts the sequences to find the longest common substring. The output is the longest common substring that is shared among all the input DNA sequences.



