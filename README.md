## genes_detective - a tool for working with nucleic acid sequences 🔎 ##
**genes_detective** allows you to make various transformations with DNA and RNA, as well
as filter fastq sequences according to specified parameters

![Alt текст](https://a.d-cd.net/DqAAAgEQDeA-960.jpg)



## Content

* [Installation](#installation)
* [Main Features](#main-features)
* [Examples](#examples)
* [Acknowledgements](#acknowledgements)



## Installation

Run the following commands:  
`git clone git@github.com:eveelinss/genes_detective.git`  
`cd genes_detective`  
You have gained access to the program `genes_detective` and the folder with modules  
Use the code as a module  
Profit!



## Main Features

`run_dna_rna_tools` returns the transformed dna and rna sequences depending 
on the operation to be performed:  
    - reverse  
    - complement  
    - reverse_complement  
    - transcribe  
    - search_start_codon_in_rna  
    - search_first_codon_in_rna  

`filter_fastq` filtering fastq sequences in file by GC composition, length, and 
quality threshold and writes the filtering result to a file in the `filtered` folder.

`convert_multiline_fasta_to_oneline` receives a fasta file as input, in which 
the sequence of nucleic acid or protein is broken down and broken down sequentially 
on different lines. The function overwrites the data in the correct format in a new file 

`parse_blast_output` receives the input path to the file containing the 
result of the BLAST operation. The function records the best matches from a file 
for a specific sequence in a new file in alphabetical order. 



## Examples

`run_dna_rna_tools("AUG", "reverse")` will return `"GUA"`

`filter_fastq("/path_to_data/data_example", "name_outfile", 60, 30)` will return 
filtered data with parametrs gc-content from 0 to 60% and length from 0 to 30 in folder filtered

`convert_multiline_fasta_to_oneline("/path_to_data/data_example.fasta")` will return correct fasta file

`parse_blast_output("/path_to_data/data_example", "/path_to_result_function/result_file")` will return 
list with the best score from data_example file



## Acknowledgements

I would like to express my gratitude to my nerve cells for their patience and work and to the 
following people: Nikita Vaulin (https://github.com/nvaulin) for understandable lectures, 
Anton Sidorin (https://github.com/SidorinAnton) for detailed seminars and Artyom Ershov 
(https://github.com/iam28th) for edits that allowed us to significantly improve the code

