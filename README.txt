----------# MUTANTGENE #----------

MutantGene is a Python library containing 13 functions that has as its primary objective the 
characterization of genes that have suffered mutations detected via alignments of tumour sample
reads. This way, through posterior analysis, it is possible to analyse each affected gene product's
functions in the cell and attempt to understand the role each of them played in neoplasm formation.

This library extracts variant information from MAF files (stored in the Genomic Data Commons data
repository) and uses the retrieved HUGO gene symbols to locate their associated Gene Ontology (GO)
IDs organized in Gene Ontology Annotation (GOA) files. From there, it converts the GO IDs into GO
terms, stored in OBO files. The end result is a dictionary containing the desired gene symbols as
keys and lists of all their retrieved GO terms as values, linking each gene to its functions and/or
locations of action in the cell.

The hope is that this library can assist doctors and researchers in a more accurate cancer diagnosis
and adequate treatment applications. This library is in its early phases and the hope is that, with
the help of a committed community, it can reach a state where it can be easily applied to scientific
research and patient diagnosis. A suggested rout for improvement is a more practical gene grouping
to isolate higher impact oncogenes, and possibly features to better characterize the actual mutations
themselves to better predict their impact on gene functions.

This library was created as an end-of-bachelor project in the course "Cellular and Molecular Biology"
(University: Universidade Nova de Lisboa, Faculty: Faculdade de Ciências e Tecnologias (Portugal)).

A special thanks to Professor Ludwig Krippahl for guiding me through this project and giving me advice
on how to make this library the best it could be in the aount of time given.



----------# Functions #----------

The functions are divided into 3 steps: (STEP 1) deals with extracting variant information from the
selected MAF file; (STEP 2) deals with retrieving the GO IDs associated to each desired gene from
the selected GOA file; (STEP 3) deals with converting the collected GO IDs into their corresponding
GO terms using an OBO file.



(STEP 1)

# mutation_samples(maf_file, limit = None, start = 0, index_list = [0,34,35,95,93]):

- Extract sample information from a MAF file.
    
    PARAMETERS:
    - maf_file (str): A MAF file directory.
    - limit (int): A desired number of samples to be processed (optional).
    - start (int): A desired starting index (optional).
    - index_list (list): A list of indexes to use for information retrieval (optional). If no list is specified,
    returns a list of lists containing the gene names (HUGO_symbol), the gene mutation suffered (HGVSc), 
    the effect on the protein produced (HGVSp), the type of mutation (VARIANT_CLASS) and how big is the predicted
    impact on protein viability (IMPACT).
    
    RETURN:
    A list of each processed sample information lists.



# maf_index(maf_file):

- Index each parameter of a MAF file.
    
    PARAMETERS:
    - maf_file (str): A MAF file directory.
    
    RETURN:
    A dictionary of attribute titles as keys and their corresponding indexes as values.



(STEP 2)

# get_genes(samples, gene_index = 0):

- Filter gene names from sample information lists.
    
    PARAMETERS:
    - samples (list): A list of sample information lists.
    - gene_index (int): The index where the gene names are stored (optional).
    
    RETURN:
    A set of unique gene names collected from the list of samples.



# all_goa_id(goa_file):

- Retrieves all genes and corresponding GO IDs from a GOA file.
    
    PARAMETERS:
    - goa_file (str): A GOA file directory.
    
    RETURN: 
    A dictionary containing all the gene names from the  selected GOA file as 
    keys and a set of their respective GO IDs as values.



# get_goa_id(genes, id_dict):

- Obtain the specified genes' GO IDs.
    
    PARAMETERS:
    - genes (list): A list of gene names.
    - id_dict (dict): A dictionary of gene names as keys and a list of GO IDs as values
    
    RETURN: 
    A dictionary containing the selected gene names as keys and set of
    their associated GO IDs as values.



(STEP 3)

# get_ontologies(obo_file = None):

- Extract all GO IDs and their corresponding GO terms from an OBO file.
    
    PARAMETERS:
    - obo_file (str): An optional text file name containing ontologies in OBO format
    or an OBO URL (optional).
    
    RETURN:
    A dictionary with all GO IDs as keys and their corresponding terms as values.



# id_to_term(gene_ids, all_go_ids):

- Get the desired gene's GO terms through their GO IDs.
    
    PARAMETERS:
    - gene_ids (dict): A dictionary with gene names as keys and a list of their associated
    GO IDs as values
    - all_go_ids (dict): A dictionary of GO IDs as keys and their corresponding terms as values.
    
    RETURN:
    A dictionary of gene names as keys and the terms of their corresponding GO IDs as values.



(UNIVERSAL AUXILIARY FUNCTIONS)

# filter_list(raw_list, filter_list, include = "y"):

- Filters a chosen list.
    
    PARAMETERS:
    - raw_list (list): A list to be filtered
    - filter_list (list): A list of elements to consider for filtering
    - include (str): A character to decide if the list to be filtered must contain
    ("y") or not ("n") the given elements in the filter list (optional).
    
    RETURN:
    A filtered version of the list according to the parameters given.



# list_to_file(list_to_save, file_name):

- Save a list to a file.
    
    PARAMETERS:
    - list_to_save (list): A list to be saved.
    - file_name (str): A text file name.
    
    RETURN:
    Saves the list elements in individual lines in a text file with the given name.



# file_to_list(file_name, extend = "n"):

- Converts lines in a file to a list.
    
    PARAMETERS:
    - file_name (str): A file name
    - extend (str): A character that will define if the final list will be a list
    of lists ("n") of if all elements will be extended into only 1 list (not "n").
    
    RETURN: 
    A list of each line in the file.



# filter_dict(raw_dict, filter_list, key_value = "k", include = "y"):

- Filters a chosen dictionary.
    
    PARAMETERS:
    - raw_dict (dict): A dictionary to be filtered
    - filter_list (str): A list of elements to consider for filtering
    - key_value (str): A character to decide if the elements to be compared are the keys
    ("k") or values (not "k")
    - include (str): A character to decide if the dictionary to be filtered must contain
    ("y") or not ("n") the given elements in the filter list.
    
    RETURN: 
    A filtered dictionary according to the parameters given.



# dict_to_file(dictionary, file_name):

- Saves a dictionary to a file.
    
    PARAMETERS:
    - dictionary (dict): A dictionary to be saved.
    - file_name (str): A text file name.
    
    RETURN:
    Saves the keys and their corresponding values (seperated by ": ") in individual
    lines of a created file with the given name.



# file_to_dict(file_name):

- Converts lines in a file to a dictionary.
    
    PARAMETERS:
    - file_name (str): A file name where each line contains a key seperated from its values by ": ".
    
    RETURN:
    A dictionary with all keys and their corresponding values saved in the file.



----------# Licence #----------

This library is under the GPL (https://en.wikipedia.org/wiki/GNU_General_Public_License)




