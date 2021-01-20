# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 18:59:24 2020

@author: Roberto Bullitta
"""

import requests


maf = 'C:\\Users\\Roberto Bullitta\\Desktop\\PROJETO BCM\\Estágio\\MAF\\MAF_file_1.txt'
goa = 'C:\\Users\\Roberto Bullitta\\Desktop\\PROJETO BCM\\Estágio\\goa_human.gaf'



# ------------------------------ STEP 1 (Mutation samples) ------------------------------

def mutation_samples(maf_file, limit = None, start = 0, index_list = [0,34,35,95,93]):
    """Extract sample information from a MAF file.
    
    PARAMETERS:
    - maf_file (str): A MAF file directory.
    - limit (int): A desired number of samples to be processed (optional).
    - start (int): A desired starting index (optional).
    - index_list (list): A list of indexes to use for information retrieval (optional). If no list is specified,
    returns a list of lists containing the gene names (HUGO_symbol), the gene mutation suffered (HGVSc), 
    the effect on the protein produced (HGVSp), the type of mutation (VARIANT_CLASS) and how big the impact 
    on protein viability is (IMPACT).
    
    RETURN:
    A list of each processed sample information lists.
    """

    # Open MAF file and grab all lines
    file = open(maf_file)
    info = file.readlines()
    
    # Delete the initial annotation lines that start with "#".
    for _ in range(len(info)):
        if info[0].startswith("#"):
            del info[0]
        else:
            break
    
    # Redifine "limit" to be the full length of the file in case "limit" is None
    if limit == None:
        limit = len(info) - 1 - start
    
    # Redifine "limit" in case it is larger than the number of lines left until the end of the file.
    if len(info) - 1 - start < limit:
        limit = len(info) - 1 - start
    

    # Extract information stored in the chosen IDs.
    all_samples = []
    for i in range(limit):
        sample = []
        extract = info[i+1+start].strip().split("\t")
        
        for i in index_list:
            sample.append(extract[i])
            
        all_samples.append(sample)
        

    return all_samples





# ((((((((((AUXILIARY FUNCTIONS)))))))))):

def maf_index(maf_file):
    """Index each parameter of a MAF file.
    
    PARAMETERS:
    - maf_file (str): A MAF file directory.
    
    RETURN:
    A dictionary of attribute titles as keys and their corresponding indexes as values.
    """
    
    # Open MAF file and grab attribute titles line
    file = open(maf_file)
    info = file.readlines()
    
    # Delete the initial annotation lines that start with "#".
    for _ in range(len(info)):
        if info[0].startswith("#"):
            del info[0]
        else:
            break
    
    titles = info[0].strip().split("\t")

    # Save titles with corresponding indexes 
    all_indexes = {}

    for i in range(len(titles)):
        all_indexes[titles[i]] = i
        
    return all_indexes





# ------------------------------ STEP 2 (GO IDs) ------------------------------

def all_goa_id(goa_file):
    """Retrieves all genes and corresponding GO IDs from a GOA file.
    
    PARAMETERS:
    - goa_file (str): A GOA file directory.
    
    RETURN: 
    A dictionary containing all the gene names from the  selected GOA file as 
    keys and a set of their respective GO IDs as values.
    """
    
    # Retrieve GOA file information
    file = open(goa_file)
    info = file.readlines()
    
    # Delete the initial annotation lines that start with "!".
    for _ in range(len(info)):
        if info[0].startswith("!"):
            del info[0]
        else:
            break
    
    # Store read genes in order to add GO IDs to their respective sets
    check_repeats = []
    
    # Prepare a list of synonyms to attribute the same GO IDs to as the main gene
    main_gene_old = None
    gene_syn_old = None
    first_time = "y"
    
    # Retrieve all genes (including synonyms) and respective GO IDs in the selected GOA file.
    complete_info = {}
    
    for i in info:
        line = i.split("\t")
        main_gene_new = line[2]
        go_id = line[4]
        gene_syn_new = line[10].split("|")  # Where synonyms are located
        
        if main_gene_new not in check_repeats:
            if first_time == "n":
                for syn in gene_syn_old:
                    complete_info[syn] = complete_info[main_gene_old]
                    
            check_repeats.append(main_gene_new)
            complete_info[main_gene_new] = set([go_id])
            
            first_time = "n"
            gene_syn_old = gene_syn_new
            main_gene_old = main_gene_new
            
        else:
            complete_info[main_gene_new].add(go_id)
    
    for syn in gene_syn_old:
        complete_info[syn] = complete_info[main_gene_old]   # Add the synonyms of the last gene to the dictionary
            
    return complete_info



def get_goa_id(genes, id_dict):
    """Obtain specified genes' GO IDs.
    
    PARAMETERS:
    - genes (list): A list of gene names.
    - id_dict (dict): A dictionary of gene names as keys and a list of GO IDs as values
    
    RETURN: 
    A dictionary containing the selected gene names as keys and set of
    their associated GO IDs as values.
    """
    
    result = {}

    for g in genes:
                
        if g in id_dict.keys():
            result[g] = id_dict[g]    

    return result
        

    
def get_genes(samples, gene_index = 0):
    """Filter gene names from sample information lists.
    
    PARAMETERS:
    - samples (list): A list of sample information lists.
    - gene_index (int): The index where the gene names are stored (optional).
    
    RETURN:
    A set of unique gene names collected from the list of samples. 
    """
    
    genes = set()
    
    for s in samples:
        genes.add(s[gene_index])
    
    return genes





# ------------------------------ STEP 3 (GO Terms) ------------------------------

def get_ontologies(obo_file = None):
    """Extract all GO IDs and their corresponding GO terms from an OBO file.
    
    PARAMETERS:
    - obo_file (str): An optional text file name containing ontologies in OBO format
    or an OBO URL (optional).
    
    RETURN:
    A dictionary with all GO IDs as keys and their corresponding terms as values.
    """
    
    # Get ontology lines
    if obo_file == None:
        obo = requests.get("http://current.geneontology.org/ontology/go.obo") 
        info = obo.text.split("\n\n")
    
    else:                                  # In case another OBO URL is given
        if obo_file.startswith("http"):
            obo = requests.get(obo_file)
            info = obo.text.split("\n\n")
        else:                              # In case an OBO file is given
            obo = open(obo_file).read()
            info = obo.split("\n\n")
        
    # Find and store each GO ID in the OBO file as keys and their associated term as values
    id_term = {}
    
    for i in info:
        if i.startswith("[Term]"):
            lines = i.split("\n")
            id_term[lines[1][4::]] = lines[2][6::]   #"[4::]" to skip "id: " and "[6::]" to skip "name: "
    
    return id_term
            


def id_to_term(gene_ids, all_go_ids):
    """Get the desired gene's GO terms through their GO IDs.
    
    PARAMETERS:
    - gene_ids (dict): A dictionary with gene names as keys and a list of their associated
    GO IDs as values
    - all_go_ids (dict): A dictionary of GO IDs as keys and their corresponding terms as values.
    
    RETURN:
    A dictionary of gene names as keys and the terms of their corresponding GO IDs as values.  
    """

    # Prepare a dictionary to store the GO terms converted from their GO IDs
    gene_terms = {}
    for g in gene_ids:
        terms = []
        
        for i in gene_ids[g]:
        
            if i in all_go_ids.keys():
                terms.append(all_go_ids[i])
        
        gene_terms[g] = terms
    
    return gene_terms





#-----------------------------------------UNIVERSAL AUXILIARY FUNCTIONS--------------------------------------

def filter_list(raw_list, filter_list, include = "y"):
    """Filters a chosen list.
    
    PARAMETERS:
    - raw_list (list): A list to be filtered
    - filter_list (list): A list of elements to consider for filtering
    - include (str): A character to decide if the list to be filtered must contain
    ("y") or not ("n") the given elements in the filter list (optional).
    
    RETURN:
    A filtered version of the list according to the parameters given.
    """
    
    # Choose elements to remove from raw list according the function parameters given.
    new_list = raw_list.copy()
    to_remove = []
    

    for raw in raw_list:
        if type(raw) is list:    # Check if there are lists within the raw list
            for r in raw:
                check = 0    #To make sure lists to be excluded are completely searched
                counter = 0
                
                for f in filter_list:
                    counter = counter + 1
                    
                    if include == "y":                            # Include lists of raw list with filter list elements
                        if f in r:
                            check = 1
                            break
                        
                        else:
                            if counter == len(filter_list):
                                to_remove.append(raw)
                                check = 1
                                continue
                            
                    else:                                         # Exclude lists of raw list with filter list elements
                        if f in r:
                            to_remove.append(raw)
                            check = 1
                            break
                        
                        else:
                            continue
                
                if check == 1:
                    break
    
        else:
            counter = 0
            for f in filter_list:
                counter = counter + 1
                
                if include == "y":                                # Include elements of raw list with filter list elements
                    if f in raw:
                        break
                        
                    else:
                        if counter == len(filter_list):
                            to_remove.append(raw)
                            continue
                        
                else:                                             # Exclude elements of raw list with filter list elements
                    if f in raw:
                        to_remove.append(raw)
                        break
        
    # Filter the raw list to obtain final filtered list
    if to_remove != []:
        for t in to_remove:
            new_list.remove(t)
            
    return new_list
                
            

def list_to_file(list_to_save, file_name):
    """Save a list to a file.
    
    PARAMETERS:
    - list_to_save (list): A list to be saved.
    - file_name (str): A text file name.
    
    RETURN:
    Saves the list elements in individual lines in a text file with the given name.
    """
    
    save_file = open(file_name, "w")
    
    for l in list_to_save:
        save_file.write(str(l)+"\n")
    
    save_file.close()



def file_to_list(file_name, extend = "n"):
    """Converts lines in a file to a list.
    
    PARAMETERS:
    - file_name (str): A file name
    - extend (str): A character that will define if the final list will be a list
    of lists ("n") of if all elements will be extended into only 1 list (not "n").
    
    RETURN: 
    A list of each line in the file.
    """
    
    line_list = []
    lines = open(file_name).readlines()
    
    for l in lines:
        l = l.strip().strip("][").strip("}{")         # Turn list or set strings into lists
        l = l.replace("'","").replace(" ","").split(",")  # Remove extra quote marks and spaces         
        
        for ix in range(len(l)):
            try:
                l[ix] = int(l[ix])
            except:
                None
        
        if extend == "n":
            line_list.append(l)
        else:
            line_list.extend(l)

    return line_list



def filter_dict(raw_dict, filter_list, key_value = "k", include = "y"):
    """Filters a chosen dictionary.
    
    PARAMETERS:
    - raw_dict (dict): A dictionary to be filtered
    - filter_list (str): A list of elements to consider for filtering
    - key_value (str): A character to decide if the elements to be compared are the keys
    ("k") or values (not "k")
    - include (str): A character to decide if the dictionary to be filtered must contain
    ("y") or not ("n") the given elements in the filter list.
    
    RETURN: 
    A filtered dictionary according to the parameters given.
    """
    
    # Create a new dictionary to be filtered and a set to use to filter genes
    new_dict = raw_dict.copy()
    to_remove = set()
    for k in new_dict:
        if key_value == "k":               # Filter genes by comparing the genes themselves to the elements of the filter list
            counter = 0
            
            for f in filter_list:
                counter = counter + 1
                
                if include == "y":                                      # Include key-value pairs with specified keys
                    if f in k:
                        break
                    
                    else:                                               # Exclude key-value pairs without specified keys
                        if counter == len(filter_list):
                            to_remove.add(k)
                
                else:
                    if f in k:                                          # Exclude key-value pairs with specified keys
                        to_remove.add(k)
                        break
            
        else:                      # Filter genes according to the values of the dictionary
            if type(new_dict[k]) is list or type(new_dict[k]) is set:   # In case the values are lists
                check = 0
                counter_f = 0
                counter_v = 0
                for f in filter_list:
                    counter_f = counter_f + 1
                    for v in new_dict[k]:
                        counter_v = counter_v + 1
                        if include == "y":
                            if f in v:                             # Include key-value pairs with specified elements in value lists
                                check = 1
                                break
                            
                            else:                                  # Exclude key-value pairs without specified elements in value lists
                                if counter_f == len(filter_list) and counter_v == len(new_dict[k]):
                                    to_remove.add(k)
                                    check = 1
                            
                        else:                                      # Exclude key-value pairs with specified elements in value lists
                            if f in v:
                                to_remove.add(k)
                                check = 1
                                break
                    if check == 1:
                        break
            else:
                counter = 0
                for f in filter_list:
                    counter = counter + 1
                            
                    if include == "y":
                        if f in raw_dict[k]:                     # Include key-value pairs with a specified value
                            break
                            
                        else:                                    # Exclude key-value pairs without a specified value
                            if counter == len(filter_list):
                                to_remove.add(k)
                            
                    else:                                        # Exclude key-value pairs with a specified value
                        if f in raw_dict[k]:
                            to_remove.add(k)
                            break
                        
    # Filter the raw list to obtain final filtered list
    
    if to_remove != []:
        for t in to_remove:
            del new_dict[t]
            
    return new_dict



def dict_to_file(dictionary, file_name):
    """Saves a dictionary to a file.
    
    PARAMETERS:
    - dictionary (dict): A dictionary to be saved.
    - file_name (str): A text file name.
    
    RETURN:
    Saves the keys and their corresponding values (seperated by ": ") in individual
    lines of a created file with the given name.
    """
    
    save_file = open(file_name, "w")
    
    for d in dictionary:
        save_file.write(str(d)+": "+str(dictionary[d])+"\n")
    
    save_file.close()



def file_to_dict(file_name):
    """Converts lines in a file to a dictionary.
    
    PARAMETERS:
    - file_name (str): A file name where each line contains a key seperated from its values by ": ".
    
    RETURN:
    A dictionary with all keys and their corresponding values saved in the file.
    """
    
    line_dict = {}
    lines = open(file_name).readlines()
    
    for l in lines:
        k,v = l.strip().split(": ")
        if (v[0] == "[" and v[-1] == "]") or (v[0] == "{" and v[-1] == "}"):   # Check to see iv values are lists
            v = v.strip("][").strip("}{")                                      # Turn list strings into lists
            v = v.replace("'","").replace(" ","").split(",")                   # Remove extra quote marks and spaces
        else:
            try:
                v = int(v)
            except:
                None
            
        line_dict[k] = v
    
    return line_dict