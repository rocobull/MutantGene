# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 18:59:24 2020

@author: Roberto Bullitta
"""

import MutantGene as mg

maf = None #DEFINE THE DIRECTORY OF A MAF FILE (downloadable here - https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22maf%22%5D%7D%7D%5D%7D)
goa = None #DEFINE THE DIRECTORY OF A GOA FILE (downloadable here - ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/)


# This file serves to show how each function in this library works and possible ways they interact with each other!


# NOTE: The results of each main function will be saved in text files using the functions "list_to_file" (to save lists)
# and "dict_to_file" (to save dictionaries) for a better understanding of the function's purpose.



# We begin by retrieving the information regarding 100 samples from a MAF file (defined at the start of this file)
# starting from the 40th sample. We'll retrieve the information using the "mutation_samples" function and can
# also rely of "maf_index" to decide which indexes retrieve what kind of information.



index = mg.maf_index(maf)
mg.dict_to_file(index, "maf_index.txt")

print(index)

# Now knowing the attribute parameter titles and their indexes, lets pick out HGVSc (34),
# HGVSp (35), VARIANT_CLASS (95), HUGO_Symbol (0) and IMPACT (93).


raw_samples = mg.mutation_samples(maf, 100, 39, [93,0,34,35,95])
mg.list_to_file(raw_samples, "mutation_samples.txt")

#print(raw_samples)

    
# The information retrieved can now be filtered with the auxiliary function "filter_list". 
# In this case, we'll filter out all mutations with an impact considered "LOW".


samples = mg.filter_list(raw_samples, ["LOW"], "n")
mg.list_to_file(samples, "filtered_samples.txt")      #77 samples



# Next, GO IDs need to be generated. The GOA file defined at the start of this file will be used
# to retrieve a dictionary with all it's stored genes as keys and their associated GO IDs as values
# using the "all_goa_id" function.

            
all_ids = mg.all_goa_id(goa)
mg.dict_to_file(all_ids, "all_goa_id.txt")


# With all GO IDs in the GOA file retrieved, we can use the genes from the collected MAF file samples
# to collect the GO IDs associated to them. Because the "get_goa_id" function requires a gene list for
# this retrieval, we need to first convert the list mutation sample information into a list with only
# the genes (contained in index 1 of each list previously created). 


genes = mg.get_genes(samples, 1)     #76 genes (meaning 1 appeared twice in the samples collected)
mg.list_to_file(genes, "genes.txt")
               


# Now we are ready to collect the gene's GO IDs using the "get_goa_id" function.

        
gene_ids = mg.get_goa_id(genes, all_ids)     #72 genes (not all genes are stored in the GOA file of choice)
mg.dict_to_file(gene_ids, "get_goa_id.txt")

# Here is a good opportunity to use the "filter_dict" function. A text file containing
# the GO IDs that are "Child Terms" to the GO ID of the cytoplasm (meaning they are present
# or their function is in some way related to elements of the cytoplasm) has been created, 
# and can be converted to a list using the "file_to_list" function and used to filter out 
# genes whose products don't interact with the cytoplasm.

cyto_ids = mg.file_to_list("cytoplasm_ids.txt", "y")
cyto_genes = mg.filter_dict(gene_ids, cyto_ids, "v", "y")    #24 genes
mg.dict_to_file(cyto_genes, "cyto_genes.txt")




# For the final steps, we will now convert our GO IDs into Terms to better understand
# the function of the collected gene products. For this we need an OBO file. If no
# parameter is specified, the function "get_ontologies" will connect to the URL of
# a general OBO file with the ontologies of all known and used GO IDs/Terms. In 
# this example, this is what we will use.
            

all_terms = mg.get_ontologies()
mg.dict_to_file(all_terms, "all_terms.txt")

# (In case the user wishes to use a more specific OBO file, they can be found here - http://geneontology.org/docs/download-ontology/)


# Finally, we will convert the GO IDs of the genes that were selected by using the
# "id_to_term" function that uses the previously created dictionary of terms for this.
# This will make the characterization of each gene product much easier.


gene_terms = mg.id_to_term(cyto_genes, all_terms)
mg.dict_to_file(gene_terms, "gene_terms.txt")





# NOTE: The "file_to_dict" wasn't used in this example, but for demonstration purposes,
# we can test it now by using the file just saved ("gene_terms.txt") to create a
# dictionary.


terms_dict = mg.file_to_dict("gene_terms.txt")




