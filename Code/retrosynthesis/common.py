import random
import pandas as pd
import json
import os
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem import AllChem
from tqdm import tqdm
from Bio import SeqIO
import pandas as pd
import itertools
import pickle
from multiprocessing import Pool
import pubchempy as pcp
import time

from concurrent.futures import ProcessPoolExecutor
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDraw2DCairo
from itertools import islice
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import cobra
from cobra import Metabolite,Reaction,Gene
from ast import literal_eval
import re
import numpy as np
import multiprocessing as mp
from functools import partial
import  ast
import multiprocessing

def normalize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        # mol = Chem.RemoveHs(mol)
        if mol is not None:
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)###False
            # canonical_smiles = Chem.MolToSmiles(mol,True)###False
            return canonical_smiles
        else:
            pass
    except:
        pass


def get_all_smiles_in_model(model_path,ymdb_path):
    model = pd.read_csv(model_path)
    model = model.dropna(subset=['standard_smiles'])
    model_smiles = model['standard_smiles'].to_list()
    ymdb = pd.read_excel(ymdb_path)
    ymdb = ymdb.dropna(subset=['standard_smiles'])
    ymdb = ymdb[ymdb['in_model'] == 1]
    ymdb_smiles = ymdb['standard_smiles'].to_list()
    total_smiles = model_smiles + ymdb_smiles
    # Normalize the SMILES strings in model_smiles
    total_smiles = [normalize_smiles(x) for x in total_smiles]
    # Remove duplicates
    total_smiles = list(set(total_smiles))
    return total_smiles
def get_GEM_all_smiles(model_path):
    model = pd.read_csv(model_path)
    model = model.dropna(subset=['standard_smiles'])
    model_smiles = model['standard_smiles'].to_list()
    # Normalize the SMILES strings in model_smiles
    model_smiles = [normalize_smiles(x) for x in model_smiles]
    # Remove duplicates
    model_smiles = list(set(model_smiles))
    return model_smiles



def check_rule(rules):
    Reactant = rules.split('>>')[0]
    Product = rules.split('>>')[1]
    Reactant = Reactant.replace('c','C').replace('[','').replace(']','').replace('@','').replace('+','').replace('-','').replace('(','').replace(')','').replace('.','').replace('>','').replace('=','')
    Product = Product.replace('c','C').replace('[','').replace(']','').replace('@','').replace('+','').replace('-','').replace('(','').replace(')','').replace('.','').replace('>','').replace('=','')
    #Use Counter to count the number of characters that appear
    Reactant_count = Counter(Reactant)
    Product_count = Counter(Product)
    if Product_count['C'] == Reactant_count['C'] and Reactant_count['P']==Product_count['P']:
        return('blance')
    else:return('unblance')
def get_most_similar_smiles(retrosys_smiles_calculate_similarity_pd,num):
    retrosys_smiles_calculate_similarity_pd[['smiles_in_mets_total_smiles', 'scores']] = retrosys_smiles_calculate_similarity_pd[['smiles_in_mets_total_smiles', 'scores']].applymap(lambda x: x[:num] if len(x)>num else x)
    # retrosys_smiles_calculate_similarity_pd = retrosys_smiles_calculate_similarity_pd.applymap(lambda x: x[:num] if len(x)>num else x) # Add closing parenthesis ")" here
    return retrosys_smiles_calculate_similarity_pd
def get_calculate_similarity(product_smiles,retrosys_smiles_calculate_similarity_pd):
    smiles_in_mets_total_smiles = retrosys_smiles_calculate_similarity_pd[retrosys_smiles_calculate_similarity_pd['rules_smiles']== product_smiles]['smiles_in_mets_total_smiles'].to_list()[0]
    scores = retrosys_smiles_calculate_similarity_pd[retrosys_smiles_calculate_similarity_pd['rules_smiles']== product_smiles]['scores'].to_list()[0]
    return(smiles_in_mets_total_smiles,scores)

def get_GEM_all_mebaolite_with_smiles(model_path):
    if 'mat' in model_path:
        model = cobra.io.load_matlab_model(model_path)
        yeast_total_smiles = []
        for i in model.metabolites:
            if 'SMILES' in i.annotation:
                try: 
                    smiles_tmp = i.annotation['SMILES'][0]
                    yeast_total_smiles.append(normalize_smiles(smiles_tmp))
                except:pass
            # if hasattr(i, 'smiles'):
        yeast_total_smiles = list(set(yeast_total_smiles))    
        return(yeast_total_smiles)
    else: print('error model file')
def get_total_smiles_lipid(ymdb_path,model_path, target_smiles, output_file):
    # Load the model from a CSV file
    model = pd.read_csv(model_path)

    model = model.dropna(subset=['standard_smiles'])

    model_smiles = model['standard_smiles'].to_list()
    model_smiles = list(set(model_smiles))
    ymdb = pd.read_excel(ymdb_path)
    ymdb = ymdb.dropna(subset=['standard_smiles'])
    ymdb = ymdb[ymdb['in_model'] == 1]
    ymdb_smiles = ymdb['standard_smiles'].to_list()
    ymdb_smiles = list(set(ymdb_smiles))
    # Combine model_smiles with target_smiles
    total_smiles = model_smiles + target_smiles + ymdb_smiles

    # Normalize the SMILES strings in total_smiles
    total_smiles = [normalize_smiles(x) for x in total_smiles]

    # total_smiles = list(set(total_smiles))

    print('total_smiles', len(total_smiles))

    dump_file(total_smiles, output_file)
def get_total_smiles(model_path, target_smiles, output_file):
    # Load the model from a CSV file
    model = pd.read_csv(model_path)

    model = model.dropna(subset=['standard_smiles'])

    model_smiles = model['standard_smiles'].to_list()
    model_smiles = list(set(model_smiles))
    # Combine model_smiles with target_smiles
    total_smiles = model_smiles + target_smiles 

    # Normalize the SMILES strings in total_smiles
    total_smiles = [normalize_smiles(x) for x in total_smiles]

    # total_smiles = list(set(total_smiles))

    print('total_smiles', len(total_smiles))

    dump_file(total_smiles, output_file)

def get_target_smiles(target_smiles, output_file):
    # Normalize the SMILES strings in target_smiles
    target_smiles = [normalize_smiles(x) for x in target_smiles]
    target_smiles = [x for x in target_smiles if x is not None]
    # target_smiles = list(set(target_smiles))
 
    print('target_smiles:', len(target_smiles))

    dump_file(target_smiles, output_file)
    return target_smiles
def process_yeast_smiles(file_path):
    # Load smiles from file
    total_smiles = load_pickle(file_path)
    print('total_smiles', len(total_smiles))
    total_smiles = [smile for smile in total_smiles if  smile is not None]
    # Filter out smiles containing '.'
    total_smiles = [smile for smile in total_smiles if '.' not in smile]
    print('total_smiles', len(total_smiles))

    # Filter out smiles that do not contain 'C' or 'c'
    total_smiles = [smile for smile in total_smiles if 'C' in smile or 'c' in smile]
    print('total_smiles', len(total_smiles))

    # Normalize the smiles and remove duplicates
    total_smiles = list(set([normalize_smiles(smile) for smile in total_smiles]))
    print('total_smiles', len(total_smiles))

    return total_smiles


def smiles2inchikey0(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchikey = Chem.MolToInchiKey(mol)
        if inchikey:
            return inchikey.split('-')[0]
        else:
            return smiles 
    except:
        return smiles





def compare_smiles_inchikey(smile1,smile2):
    
    inchikey1 = smiles2inchikey0(smile1)
    inchikey2 = smiles2inchikey0(smile2)
    
    return inchikey1 == inchikey2


def print_first_ten(dictionary):
    print(dict(islice(dictionary.items(), 10)))


def fasta_to_dataframe(fasta_file):
    sequences = []
    genes = []

    # Read the FASTA file and extract sequences and headers
    for record in SeqIO.parse(fasta_file, "fasta"):
        genes.append(record.id)
        sequences.append(str(record.seq))

    # Create a DataFrame from the sequences and headers
    df = pd.DataFrame({'gene': genes, 'sequence': sequences})
    return df


def dump_file(dictionary, filename):
    with open(filename, 'wb') as file:
        pickle.dump(dictionary, file)

def load_pickle(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)
    
def fasta2dataframe(fasta_file):
    records = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        record_data = {
            'ID': record.id,
            # 'Description': record.description,
            'Sequence': str(record.seq)
        }
        records.append(record_data)

    # create DataFrame
    df = pd.DataFrame(records)
    return(df)

def expand_list_column(df, column_name):
    expanded_rows = pd.DataFrame()

    for index, row in tqdm(df.iterrows(),total=len(df)):
        values = row[column_name]
        if isinstance(values, list):
            if len(values) > 1:
                for value in values:
                    new_row = row.copy()
                    new_row[column_name] = value
                    expanded_rows = expanded_rows.append(new_row, ignore_index=True)
            elif len(values) == 1:
                # For lists with length one, convert the list to a single string
                new_row = row.copy()
                new_row[column_name] = values[0]  # Take the first (and only) element as a string
                expanded_rows = expanded_rows.append(new_row, ignore_index=True)
            else:pass
        else:
            expanded_rows = expanded_rows.append(row, ignore_index=True)
    return expanded_rows



def compare_smiles(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 and mol2:
        mol1_fp = Chem.RDKFingerprint(mol1)
        mol2_fp = Chem.RDKFingerprint(mol2)
        similarity = TanimotoSimilarity(mol1_fp, mol2_fp)
        return similarity
    else:
        return 0


def are_atom_counts_equal(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    # Check if molecule creation was successful
    if mol1 is None or mol2 is None:
        return None
    
    atom_count1 = Counter([atom.GetAtomicNum() for atom in mol1.GetAtoms()])
    atom_count2 = Counter([atom.GetAtomicNum() for atom in mol2.GetAtoms()])
    
    # Extract the counts of specific elements (e.g., carbon, nitrogen, oxygen, phosphorus)
    c_count1 = atom_count1[6]  # 6 represents carbon
    c_count2 = atom_count2[6]
    
    n_count1 = atom_count1[7]  # 7 represents nitrogen
    n_count2 = atom_count2[7]
    
    o_count1 = atom_count1[8]  # 8 represents oxygen
    o_count2 = atom_count2[8]
    
    p_count1 = atom_count1[15]  # 15 represents phosphorus
    p_count2 = atom_count2[15]
    
    # Check if the counts of carbon, nitrogen, oxygen, and phosphorus atoms in both molecules are equal
    if c_count1 == c_count2 and n_count1 == n_count2 and o_count1 == o_count2 and p_count1 == p_count2:
        return 1
    else:
        return 0

def calculate_similarity(smiles1, smiles2, fingerprint_type="MACCS"):
    # Convert the SMILES strings to RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Check if the conversion to molecule objects was successful
    if mol1 is None or mol2 is None:
        raise ValueError("Unable to parse the provided SMILES strings")

    if fingerprint_type == "MACCS":
        fingerprint1 = AllChem.GetMACCSKeysFingerprint(mol1)
        fingerprint2 = AllChem.GetMACCSKeysFingerprint(mol2)
    else:
        raise ValueError("Unsupported fingerprint type")

    # Calculate the Tanimoto similarity
    similarity = FingerprintSimilarity(fingerprint1, fingerprint2)

    return similarity

def draw_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    # Convert the RDKit image object to a format that Matplotlib can use, and set the DPI and image size
    img_array = Draw.MolToImage(mol, size=(400, 300)).convert('CMYK')

    fig, ax = plt.subplots(figsize=(1,1), dpi=300)
    # Display the image
    ax.imshow(img_array)
    ax.axis('off')  # Turn off the axes
    # Save the image as a PDF file
    plt.show()
def are_dicts_equal(dict1, dict2):
    # Determine whether the key value pairs of the two dictionaries are exactly the same and are not sensitive to the order.
    return sorted(dict1.items()) == sorted(dict2.items())

    # # example
    # dict1 = {'a': 1, 'b': 2, 'c': 3}
    # dict2 = {'b': 2, 'a': 1, 'c': 3}
    # dict3 = {'a': 1, 'b': 2, 'c': 4}

    # print(are_dicts_equal(dict1, dict2))  # output True
    # print(are_dicts_equal(dict1, dict3))  # output False

def get_yeast8U(new_met_info_to_GEM_path,rxndb_to_model_total_info_path,yeast8_reaction_in_rxndb_json,yeast870_path,yeast8U_path):
    met_info = pd.read_csv(new_met_info_to_GEM_path)
    print(met_info.shape)
    rxndb_to_model = pd.read_csv(rxndb_to_model_total_info_path)
    print(rxndb_to_model.shape)
    rxndb_to_model['equation_dict'] = rxndb_to_model['equation_dict'].apply(lambda x:literal_eval(x))
    rxndb_to_model = rxndb_to_model[~rxndb_to_model['equation'].str.startswith(' =>')]
    rxndb_to_model = rxndb_to_model[~rxndb_to_model['equation'].str.endswith(' => ')]
    rxndb_to_model = rxndb_to_model.drop_duplicates(['equation'], keep='first')
    print(rxndb_to_model.shape)

    with open(yeast8_reaction_in_rxndb_json, 'r') as file:
        yeast8_reaction_in_rxndb_dict = json.load(file)
        
    # yeast8_reaction_in_rxndb = list(set(yeast8_reaction_in_rxndb.values()))
    yeast8_reaction_in_rxndb = list(set([x for lst in yeast8_reaction_in_rxndb_dict.values() for x in lst]))
    rxndb_to_model = rxndb_to_model[~rxndb_to_model['NO'].isin(yeast8_reaction_in_rxndb)]
    print(rxndb_to_model.shape)

    model = cobra.io.load_yaml_model(yeast870_path)
    ## new metabolite
    metabolites = {}
    for index, row in met_info.iterrows():
        met_id = row['ID']
        met_compartment = row['compartment']
        met_smiles = row['new_met_smiles']
        metabolites[met_id] = Metabolite(met_id, formula='', name='', compartment=met_compartment)
    ## new reaction
    for index,row in tqdm(rxndb_to_model.iterrows(),total=len(rxndb_to_model)):
        reaction = Reaction(row['NO'])
        reaction.name = row['NO']
        reaction.subsystem = ''
        reaction.lower_bound = 0.  # This is the default
        reaction.upper_bound = 1000.  # This is the default

        reactant_met_num = row['equation_dict']
        for met_id, coeff in reactant_met_num.items():
            met = metabolites.get(met_id)
            if met:
                reaction.add_metabolites({met: coeff})
        if  row['GPR']:
            reaction.gene_reaction_rule = row['GPR'].replace('(','').replace(')','')
        model.add_reactions([reaction])

    # cobra.io.save_matlab_model(model, yeast8U_path)
    cobra.io.save_yaml_model(model, yeast8U_path)

def compare_atoms_between_formula_smiles(smiles, formula):
    try:
        mol = Chem.MolFromSmiles(smiles)
        atom_count = Counter([atom.GetAtomicNum() for atom in mol.GetAtoms()])
        
        c_count = atom_count[6]  # Carbon
        o_count = atom_count[8]  # Oxygen
        n_count = atom_count[7]  # Nitrogen
        p_count = atom_count[15] # Phosphorus
        
        # Extract counts from formula
        carbon_count_formula_match = re.search(r'C(\d*)', formula)
        oxygen_count_formula_match = re.search(r'O(\d*)', formula)
        nitrogen_count_formula_match = re.search(r'N(\d*)', formula)
        phosphorus_count_formula_match = re.search(r'P(\d*)', formula)
        
        carbon_count_formula = int(carbon_count_formula_match.group(1) or 1) if carbon_count_formula_match else 0
        oxygen_count_formula = int(oxygen_count_formula_match.group(1) or 1) if oxygen_count_formula_match else 0
        nitrogen_count_formula = int(nitrogen_count_formula_match.group(1) or 1) if nitrogen_count_formula_match else 0
        phosphorus_count_formula = int(phosphorus_count_formula_match.group(1) or 1) if phosphorus_count_formula_match else 0
        
        # compare counts
        if (c_count == carbon_count_formula and
            o_count == oxygen_count_formula and
            n_count == nitrogen_count_formula and
            p_count == phosphorus_count_formula):
            return 1
        else:
            return 0
    except:
        return 0
    
def get_gene2ec_dict_clean(sce_gene_clean_ec):
    with open(sce_gene_clean_ec, 'r') as file:
        csv_data = file.read()
    rows = csv_data.split('\n')
    gene2ec_dict = {}
    for row in rows:
        columns = row.split(',')
        key = columns[0]
        values = [v.split('/')[0].replace('EC:','') for v in columns[1:]]
        if len(values)>0:
            gene2ec_dict[key] = values
    for key, values in gene2ec_dict.items():
        gene2ec_dict[key] = list(set([".".join(value.split(".")[:3]) for value in values]))
    return gene2ec_dict

def get_ec2gene_dict_clean(gene2ec_dict):
    clean_gene_list = list(gene2ec_dict.keys())
    ## gene2EC 转换为 ec2gene
    clean_ec2gene_dict = {}
    for gene, ec_list in gene2ec_dict.items():
        for ec in ec_list:
            if ec not in clean_ec2gene_dict:
                clean_ec2gene_dict[ec] = []
            clean_ec2gene_dict[ec].append(gene)
    return clean_ec2gene_dict

def get_gene2ec_dict_DeepEC(DeepEC_path):
    DeepECv2_res = pd.read_csv(DeepEC_path,sep='\t')
    DeepECv2_res = DeepECv2_res.dropna(subset=['prediction'])
    DeepECv2_res = DeepECv2_res[DeepECv2_res['prediction']!='None']
    DeepECv2_res['prediction'] = DeepECv2_res['prediction'].apply(lambda x:x.split(':')[1])
    gene2ec_dict = {}
    for index,row in DeepECv2_res.iterrows():
        if row['sequence_ID'] not in gene2ec_dict:
            gene2ec_dict[row['sequence_ID']] = []
        gene2ec_dict[row['sequence_ID']].append(row['prediction'])
    for key, values in gene2ec_dict.items():
        gene2ec_dict[key] = list(set([".".join(value.split(".")[:3]) for value in values]))

    return gene2ec_dict 


def get_ec2gene_dict_DeepEC(gene2ec_dict):
    # DeepEC_gene_list = list(gene2ec_dict.keys())
    DeepEC_ec2gene_dict = {}
    for gene, ec_list in gene2ec_dict.items():
        for ec in ec_list:
            if ec not in DeepEC_ec2gene_dict:
                DeepEC_ec2gene_dict[ec] = []
            DeepEC_ec2gene_dict[ec].append(gene)
            
    return DeepEC_ec2gene_dict 

def list_to_frequency_dict(input_list):
    frequency_dict = {}
    for item in input_list:
        if item in frequency_dict:
            frequency_dict[item] += 1
        else:
            frequency_dict[item] = 1
    return frequency_dict

def calculate_gene_frequency(yeast):
    gene_frequency_dict = {}
    gene_list = []
    for reaction in yeast.reactions:
        tmp = str(reaction.gpr).replace(')','').replace('(','').replace(' and ',' ').replace(' or ',' ').split(' ')
        gene_list += tmp
    gene_list = [x for x in gene_list if x != '']
    print(len(set(gene_list)))
    gene_frequency_dict = list_to_frequency_dict(gene_list)
    return gene_frequency_dict

def gene_list2_ec_list(gene_lst,gene2ec_dict):
    ec_list = []
    for gene in gene_lst:
        if gene in gene2ec_dict.keys():
            ec_list+=gene2ec_dict[gene]
        else:
            pass
    ec_list = list(set(ec_list))
    return ec_list
def get_exchange_reaction(Target_met,model):
    if Target_met!='':
        exchange_reaction = ''
        for i in model.reactions:
            if Target_met in i.reaction and len(i.metabolites) == 1:
                exchange_reaction = i.id
        if exchange_reaction == '':
            reaction_name = 'DM_'+ Target_met
            reaction = Reaction(reaction_name)
            reaction.name = reaction_name
            reaction.subsystem = ''
            reaction.lower_bound = 0.  # This is the default
            reaction.upper_bound = 1000.  # This is the default

            reactant_met_num = {Target_met:-1}
            for met_id, coeff in reactant_met_num.items():
                met = model.metabolites.get_by_id(met_id)
                if met:
                    reaction.add_metabolites({met: coeff})
            model.add_reactions([reaction])  
            exchange_reaction = reaction_name     
        return(exchange_reaction)
    else:
        pass

def rank_ec_number_dict(ec_numbers):
    ec1 = {}
    ec2 = {}
    ec3 = {}
    ec4 = {}
    ec5 = {}
    ec6 = {}
    ec7 = {}

    # Assign values to different dictionaries based on the key prefix
    for k, v in ec_numbers.items():
        if k.startswith('1.'):
            ec1[k] = v
        elif k.startswith('2.'):
            ec2[k] = v
        elif k.startswith('3.'):
            ec3[k] = v
        elif k.startswith('4.'):
            ec4[k] = v
        elif k.startswith('5.'):
            ec5[k] = v
        elif k.startswith('6.'):
            ec6[k] = v
        elif k.startswith('7.'):
            ec7[k] = v

    # Sort the dictionaries in descending order by value
    ec1 = dict(sorted(ec1.items(), key=lambda x: x[1], reverse=True))
    ec2 = dict(sorted(ec2.items(), key=lambda x: x[1], reverse=True))
    ec3 = dict(sorted(ec3.items(), key=lambda x: x[1], reverse=True))
    ec4 = dict(sorted(ec4.items(), key=lambda x: x[1], reverse=True))
    ec5 = dict(sorted(ec5.items(), key=lambda x: x[1], reverse=True))
    ec6 = dict(sorted(ec6.items(), key=lambda x: x[1], reverse=True))
    ec7 = dict(sorted(ec7.items(), key=lambda x: x[1], reverse=True))

    # Merge the dictionaries
    tmp_dict = {**ec1, **ec2, **ec3, **ec4, **ec5, **ec6, **ec7}
    return tmp_dict 


def calculate_ec_frequency(yeast, gene2ec_dict):
    ec_frequency_dict = {}
    ec_list = []
    for reaction in yeast.reactions:
        tmp = str(reaction.gpr).replace(')','').replace('(','').replace(' and ',' ').replace(' or ',' ').split(' ')
        ec_list += gene_list2_ec_list(tmp, gene2ec_dict)
    ec_frequency_dict = list_to_frequency_dict(ec_list)
    ec_frequency_dict = rank_ec_number_dict(ec_frequency_dict)
    return ec_frequency_dict



def create_frequency_df(yeast8U_ec_frequency_dict,yeast8_ec_frequency_dict):
    yeast8U_yeast8_ec_frequency = {'EC':[],
                                   'yeast8':[],
                                   'yeast8U':[]}
    for k,v in yeast8U_ec_frequency_dict.items():
        yeast8U_yeast8_ec_frequency['EC'].append(k)
        yeast8U_yeast8_ec_frequency['yeast8U'].append(v)
        if k in yeast8_ec_frequency_dict.keys():
            yeast8U_yeast8_ec_frequency['yeast8'].append(yeast8_ec_frequency_dict[k])
        else:
            yeast8U_yeast8_ec_frequency['yeast8'].append(0)

    yeast8U_yeast8_ec_frequency_df = pd.DataFrame(yeast8U_yeast8_ec_frequency)
    yeast8U_yeast8_ec_frequency_df['ratio'] = yeast8U_yeast8_ec_frequency_df.apply(
    lambda row: row['yeast8U'] / row['yeast8'] if row['yeast8'] != 0 else 0,
    axis=1
)
    return yeast8U_yeast8_ec_frequency_df

def get_exchange_reaction(Target_met,model):
    if Target_met!='':
        exchange_reaction = ''
        for i in model.reactions:
            if Target_met in i.reaction and len(i.metabolites) == 1:
                exchange_reaction = i.id
        if exchange_reaction == '':
            reaction_name = 'DM_'+ Target_met
            reaction = Reaction(reaction_name)
            reaction.name = reaction_name
            reaction.subsystem = ''
            reaction.lower_bound = 0.  # This is the default
            reaction.upper_bound = 1000.  # This is the default

            reactant_met_num = {Target_met:-1}
            for met_id, coeff in reactant_met_num.items():
                met = model.metabolites.get_by_id(met_id)
                if met:
                    reaction.add_metabolites({met: coeff})
            model.add_reactions([reaction])  
            exchange_reaction = reaction_name     
        return(exchange_reaction)
    else:
        pass

def calculate_carbon_count(smiles):
    mol = Chem.MolFromSmiles(smiles)
    atom_count1 = Counter([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    c_count1 = atom_count1[6]  # 6 represents carbon
    return c_count1

def process_rule(rule_smiles, mets_total_smiles,cutoff):
    score = []
    smiles = []
    if 'C' in rule_smiles or 'c' in rule_smiles:
        for j in mets_total_smiles:
            current_score = calculate_similarity(rule_smiles, j)
            if current_score > cutoff:
                score.append(current_score)
                smiles.append(j)
    else:
        smiles.append(rule_smiles)
        score.append(1)

    if len(score) == 0:
        return rule_smiles, [], []
    # elif len(score) > 20:
    #     top20_data = sorted(zip(score, smiles), reverse=True)[:20]
    else:
        data = sorted(zip(score, smiles), reverse=True)
    
    scores, smiles = zip(*data)
    
    return rule_smiles, list(smiles), list(scores)

def process_retrorules(retrorules,s):
    # Load retrorules from a CSV file

    # Initialize an empty list to store product smiles
    rules_product_smiles_lst = []

    for index, row in tqdm(retrorules.iterrows()):
        tmp = row[s].split('.')
        for i in tmp:
            if i not in rules_product_smiles_lst:
                rules_product_smiles_lst.append(i)

    # Normalize the smiles and remove duplicates
    if s == 'product_smiles':
        rules_product_smiles_lst = list(set([normalize_smiles(x) for x in rules_product_smiles_lst]))

    print(len(rules_product_smiles_lst))

    return rules_product_smiles_lst


def process_rules_with_multiprocessing(rules_product_smiles_lst, yeast_total_smiles, num_processes=60,cutoff = 0.3  ):
    pool = mp.Pool(num_processes)

    # Initialize a dictionary to store the results
    retrosys_smiles_calculate_similarity = {'rules_smiles': [],
                                            'smiles_in_mets_total_smiles': [],
                                            'scores': []}

    # Create a partial function with the fixed argument mets_total_smiles
    process_rule_partial = partial(process_rule, mets_total_smiles=yeast_total_smiles,cutoff = cutoff)

    # Use the multiprocessing pool to apply process_rule_partial to each element in rules_product_smiles_lst
    for result in tqdm(pool.imap(process_rule_partial, rules_product_smiles_lst), total=len(rules_product_smiles_lst)):
        if result[0] not in retrosys_smiles_calculate_similarity['rules_smiles']:
            retrosys_smiles_calculate_similarity['rules_smiles'].append(result[0]) 
            retrosys_smiles_calculate_similarity['smiles_in_mets_total_smiles'].append(result[1]) 
            retrosys_smiles_calculate_similarity['scores'].append(result[2]) 

    pool.close()
    pool.join()

    return retrosys_smiles_calculate_similarity



def merge_smiles_similarity_rule(lipid_retrorules, lipid_retrosys_smiles_calculate_similarity_pd):
    lipid_retrorules['smiles_similarity_total'] = None
    for index, row in tqdm(lipid_retrorules.iterrows(),total=len(lipid_retrorules)):
        rxn_id, _, reactants, products = row['MNX_ID'], row['classifs'], row['substrate_smiles'], row['product_smiles']
        if rxn_id not in ['MNXR101882', 'MNXR101884', 'MNXR101879', 'MNXR101885', 'MNXR101886', 'MNXR101887',
                         'MNXR103528', 'MNXR103575']:  
            productSMARTs = row['ProductSMARTs'].split('.')
            products = products.split('.')
            new_products_smiles, scores = [], []
            similar_smiles_dict = {}
            for x in range(len(products)):  # prepare mets with substructure and sort the reactants combination
                new_products_smile_one,score_one = get_calculate_similarity(normalize_smiles(products[x]),lipid_retrosys_smiles_calculate_similarity_pd)
                similar_smiles_dict[normalize_smiles(products[x])] = {'smiles_in_mets_total_smiles': new_products_smile_one, 'scores': score_one}
            lipid_retrorules.at[index, 'smiles_similarity_total'] = similar_smiles_dict
    return lipid_retrorules

# def filter_smiles_parallel(new_products_smiles, productSMARTs, n_jobs=4):
#     productSMARTs_mol = [Chem.MolFromSmarts(x) for x in productSMARTs]
#     with mp.Pool(n_jobs) as p:
#         new_products_smiles = list(p.map(process_smi, [(smiles, productSMARTs_mol) for smiles in new_products_smiles]))
#     # new_products_smiles = list(map(process_smi, [(smiles, productSMARTs_mol) for smiles in new_products_smiles]))


#     return new_products_smiles
def process_smi(args):
    try:
        smiles, productSMARTs_mol = args
        mol = Chem.MolFromSmiles(smiles)
        result = [mol.HasSubstructMatch(x) for x in productSMARTs_mol]
        if sum(result) == 0:
            return None
        else:
            return smiles
    except:
        return None

def filter_smiles(new_products_smiles, productSMARTs):
    productSMARTs_mol = [Chem.MolFromSmarts(x) for x in productSMARTs]
    new_products_smiles = list(map(process_smi, [(smiles, productSMARTs_mol) for smiles in new_products_smiles]))

    return new_products_smiles
def process_retrorules_and_save(rxndb_path,failedrxn_path,retrorules,heterologous_met_smiles=None,num_process=50):
    failedrxn = []
    process_retrorule_new_partial = partial(process_retrorule_new,heterologous_met_smiles = heterologous_met_smiles)
    with multiprocessing.Pool(num_process) as pool:

        for result in tqdm(pool.imap_unordered(process_retrorule_new_partial, retrorules.iterrows(),chunksize=1), total=len(retrorules)):
        # for result in tqdm(map(process_retrorule_new_partial, retrorules.iterrows()), total=len(retrorules)):
            rxndb_list = []
            if result[0]:  # 
                rxndb_list.extend(result[0]) 
            
            if result[1] is not None:
                failedrxn.extend(result[1]) 
            index = result[2] 
            rxndb = {}
            for i, rxn_dict in enumerate(rxndb_list):
                key = f'rxn{i+1}'
                rxndb[key] = rxn_dict


            # if len(rxndb) > 0:
            with open(rxndb_path + str(index) + '.json', 'w') as json_file:
                json.dump(rxndb, json_file, indent=4)
    dump_file(failedrxn,failedrxn_path)
    # print('failed_rxn',len(failedrxn))
# def process_retrorule_new(index_row,heterologous_met_smiles = None):
#     rxndb_list_tmp = []
#     # newdbSmiles_tmp = []
#     failedrxn_tmp = []
#     index, row = index_row
#     rxn_id, ECnumber, reactants, products = row['MNX_ID'], row['classifs'], row['substrate_smiles'], row['product_smiles']
#     deprecated_equ_smiles, reactantSMARTs, productSMARTs, rule = row['deprecated_equ_smiles'], row['ReactantsSMARTs'], row['ProductSMARTs'], row['RetroRules']
#     deprecated_equ_smiles = deprecated_equ_smiles.split('>>')
#     if rxn_id not in ['MNXR101882', 'MNXR101884', 'MNXR101879', 'MNXR101885', 'MNXR101886', 'MNXR101887',
#                         'MNXR103528', 'MNXR103575']:  # ATP NADH conversion stuff
#         productSMARTs = productSMARTs.split('.')
#         products = products.split('.')
#         new_products_smiles = []
#         scores = []
#         for x in range(len(products)):  # prepare mets with substructure and sort the reactants combination
#             new_products_smile_one,score_one = get_calculate_similarity_new(row['smiles_similarity_total'][normalize_smiles(products[x])])
#             new_products_smiles.append(new_products_smile_one)
#             scores.append(score_one)
#         if  heterologous_met_smiles:
#             product_prepared = []
#             scores_prepared = []
#             for lst, score_lst in zip(new_products_smiles, scores):
#                 if heterologous_met_smiles in lst:
#                     if len(new_products_smiles) == 1:  # If there is only one sublist
#                         product_prepared.append([heterologous_met_smiles])
#                         scores_prepared.append(score_lst)
#                     else:
#                         other_lists = [sub_lst for sub_lst in new_products_smiles if sub_lst != lst]
#                         tmp_product = list(itertools.product([heterologous_met_smiles], *other_lists))
#                         product_prepared += tmp_product
#                         other_lists_scores = [sub_lst for sub_lst in scores if sub_lst != score_lst]
#                         heterologous_met_index = lst.index(heterologous_met_smiles)
#                         tmp_scores = list(itertools.product([score_lst[heterologous_met_index]], *other_lists_scores))
#                         # tmp_scores = list(itertools.product([score_lst[new_products_smiles.index(lst)]], *scores))
#                         scores_prepared += tmp_scores
#         else:
#             estimated_size = 1
#             for lst in new_products_smiles:
#                 estimated_size *= len(lst)
                
#             if estimated_size > 2000000:
#                 new_products_smiles_ = []
#                 scores_ = []
#                 for lst, score_lst in zip(new_products_smiles, scores):
                    
#                     new_products_smiles_.append(lst[:40])
#                     scores_.append(score_lst[:40])
#                 product_prepared = list(itertools.product(*new_products_smiles_))
#                 scores_prepared = list(itertools.product(*scores_))
#             else:
#                 product_prepared = list(itertools.product(*new_products_smiles))
#                 scores_prepared = list(itertools.product(*scores))

#         for m in range(0, len(product_prepared)):  # apply the rule and output the rxn
#             smile = product_prepared[m]
#             score = scores_prepared[m]
            
#             if len(smile) > 1:
#                 product_smile = '.'.join(smile)
#             else:
#                 product_smile = smile[0]

#             retrorule_tmp = rule.split('>>')
#             retrorule_tmp = '(' + retrorule_tmp[0] + ')>>' + '(' + retrorule_tmp[1] + ')'            
#             product_smile_tmp = product_smile
#             try:
#                 # print('product_smile_tmp',product_smile_tmp)
#                 product_smile_tmp = product_smile_tmp.split('.')
#                 product_smile_tmp = '.'.join([neutralize_charge(x) for x in product_smile_tmp])
#                 product_mol = Chem.MolFromSmiles(product_smile_tmp)
#                 reaction = AllChem.ReactionFromSmarts(retrorule_tmp)
#                 reactant_smiles = reaction.RunReactants((product_mol,))
#                 smiles = [Chem.MolToSmiles(x[0]) for x in reactant_smiles]
#                 reactant_smiles = list(set(smiles))
#                 for reactant_smile in reactant_smiles:
#                     if len(deprecated_equ_smiles) > 0:
#                         if deprecated_equ_smiles[0] == '':
#                             reactant_smile_final = reactant_smile
#                         else:
#                             reactant_smile_final = reactant_smile + '.' + deprecated_equ_smiles[0]
#                         if deprecated_equ_smiles[1] == '':
#                             product_smile_final = product_smile_tmp
#                         else:
#                             product_smile_final = product_smile_tmp + '.' + deprecated_equ_smiles[1]
#                     rxn_smiles = reactant_smile + '>>' + product_smile_tmp
#                     rxn_smiles_final = reactant_smile_final + '>>' + product_smile_final
                    
#                     tmp = {
#                         'EC number': ECnumber,
#                         'rule': rule,
#                         'templateID': rxn_id,
#                         'templateSubstrate': products,###
#                         'rxn_smiles_basic': rxn_smiles,
#                         'rxn_smiles_final': rxn_smiles_final,
#                         'reactant_smile': reactant_smile,
#                         'productsmile': product_smile_tmp,
#                         'similarity': score
#                     }
#                     if check_rule(tmp['rxn_smiles_final'])=='blance':
#                         rxndb_list_tmp.append(tmp)  
#                         # mets_tmp = reactant_smile.split('.')      
#                         # for reactant in mets_tmp:             
#                         #     if reactant not in mets_total_smiles and reactant not in newdbSmiles_tmp:         
#                         #         newdbSmiles_tmp.append(reactant)        
#             except:
#                 failedrxn_tmp.append(rxn_id)
#     return rxndb_list_tmp, failedrxn_tmp, index
def process_retrorule_new(index_row,heterologous_met_smiles = None):
    rxndb_list_tmp = []
    # newdbSmiles_tmp = []
    failedrxn_tmp = []
    index, row = index_row
    rxn_id, ECnumber, reactants, products = row['MNX_ID'], row['classifs'], row['substrate_smiles'], row['product_smiles']
    deprecated_equ_smiles, reactantSMARTs, productSMARTs, rule = row['deprecated_equ_smiles'], row['ReactantsSMARTs'], row['ProductSMARTs'], row['RetroRules']
    deprecated_equ_smiles = deprecated_equ_smiles.split('>>')
    if rxn_id not in ['MNXR101882', 'MNXR101884', 'MNXR101879', 'MNXR101885', 'MNXR101886', 'MNXR101887',
                        'MNXR103528', 'MNXR103575']:  # ATP NADH conversion stuff
        productSMARTs = productSMARTs.split('.')
        products = products.split('.')
        new_products_smiles = []
        scores = []
        for x in range(len(products)):  # prepare mets with substructure and sort the reactants combination
            new_products_smile_one,score_one = get_calculate_similarity_new(row['smiles_similarity_total'][normalize_smiles(products[x])])
            new_products_smiles.append(new_products_smile_one)
            scores.append(score_one)
        if  heterologous_met_smiles:
            product_prepared = []
            scores_prepared = []
            for lst, score_lst in zip(new_products_smiles, scores):
                if heterologous_met_smiles in lst:
                    if len(new_products_smiles) == 1:  # If there is only one sublist
                        product_prepared.append([heterologous_met_smiles])
                        scores_prepared.append(score_lst)
                    else:
                        other_lists = [sub_lst for sub_lst in new_products_smiles if sub_lst != lst]
                        tmp_product = list(itertools.product([heterologous_met_smiles], *other_lists))
                        product_prepared += tmp_product
                        other_lists_scores = [sub_lst for sub_lst in scores if sub_lst != score_lst]
                        heterologous_met_index = lst.index(heterologous_met_smiles)
                        tmp_scores = list(itertools.product([score_lst[heterologous_met_index]], *other_lists_scores))
                        # tmp_scores = list(itertools.product([score_lst[new_products_smiles.index(lst)]], *scores))
                        scores_prepared += tmp_scores
        else:
            estimated_size = 1
            for lst in new_products_smiles:
                estimated_size *= len(lst)
                
            if estimated_size > 2000000:
                new_products_smiles_ = []
                scores_ = []
                for lst, score_lst in zip(new_products_smiles, scores):
                    
                    new_products_smiles_.append(lst[:40])
                    scores_.append(score_lst[:40])
                product_prepared = list(itertools.product(*new_products_smiles_))
                scores_prepared = list(itertools.product(*scores_))
            else:
                product_prepared = list(itertools.product(*new_products_smiles))
                scores_prepared = list(itertools.product(*scores))

        for m in range(0, len(product_prepared)):  # apply the rule and output the rxn
            smile = product_prepared[m]
            score = scores_prepared[m]
            
            if len(smile) > 1:
                product_smile = '.'.join(smile)
            else:
                product_smile = smile[0]

            retrorule_tmp = rule.split('>>')
            retrorule_tmp = '(' + retrorule_tmp[0] + ')>>' + '(' + retrorule_tmp[1] + ')'            
            product_smile_tmp = product_smile
            try:
                # print('product_smile_tmp',product_smile_tmp)
                product_mol = Chem.MolFromSmiles(product_smile_tmp)
                reaction = AllChem.ReactionFromSmarts(retrorule_tmp)
                reactant_smiles = reaction.RunReactants((product_mol,))
                smiles = [Chem.MolToSmiles(x[0]) for x in reactant_smiles]
                reactant_smiles = list(set(smiles))
                for reactant_smile in reactant_smiles:
                    if len(deprecated_equ_smiles) > 0:
                        if deprecated_equ_smiles[0] == '':
                            reactant_smile_final = reactant_smile
                        else:
                            reactant_smile_final = reactant_smile + '.' + deprecated_equ_smiles[0]
                        if deprecated_equ_smiles[1] == '':
                            product_smile_final = product_smile_tmp
                        else:
                            product_smile_final = product_smile_tmp + '.' + deprecated_equ_smiles[1]
                    rxn_smiles = reactant_smile + '>>' + product_smile_tmp
                    rxn_smiles_final = reactant_smile_final + '>>' + product_smile_final
                    # print('rxn_smiles_final',rxn_smiles_final)
                    tmp = {
                        'EC number': ECnumber,
                        'rule': rule,
                        'templateID': rxn_id,
                        'templateSubstrate': products,###
                        'rxn_smiles_basic': rxn_smiles,
                        'rxn_smiles_final': rxn_smiles_final,
                        'reactant_smile': reactant_smile,
                        'productsmile': product_smile_tmp,
                        'similarity': score
                    }
                    if check_rule(tmp['rxn_smiles_final'])=='blance':
                        rxndb_list_tmp.append(tmp)  
                        # mets_tmp = reactant_smile.split('.')      
                        # for reactant in mets_tmp:             
                        #     if reactant not in mets_total_smiles and reactant not in newdbSmiles_tmp:         
                        #         newdbSmiles_tmp.append(reactant)        
            except:
                failedrxn_tmp.append(rxn_id)
    return rxndb_list_tmp, failedrxn_tmp, index
def get_calculate_similarity_new(smiles_similarity_total_dict):
    new_products_smiles = smiles_similarity_total_dict['smiles_in_mets_total_smiles']
    scores = smiles_similarity_total_dict['scores']
    return new_products_smiles, scores
def filter_smiles_parallel(df):
    # newdbSmiles_tmp = []
    for index, row in df.iterrows():
        rxn_id, ECnumber, reactants, products = row['MNX_ID'], row['classifs'], row['substrate_smiles'], row['product_smiles']
        deprecated_equ_smiles, reactantSMARTs, productSMARTs, rule = row['deprecated_equ_smiles'], row['ReactantsSMARTs'], row['ProductSMARTs'], row['RetroRules']
        if rxn_id not in ['MNXR101882', 'MNXR101884', 'MNXR101879', 'MNXR101885', 'MNXR101886', 'MNXR101887',
                            'MNXR103528', 'MNXR103575']:  # ATP NADH conversion stuff
            productSMARTs = productSMARTs.split('.')
            products = products.split('.')
            similar_smiles_dict = {}
            for x in range(len(products)):  # prepare mets with substructure and sort the reactants combination
                new_products_smile_one,score_one = get_calculate_similarity_new(row['smiles_similarity_total'][normalize_smiles(products[x])])
                new_products_smile_one = filter_smiles(new_products_smile_one, productSMARTs)
                combined = list(zip(new_products_smile_one, score_one))
                combined = [(a, b) for a, b in combined if a is not None]
                if combined:
                    new_products_smile_one, score_one = zip(*combined)
                else:
                    new_products_smile_one, score_one = [], []
                similar_smiles_dict[normalize_smiles(products[x])] = {'smiles_in_mets_total_smiles': new_products_smile_one, 'scores': score_one}
            df.at[index, 'smiles_similarity_total'] = similar_smiles_dict
    return df
def filter_smiles_muti(retrorules,retrosys_smiles_calculate_similarity_filter_file_path,num_process=60):
    chunks = np.array_split(retrorules, num_process) 
    with mp.Pool(num_process) as pool:
        results = list(tqdm(pool.imap(filter_smiles_parallel, chunks), total=len(chunks)))
        new_retrorules = pd.concat(results) 
        retrorules_dict = new_retrorules.to_dict(orient='records')
        with open(retrosys_smiles_calculate_similarity_filter_file_path,'w') as f:
            json.dump(retrorules_dict,f)
    # return retrorules_dict
def smiles_has_carbon(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        has_carbon = any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
        return has_carbon
    except:
        return None
def smiles_has_phosphorus(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    has_phosphorus = any(atom.GetSymbol() == 'P' for atom in mol.GetAtoms())
    return has_phosphorus

def process_retrorules(retrorules,s):
    # Load retrorules from a CSV file

    # Initialize an empty list to store product smiles
    rules_product_smiles_lst = []

    for index, row in tqdm(retrorules.iterrows()):
        tmp = row[s].split('.')
        for i in tmp:
            if i not in rules_product_smiles_lst:
                rules_product_smiles_lst.append(i)

    # Normalize the smiles and remove duplicates
    if s == 'product_smiles':
        rules_product_smiles_lst = list(set([normalize_smiles(x) for x in rules_product_smiles_lst]))

    print(len(rules_product_smiles_lst))

    return rules_product_smiles_lst
def process_file(file,rxndb_path,inchikey0,rxndb_drop_path):
    if not os.path.exists(rxndb_drop_path + file):

        try:
            with open(rxndb_path + file, 'r') as f:
                rxndb = json.load(f)    
            remove = []
            for k,v in rxndb.items():
                try:
                    reaction = v['rxn_smiles_basic']
                    reactant_smiles = reaction.split('>>')[0].split('.')
                    for reactants in reactant_smiles:
                        assert smiles2inchikey0(reactants) in inchikey0
                    
                except:
                    #remove the reaction that contains the unknown compound
                    remove.append(k)
            for k in remove:
                rxndb.pop(k)
            with open(rxndb_drop_path + file, 'w') as f:
                json.dump(rxndb, f,indent=4)
        except:
            pass

def drop_rxndb(rxndb_path, rxndb_drop_path, total_met,num_processes=60):
    inchikey0 = load_pickle(total_met)
    files = os.listdir(rxndb_path)
    process_file_partial = partial(process_file,rxndb_path=rxndb_path,inchikey0=inchikey0,rxndb_drop_path=rxndb_drop_path)
    with mp.Pool(num_processes) as p:
        list(tqdm(p.imap_unordered(process_file_partial, files), total=len(files)))
def load_and_count(file_path):
    with open(file_path, 'r') as f:
        rxndb = json.load(f)
    return len(rxndb)

def statistic_reaction_num(rxndb_path, rxndb_drop_path):
    files = [os.path.join(rxndb_path, file) for file in os.listdir(rxndb_path)]
    files_drop = [os.path.join(rxndb_drop_path, file) for file in os.listdir(rxndb_drop_path)]

    with ProcessPoolExecutor(max_workers=50) as executor:
        num = sum(tqdm(executor.map(load_and_count, files), total=len(files)))
        num_drop = sum(tqdm(executor.map(load_and_count, files_drop), total=len(files_drop)))

    print('rxndb:', num)
    print('rxndb_drop:', num_drop)

def merge_rxndb(rxndb_drop_path,rxndb_all_path,name='not_lipid'):
    rxndb_all = []
    files = os.listdir(rxndb_drop_path)

    for file in files:
        with open(rxndb_drop_path + file, 'r') as f:
            rxndb = json.load(f)
            rxndb_all.extend(list(rxndb.values()))
    rxndb = {}
    for i, rxn_dict in enumerate(rxndb_all):
        if name == 'not_lipid':
            key = f'rxn{i+1}'
            rxndb[key] = rxn_dict
        else:
            key = f'rxnl{i+1}'
            rxndb[key] = rxn_dict
    with open(rxndb_all_path, 'w') as f:
        json.dump(rxndb, f,indent=4)
    return rxndb

def process_smiles(i,smile_max_score,yeast_total_smiles):    ###The highest similarity scores of all the reactants and model metabolites of rxn were obtained.
    if i not in smile_max_score['smile']:
        result = {'smile': i, 'score': 0, 'sim_smile': ''}
        score_init = 0
        sim_smile = ' '
        for j in yeast_total_smiles:
            tmp_score = compare_smiles(i, j)  
            if compare_smiles_inchikey(i,j) or (tmp_score == 1 and are_atom_counts_equal(i,j)):##for the compounds without carbon atom
                score_init = 1
                sim_smile = j
        result['score'] = score_init
        result['sim_smile'] = sim_smile
        return result
    else: return {'smile': '', 'score': '', 'sim_smile': ''}

def calculate_smile_max_score(rxndb_all_reactant_smile,yeast_total_smiles,num_processes=60):
    smile_max_score = {'smile': [], 'score': [], 'sim_smile': []}

    pool = multiprocessing.Pool(num_processes)

    partial_process_smiles = partial(process_smiles, smile_max_score=smile_max_score,yeast_total_smiles=yeast_total_smiles)
    results = list(tqdm(pool.imap(partial_process_smiles, rxndb_all_reactant_smile), total=len(rxndb_all_reactant_smile)))
    # results = list(tqdm(pool.imap(process_smiles, rxndb_all_reactant_smile), total=len(rxndb_all_reactant_smile)))

    pool.close()
    pool.join()

    for result in results:
        smile_max_score['smile'].append(result['smile'])
        smile_max_score['score'].append(result['score'])
        smile_max_score['sim_smile'].append(result['sim_smile'])
        
    return(smile_max_score)

def process_chunk(chunk_df):
    chunk_all_smile = []
    for i in chunk_df['rxn_smiles_basic']:
        tmp_list = i.split('>>')[0].split('.')
        for j in tmp_list:
            if j not in chunk_all_smile:
                chunk_all_smile.append(j)
        tmp_list_product = i.split('>>')[1].split('.')
        for j in tmp_list_product:
            if j not in chunk_all_smile:
                chunk_all_smile.append(j)
    return chunk_all_smile
def parallel_process_rxn_smiles(rxndb_df_test, n_splits=80,num_process=50):
    chunks = np.array_split(rxndb_df_test, n_splits)  # split the DataFrame into n parts

    with multiprocessing.Pool(num_process) as pool:
        results = list(tqdm(pool.imap(process_chunk, chunks), total=len(chunks)))

    # merging all processes
    rxndb_all_smile = [item for sublist in results for item in sublist]
    rxndb_all_smile = list(set(rxndb_all_smile))  # deduplication ( if not done before )
    print('number of metabolites in RXNDB:', len(rxndb_all_smile))
    return rxndb_all_smile

def calculate_save_smiles_max_score(rxndb_all_smiles, yeast8_total_smiles,rxndb_met_max_score_file,num_processes=60):
    smiles_max_score = calculate_smile_max_score(rxndb_all_smiles,yeast8_total_smiles,num_processes=num_processes)
    smiles_max_score_pd = pd.DataFrame(smiles_max_score)
    smiles_max_score_pd.to_csv(rxndb_met_max_score_file,index=None)
    return smiles_max_score_pd
def get_score_from_smiles(input_smiles,smile_max_score_pd): ###return highest similarity score
    # input_smiles = normalize_smiles(input_smiles)
    row = smile_max_score_pd[smile_max_score_pd['smile'] == input_smiles]
    if not row.empty:
        return row['score'].max()
    else:
        return 0
def process_reaction_product(index_row,smile_max_score_pd):
    index, row = index_row
    reactant_smiles = row['rxn_smiles_basic'].split('>>')[0].split('.')
    product_smiles = row['rxn_smiles_basic'].split('>>')[1].split('.')
    scores = []

    for i in product_smiles:
        scores.append(get_score_from_smiles(i, smile_max_score_pd))

    scores_all_0 = all(score == 1 for score in scores)

    if scores_all_0 and len(reactant_smiles)>0:
        return reactant_smiles, row['NO']
    else:
        return [], None
def process_reaction_reactant(index_row,smile_max_score_pd):
    index, row = index_row
    reactant_smiles = row['rxn_smiles_basic'].split('>>')[0].split('.')
    product_smiles = row['rxn_smiles_basic'].split('>>')[1].split('.')
    scores = []

    for i in reactant_smiles:
        scores.append(get_score_from_smiles(i, smile_max_score_pd))

    scores_all_0 = all(score == 1 for score in scores)

    if scores_all_0 and len(product_smiles)>0:
        return product_smiles, row['NO']
    else:
        return [], None
def process_chunk_reactant(chunk_df,smile_max_score_pd):
    smiles_success = []
    success_rxndbid = []
    for result in map(partial(process_reaction_reactant,smile_max_score_pd=smile_max_score_pd), chunk_df.iterrows()):
        if result[0]:  # check whether the result is non-empty
            smiles_success.extend(result[0])
        if result[1] is not None:
            success_rxndbid.append(result[1])
    return smiles_success, success_rxndbid
def process_chunk_product(chunk_df,smile_max_score_pd):
    smiles_success = []
    success_rxndbid = []
    for result in map(partial(process_reaction_product,smile_max_score_pd=smile_max_score_pd), chunk_df.iterrows()):
        if result[0]:  # check whether the result is non-empty
            smiles_success.extend(result[0])
        if result[1] is not None:
            success_rxndbid.append(result[1])
    return smiles_success, success_rxndbid
def process_reactions_in_parallel_reactant(rxndb,origin_smile_max_score_pd, num_processes=100, num_iterations=1):
    num = 0
    tmp_smile_max_score_pd = origin_smile_max_score_pd
    while num < num_iterations:
        num+=1
        # pool = multiprocessing.Pool(num_processes)
        smiles_success = []
        success_rxndbid = []
        chunks = np.array_split(rxndb, 100)
        with multiprocessing.Pool(num_processes) as pool:
            for result in tqdm(pool.imap(partial(process_chunk_reactant,smile_max_score_pd=tmp_smile_max_score_pd), chunks), total=len(chunks)):
                smiles_success.extend(result[0])
                success_rxndbid.extend(result[1])
            for result in tqdm(pool.imap(partial(process_chunk_product,smile_max_score_pd=tmp_smile_max_score_pd), chunks), total=len(chunks)):
                smiles_success.extend(result[0])
                success_rxndbid.extend(result[1])
        # pool.close()
        # pool.join()
        
        smiles_success = list(set(smiles_success))
        success_rxndbid = list(set(success_rxndbid))
        # Process cumulative successful SMILES
        for smile in smiles_success:
            if get_score_from_smiles(smile, tmp_smile_max_score_pd) < 1:
                new_row = {'smile': smile, 'score': 1, 'sim_smile': 'sys'}
                # tmp_smile_max_score_pd = tmp_smile_max_score_pd._append(new_row, ignore_index=True)
                tmp_smile_max_score_pd = tmp_smile_max_score_pd.append(new_row, ignore_index=True)

        # smiles_success = [normalize_smiles(met) for met in smiles_success]
        # smiles_success = list(set(smiles_success))

        print(f'Iteration {num} - Current success count: {len(smiles_success)}')
        print(f'Iteration {num} - Current success_rxndbid count: {len(success_rxndbid)}')
        print('============================================================================')

    print('final success', len(smiles_success))
    print('final success_rxndbid', len(success_rxndbid))
    # print(success_rxndbid)
    return smiles_success, success_rxndbid,tmp_smile_max_score_pd
def save_rxndb_to_model(rxndb_df,success_rxndbid,rxndb_to_model_path):
    rxndb_to_model = rxndb_df[rxndb_df['NO'].isin(success_rxndbid)].reset_index(drop=True)
    print(rxndb_to_model.shape)
    rxndb_to_model = rxndb_to_model.drop_duplicates(subset=['templateID', 'rxn_smiles_final'], keep='first').reset_index(drop=True)
    print(rxndb_to_model.shape)
    rxndb_to_model = rxndb_to_model[['NO','EC number','templateID','rxn_smiles_basic','rxn_smiles_final']]
    print(rxndb_to_model.shape)
    rxndb_to_model.to_csv(rxndb_to_model_path,index=None)
    print(rxndb_to_model.shape)
    # return rxndb_to_model

def save_success_fail_target_smiles(target_smiles_file,smiles_success,YMDB_success_met_smile_file,YMDB_fail_met_smile_file,uptake_smiles=None):
    target_smiles = load_pickle(target_smiles_file)
    # uptake_smiles_inchikey = [smiles2inchikey0(i) for i in uptake_smiles]
    # target_uptake = [i for i in target_smiles if smiles2inchikey0(i) in uptake_smiles_inchikey]
    # print('target_uptake:',len(target_uptake))
    # target_no_uptake = [i for i in target_smiles if i not in target_uptake]
    inchikey_success = [smiles2inchikey0(i) for i in smiles_success]
    success_target_smiles = []
    for i in target_smiles:
        if '.' not in i:
            if smiles2inchikey0(i) in inchikey_success:
                success_target_smiles.append(i)
        else:
            smi_lst = i.split('.')
            smi_lst_with_carbon = [j for j in smi_lst if smiles_has_carbon(j)]
            if all(smiles2inchikey0(j) in inchikey_success for j in smi_lst_with_carbon):
                success_target_smiles.append(i)
            
    print('success number:',len(success_target_smiles))
    dump_file(success_target_smiles,YMDB_success_met_smile_file)
    fail_target_smiles = [i for i in target_smiles if i not in success_target_smiles]
    print('fail number:',len(fail_target_smiles))
    dump_file(fail_target_smiles,YMDB_fail_met_smile_file)
# def save_success_fail_target_smiles(target_smiles_file,smiles_success,YMDB_success_met_smile_file,YMDB_fail_met_smile_file):
#     target_smiles = load_pickle(target_smiles_file)
#     # smiles_success = smiles_success + uptake_smiles
#     inchikey_success = [smiles2inchikey0(i) for i in smiles_success]
#     success_target_smiles = []
#     for i in target_smiles:
#         if smiles2inchikey0(i) in inchikey_success:
#             success_target_smiles.append(i)
            
#     print('success number:',len(success_target_smiles))
#     dump_file(success_target_smiles,YMDB_success_met_smile_file)
#     fail_target_smiles = [i for i in target_smiles if i not in success_target_smiles]
#     print('fail number:',len(fail_target_smiles))
#     dump_file(fail_target_smiles,YMDB_fail_met_smile_file)
def rxndb_gene_annotation(rxndb_to_model_path,ec2gene_dict):
    rxndb = pd.read_csv(rxndb_to_model_path,index_col=None).astype(str)
    rxndb['EC number'] = rxndb['EC number'].apply(lambda x:x.split(';') if x else [])
    rxndb['EC number'] = rxndb['EC number'].apply(lambda x:list(set([".".join(value.split(".")[:3]) for value in x])))
    rxndb = rxndb[rxndb['EC number'].apply(lambda x: len(x) > 0)]  #Eliminate rows with ec empty
    rxndb['GENE'] = rxndb['EC number'].apply(lambda x: [gene for ec in x if ec in ec2gene_dict for gene in ec2gene_dict[ec]] if x else []) 
    ####all EC
    all_EC_list = []
    for i in rxndb['EC number']:
        all_EC_list+=i
    all_EC_list = list(set(all_EC_list))
    print('rxndb all ec num',len(all_EC_list))    
    return rxndb

def get_all_met_smile(rxndb_GPR_to_model_path):
    rxndb_to_model = pd.read_csv(rxndb_GPR_to_model_path,index_col=None)
    # rxndb_to_model = rxndb_to_model.head(10000)
    rxndb_all_smile = []
    for i in tqdm(rxndb_to_model['rxn_smiles_final']):
        smiles_list = i.replace('>>','.').split('.')
        for j in smiles_list:
            if j not in rxndb_all_smile:
                rxndb_all_smile.append(j)
    rxndb_all_smile = list(set(rxndb_all_smile))
    print(len(set(rxndb_all_smile)))
    return rxndb_all_smile

def process_smiles_annotation(i,yeast_total_smiles):    ###The highest similarity scores of all the reactants and model metabolites of rxn were obtained
    result = {'smile': i, 'score': 0, 'sim_smile': ''}
    score_init = 0
    sim_smile = ' '
    for j in yeast_total_smiles:
        tmp_score = compare_smiles(i, j)  
        if compare_smiles_inchikey(i,j) or (tmp_score == 1 and are_atom_counts_equal(i,j)):##后面是针对非碳原子
            score_init = 1
            sim_smile = j
            break
    result['score'] = score_init
    result['sim_smile'] = sim_smile
    return result
        

def calculate_smile_max_score_annotation(rxndb_all_reactant_smile,yeast_total_smiles,num_process=30):
    smile_max_score = {'smile': [], 'score': [], 'sim_smile': []}

    pool = multiprocessing.Pool(num_process)

    partial_process_smiles = partial(process_smiles_annotation,yeast_total_smiles=yeast_total_smiles)
    results = list(tqdm(pool.imap(partial_process_smiles, rxndb_all_reactant_smile), total=len(rxndb_all_reactant_smile)))
    # results = list(tqdm(pool.imap(process_smiles, rxndb_all_reactant_smile), total=len(rxndb_all_reactant_smile)))

    pool.close()
    pool.join()

    for result in results:
        smile_max_score['smile'].append(result['smile'])
        smile_max_score['score'].append(result['score'])
        smile_max_score['sim_smile'].append(result['sim_smile'])
        
    return(smile_max_score)
def get_yeast8_id_smiles_mapping(model_path):
    model = pd.read_csv(model_path)
    yeast_id_mapping = model[['REPLACEMENT ID','COMPARTMENT','standard_smiles']]
    yeast_id_mapping.columns = ['met_id','compartment','smiles']
    yeast_id_mapping['smiles'] = yeast_id_mapping['smiles'].apply(normalize_smiles)
    yeast_id_mapping['inchikey0'] = yeast_id_mapping['smiles'].apply(smiles2inchikey0)
    return yeast_id_mapping


def get_gem_met_id_for_smile(smiles,yeast_id_smiles_mapping):
    ID ='error'
    compartment= 'c'
    tmp_df = yeast_id_smiles_mapping[(yeast_id_smiles_mapping['smiles']==smiles)]
    tmp_df_compartment_c = tmp_df[tmp_df['compartment']=='c']
    if len(tmp_df_compartment_c)>0:
        ID = tmp_df_compartment_c['met_id'].to_list()[0]
    else:
        ID = tmp_df['met_id'].to_list()[0]
        compartment = tmp_df['compartment'].to_list()[0]
    if ID == 'error':
        print('error')
        return(ID)
    else:
        return(ID,compartment)
    
def assign_ID_to_RXNDB_met(smile_max_score_pd,yeast_id_smiles_mapping,s = 'lipid'):

    sim_smile_met_id = []
    smiles_id_dict = {}
    num = 0
    for index,row in tqdm(smile_max_score_pd.iterrows(),total = smile_max_score_pd.shape[0]):
        if row['score'] == 1:
            # sim_smile_met_id.append(get_gem_met_id_for_smile(row['sim_smile'],yeast8_MNX_CHBEI_id_mapping))
            sim_smile_met_id.append(get_gem_met_id_for_smile(row['sim_smile'],yeast_id_smiles_mapping))
        else:
            inchikey0 = smiles2inchikey0(row['smile'])
            if inchikey0 in smiles_id_dict:
                id = smiles_id_dict[inchikey0]
            else:
                num +=1
                if s == 'lipid':
                    id = 'sl_' +str(num)
                else:
                    id = 'sn_' +str(num)
                smiles_id_dict[inchikey0] = id
            sim_smile_met_id.append((id,'c'))
    smile_max_score_pd['sim_smile_met_id'] = sim_smile_met_id

    return smile_max_score_pd    


def get_new_met_smile_list(smile_max_score_pd_with_ID):
    smile_max_score_pd_with_ID['ID'] = smile_max_score_pd_with_ID['sim_smile_met_id'].apply(lambda x:x[0])
    smile_max_score_pd_with_ID['compartment'] = smile_max_score_pd_with_ID['sim_smile_met_id'].apply(lambda x:x[1])
    smile_max_score_pd_with_ID.rename(columns={'smile':'new_met_smiles'},inplace=True)
    smile_max_score_pd_with_ID = smile_max_score_pd_with_ID[['new_met_smiles','sim_smile','ID','compartment']]
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[NH4+]','ID'] = 's_0419'
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[H+]','ID'] = 's_0794'
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[H]O[H]','ID'] = 's_0803'
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='O','ID'] = 's_0803'
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[OH]','ID'] = 's_0803'
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[Fe+3]','ID'] = 's_3855'
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[Fe+2]','ID'] = 's_0924'
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[Cl-]','ID'] = 's_3778'
    # smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[Cu+]','ID'] = 's_4019'
    smile_max_score_pd_with_ID.loc[smile_max_score_pd_with_ID['new_met_smiles']=='[Cu+2]','ID'] = 's_4019'

    return smile_max_score_pd_with_ID

def get_smiles_to_id_mapping(metabolites_info_to_GEM_path):
    met_info = pd.read_csv(metabolites_info_to_GEM_path)
    smiles_to_id_mapping = {}
    for index, row in met_info.iterrows():
        new_met_smiles = row['new_met_smiles']
        key = row['ID']
        smiles_to_id_mapping[new_met_smiles] = key
    return smiles_to_id_mapping

def parse_reaction_smiles(reaction_smiles):
    reactant_smiles, product_smiles = reaction_smiles.split(">>")
    reactants = reactant_smiles.split(".")
    products = product_smiles.split(".")
    return reactants, products

def count_compounds(compounds, coefficients):
    compound_counts = defaultdict(int)
    for compound, coefficient in zip(compounds, coefficients):
        compound_counts[compound] += coefficient
    return dict(compound_counts)

def map_smiles_to_id(smiles, mapping):
    return mapping[smiles]  
    # return mapping.get(smiles, smiles)  # return the original SMILES if the map does not exist

def process_reaction_smiles(reaction_smiles, smiles_to_id_mapping):
    reactants, products = parse_reaction_smiles(reaction_smiles)

    reactant_mols = [Chem.MolFromSmiles(reactant) for reactant in reactants]
    product_mols = [Chem.MolFromSmiles(product) for product in products]

    reactant_coefficients = [-1] * len(reactant_mols)
    product_coefficients = [1] * len(product_mols)

    reactant_counts = count_compounds(reactants, reactant_coefficients)
    product_counts = count_compounds(products, product_coefficients)

    final_counts = {}
    for compound, count in reactant_counts.items():
        mapped_compound = map_smiles_to_id(compound, smiles_to_id_mapping)
        final_counts[mapped_compound] = final_counts.get(mapped_compound, 0) + count

    for compound, count in product_counts.items():
        mapped_compound = map_smiles_to_id(compound, smiles_to_id_mapping)
        final_counts[mapped_compound] = final_counts.get(mapped_compound, 0) + count

    # remove the term with a coefficient of 0
    final_counts = {compound: count for compound, count in final_counts.items() if count != 0}

    return final_counts


def convert_reaction_format(input_string):
    input_string = str(input_string)
    # Remove curly braces and single quotes to get a clean dictionary string
    dict_string = re.sub(r"[{}']", "", input_string)

    # Split the dictionary string into individual key-value pairs
    pairs = dict_string.split(", ")

    # Initialize dictionaries to store coefficients for reactants and products
    reactants = {}
    products = {}

    # Parse each key-value pair and populate the reactants and products dictionaries
    for pair in pairs:
        key, value = pair.split(": ")
        if int(value) < 0:
            reactants[key] = -int(value)
        elif int(value) > 0:
            products[key] = int(value)

    # Create the reaction string in the desired format
    reaction_string = " + ".join([f"{coeff if coeff != 1 else ''} {species}" for species, coeff in reactants.items()])
    reaction_string += " => "
    reaction_string += " + ".join([f"{coeff if coeff != 1 else ''} {species}" for species, coeff in products.items()])

    return reaction_string

def get_rxndb_to_model_total_info(rxndb_GPR_to_model_path,smiles_to_id_mapping):
    rxndb_to_model = pd.read_csv(rxndb_GPR_to_model_path,index_col=None)
    # rxndb_to_model = rxndb_to_model.head(10000)
    rxndb_to_model['GPR'] = rxndb_to_model['GENE'].apply(lambda x:x.replace('[]','nogene').replace("['",'(').replace("']",')').replace("', '",') or ('))
    rxndb_to_model['equation_dict'] = rxndb_to_model['rxn_smiles_final'].apply(lambda x:process_reaction_smiles(x, smiles_to_id_mapping))

    rxndb_to_model['equation'] = rxndb_to_model['equation_dict'].apply(lambda x:convert_reaction_format(x) if x else ' =>')
    rxndb_to_model['EC_number'] = rxndb_to_model['EC number'].apply(lambda x:x.replace('[]','noec').replace("', '",', ').replace("['",'').replace("']",''))

    rxndb_to_model = rxndb_to_model[~rxndb_to_model['equation_dict'].apply(lambda x:len(x)<2)]
    rxndb_to_model = rxndb_to_model[~rxndb_to_model['equation'].str.startswith(' =>')]
    return rxndb_to_model


def compare_carbon_atoms_between_formula_smiles(smiles, formula):
    mol = Chem.MolFromSmiles(smiles)
    atom_count1 = Counter([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    c_count1 = atom_count1[6]  # 6 represents carbon
    carbon_count_formula_match = re.search(r'C(\d*)', formula)
    if carbon_count_formula_match:
        carbon_count_formula = int(carbon_count_formula_match.group(1) or 1)
    else:
        carbon_count_formula = 0
    if c_count1 == carbon_count_formula:
        return 1
    else: return 0 


def check_met_and_reaction_formula(yeast870_path,metabolites_info_to_GEM_path,rxndb_to_model_total_info_path):
    ## model
    yeast8 = pd.read_csv(yeast870_path)

    # yeast8 = cobra.io.load_yaml_model(yeast870_path)
    # ## error met
    new_met_info_to_GEM_tmp = pd.read_csv(metabolites_info_to_GEM_path)
    new_met_info_to_GEM_tmp = new_met_info_to_GEM_tmp[new_met_info_to_GEM_tmp['ID'].str.contains('s_')]
    # new_met_info_to_GEM_tmp['formula'] = new_met_info_to_GEM_tmp['ID'].apply(lambda x:yeast8.metabolites.get_by_id(x).formula)
    # new_met_info_to_GEM_tmp = new_met_info_to_GEM_tmp[new_met_info_to_GEM_tmp['formula'].notna()]
    # new_met_info_to_GEM_tmp = new_met_info_to_GEM_tmp[['ID','new_met_smiles','formula']]
    # new_met_info_to_GEM_tmp
    yeast8 = yeast8[yeast8['REPLACEMENT ID'].isin(new_met_info_to_GEM_tmp['ID'])]
    yeast8['res'] = 0 
    for index,row in yeast8.iterrows():
        # new_met_info_to_GEM_tmp.at[index, 'res'] = compare_carbon_atoms_between_formula_smiles(row['new_met_smiles'], row['formula'])
        yeast8.at[index, 'res'] = compare_atoms_between_formula_smiles(row['standard_smiles'], row['COMPOSITION'])
    #new_met_info_to_GEM_tmp

    error_met = yeast8[yeast8['res']==0]['REPLACEMENT ID'].to_list()
    # error_met += ['s_3888']###########################################
    # error_met += ['s_1184','s_0123','s_4000', 's_3888']###########################################
    error_met = set(error_met)
    error_met.discard('s_3778')  #chloride, to avoid mistaken deletion
    print('formula error_met',error_met)

    ## reaction
    rxndb_to_model_total_info = pd.read_csv(rxndb_to_model_total_info_path)
    error_reaction_list = []
    for met in tqdm(error_met,total=len(error_met)):
        for index,row in rxndb_to_model_total_info.iterrows():
            if met in row['equation']:
                if row['NO'] not in error_reaction_list:
                    error_reaction_list.append(row['NO'])
    rxndb_to_model_total_info = rxndb_to_model_total_info[~rxndb_to_model_total_info['NO'].isin(error_reaction_list)]

    ### met
    rxndb_to_model_total_info['equation_dict'] = rxndb_to_model_total_info['equation_dict'].apply(lambda x:literal_eval(x))
    new_met_list = []
    for i in rxndb_to_model_total_info['equation_dict']:
        new_met_list+=list(i.keys())
    new_met_list = list(set(new_met_list))

    new_met_info_to_GEM = pd.read_csv(metabolites_info_to_GEM_path)
    new_met_info_to_GEM = new_met_info_to_GEM[new_met_info_to_GEM['ID'].isin(new_met_list)]

    return new_met_info_to_GEM,rxndb_to_model_total_info


def get_same_rxndb_reaction_mapping(rxndb_total_info_to_model_path):
    rxndb_total_info_to_model = pd.read_csv(rxndb_total_info_to_model_path)
    print(rxndb_total_info_to_model.shape)
    rxndb_total_info_to_model['equation_dict'] = rxndb_total_info_to_model['equation_dict'].apply(lambda x:literal_eval(x))
    ##
    rea_ID2dice = {}
    for index,row in tqdm(rxndb_total_info_to_model.iterrows(),total = len(rxndb_total_info_to_model)):
        rea_ID2dice[row['NO']] = row['equation_dict']    
    ##
    same_rxndb_reaction = {}
    for index,row in tqdm(rxndb_total_info_to_model.iterrows(),total = len(rxndb_total_info_to_model)):
        # tmp_met_dict = sorted(row['equation_dict'])
        tmp_met_dict = row['equation_dict']
        tmp_met_dict = {key: abs(value) for key, value in tmp_met_dict.items()}
        present_ID = ''
        for rea_ID in same_rxndb_reaction.keys():
            # present_reaction_dict = sorted(rxndb_total_info_to_model.loc[rxndb_total_info_to_model['NO']==rea_ID,'equation_dict'])
            present_reaction_dict = rea_ID2dice[rea_ID]
            present_reaction_dict = {key: abs(value) for key, value in present_reaction_dict.items()}
            # if tmp_met_dict == present_reaction_dict:
            if are_dicts_equal(tmp_met_dict,present_reaction_dict):
                present_ID = rea_ID
                break
        if present_ID =='':
            same_rxndb_reaction[row['NO']] = []
        else:
            same_rxndb_reaction[present_ID].append(row['NO'])  
    return  same_rxndb_reaction


def get_no_same_rxndb_reaction_parallelize(rxndb_total_info_to_model_path, func,same_rxndb_reaction, n_cores=4):
    rxndb_total_info_to_model = pd.read_csv(rxndb_total_info_to_model_path)
    rxndb_total_info_to_model_split = np.array_split(rxndb_total_info_to_model, n_cores)
    pool = mp.Pool(n_cores)
    get_no_same_rxndb_reaction_partial = partial(func, same_rxndb_reaction=same_rxndb_reaction,rxndb_total_info_to_model_all=rxndb_total_info_to_model)
    rxndb_total_info_to_model = pd.concat(pool.imap_unordered(get_no_same_rxndb_reaction_partial, rxndb_total_info_to_model_split))
    pool.close()
    pool.join()
    indexes_to_remove = []
    for index, row in tqdm(rxndb_total_info_to_model.iterrows(), total=len(rxndb_total_info_to_model)):
        if row['NO'] not in same_rxndb_reaction.keys():
            indexes_to_remove.append(index)
    rxndb_total_info_to_model.drop(indexes_to_remove, inplace=True)
    return rxndb_total_info_to_model

def get_no_same_rxndb_reaction(rxndb_total_info_to_model,same_rxndb_reaction,rxndb_total_info_to_model_all):
    for index,row in rxndb_total_info_to_model.iterrows():
        if row['NO'] in same_rxndb_reaction.keys():
            GPR = row['GPR'].replace('(','').replace(')','').split(' or ')
            for i in same_rxndb_reaction[row['NO']]:
                try:
                    tmp_GPR = rxndb_total_info_to_model_all.loc[rxndb_total_info_to_model_all['NO']==i,'GPR'].values[0].replace('(','').replace(')','').split(' or ')
                    GPR += [gene for gene in tmp_GPR if gene not in GPR and gene]
                except:
                    continue
            GPR = '(' + ') or ('.join(GPR)+ ')'
            rxndb_total_info_to_model.loc[rxndb_total_info_to_model['NO']==row['NO'],'GPR'] = GPR

    

    return rxndb_total_info_to_model

def process_model(model_path, miss_met_id):
    model_ = pd.read_csv(model_path)
    model_ = model_.dropna(subset=['standard_smiles'])
    model_['normalize_smiles'] = model_['standard_smiles'].apply(normalize_smiles)
    metID_with_smiles = model_['REPLACEMENT ID'].tolist()
    print(len(metID_with_smiles))
    print(metID_with_smiles)
    yeast8_smiles_lst = model_['standard_smiles'].tolist()
    metID_with_smiles = [x for x in metID_with_smiles if x not in miss_met_id]
    yeast8_smiles_lst = list(set(yeast8_smiles_lst))
    print(len(yeast8_smiles_lst))
    print(yeast8_smiles_lst)
    return model_, metID_with_smiles, yeast8_smiles_lst

def get_rxndb_met_with_s_smiles_lst(rxndb_total_info_to_model_path,model_,miss_met_id):
    rxndb_to_model_total_info = pd.read_csv(rxndb_total_info_to_model_path)
    rxndb_to_model_total_info['equation_dict'] = rxndb_to_model_total_info['equation_dict'].apply(lambda x:literal_eval(x))
    rxndb_to_model_total_info['equation_dict'] = rxndb_to_model_total_info['equation_dict'].apply(lambda x:{k:v for k, v in x.items() if k not in miss_met_id})
    # print(rxndb_to_model_total_info.shape)
    # rxndb_to_model_total_info = rxndb_to_model_total_info[rxndb_to_model_total_info['equation_dict'].apply(lambda x : all('s_' in k for k, v in x.items()))]
    print(rxndb_to_model_total_info.shape)
    rxndb_met_with_s_smiles_lst = []
    for i in tqdm(rxndb_to_model_total_info['equation_dict'],total=rxndb_to_model_total_info.shape[0]):
        # print(i)
        for k,v in i.items():
            # print(normalize_smiles(model.metabolites.get_by_id(k).smiles))
            if 's_' in k and model_[model_['REPLACEMENT ID'] == k]['normalize_smiles'].to_list()[0] not in rxndb_met_with_s_smiles_lst:
                rxndb_met_with_s_smiles_lst.append(model_[model_['REPLACEMENT ID'] == k]['normalize_smiles'].to_list()[0])
        # rxndb_met_with_s.append(normalize_smiles(model.metabolites.get_by_id(k).smiles) for k,v in i.items() if 's_' in k and k not in miss_met_id)
    print(len(rxndb_met_with_s_smiles_lst))
    print(rxndb_met_with_s_smiles_lst[:10])
    return rxndb_met_with_s_smiles_lst


def get_reaID_with_smiles_dict(model, model_, miss_met_id,rxndb_met_with_s_smiles_lst,metID_with_smiles):
    reaID_allwith_smiles = []
    for i in model.reactions:
        tmp = {k:v for k,v in i.metabolites.items() if k.id not in miss_met_id}
        if all(k.id in metID_with_smiles for k,v in tmp.items()):
            reaID_allwith_smiles.append(i.id)
    print(len(reaID_allwith_smiles))

    reaID_with_smiles_dict = {}
    for i in reaID_allwith_smiles:
        smiles_list = [model_[model_['REPLACEMENT ID'] == k.id]['normalize_smiles'].to_list()[0] for k, v in model.reactions.get_by_id(i).metabolites.items() if k.id not in miss_met_id]
        # smiles_list = [normalize_smiles(model_gem.metabolites.get_by_id(k.id).smiles) for k, v in model_gem.reactions.get_by_id(i).metabolites.items() if k.id not in miss_met_id]
        if all(x in rxndb_met_with_s_smiles_lst for x in smiles_list) and len(list(set(smiles_list)))>1:
        # if len(list(set(smiles_list)))>1:
            reaID_with_smiles_dict[i] = list(set(smiles_list))
        else:
            pass
            # print(model.reactions.get_by_id(i).reaction)
    print(len(reaID_with_smiles_dict))
    print_first_ten(reaID_with_smiles_dict)
    return reaID_with_smiles_dict


def get_yeast8_reaction_in_rxndb(item,rxndb_to_model_total_info,model_):
    reaID, smiles_list = item
    result = []
    for index, row in rxndb_to_model_total_info.iterrows():
        rxndb_smiles_lst = [model_[model_['REPLACEMENT ID'] == k]['normalize_smiles'].to_list()[0] for k,v in row['equation_dict'].items()]
        if sorted(list(set(rxndb_smiles_lst))) == sorted(list(set(smiles_list))):
            result.append(row['NO'])
    return reaID, result

def GEM_gene_essential(model_path):
    model = cobra.io.load_yaml_model(model_path)
    model.solver = 'gurobi'
    model.optimize()
    initial_biomass = model.reactions.get_by_id('r_2111').flux
    essential_predict_list = []
    notessential_predict_list = []
    gene_lst = [x.id for x in model.genes]
    for i in tqdm(gene_lst,total=len(gene_lst)):
        with model:
            model.genes.get_by_id(i).knock_out()
            model.solver = 'gurobi'
            model.optimize()            
            flux_r_2111 = model.reactions.get_by_id('r_2111').flux
            if flux_r_2111 < 0.1 * initial_biomass:
            # if flux_r_2111 < 0.000001:
                essential_predict_list.append(i)
            else:
                notessential_predict_list.append(i)
    return essential_predict_list,notessential_predict_list
# results is a dictionary where the keys are the reaIDs and the values are the lists of NOs
def get_yeast8_reaction_in_rxndb_parallel(rxndb_total_info_to_model_path,reaID_with_smiles_dict,model_,miss_met_id,num_process):
    yeast8_reaction_in_rxndb = {}
    rxndb_to_model_total_info = pd.read_csv(rxndb_total_info_to_model_path)
    rxndb_to_model_total_info['equation_dict'] = rxndb_to_model_total_info['equation_dict'].apply(lambda x:literal_eval(x))
    # print(rxndb_to_model_total_info.shape)
    rxndb_to_model_total_info['equation_dict'] = rxndb_to_model_total_info['equation_dict'].apply(lambda x:{k:v for k, v in x.items() if k not in miss_met_id})
    rxndb_to_model_total_info = rxndb_to_model_total_info[rxndb_to_model_total_info['equation_dict'].apply(lambda x : all('s_' in k for k, v in x.items()))]
    # print(rxndb_to_model_total_info.shape)
    get_yeast8_reaction_in_rxndb_partial = partial(get_yeast8_reaction_in_rxndb,rxndb_to_model_total_info=rxndb_to_model_total_info,model_=model_)
    with Pool(num_process) as p:
        yeast8_reaction_in_rxndb = dict(p.imap(get_yeast8_reaction_in_rxndb_partial, reaID_with_smiles_dict.items()))
   
    return yeast8_reaction_in_rxndb

def get_sce_gene2ec_dict(sce_uniprot_path):
    sce_uniprot = pd.read_csv(sce_uniprot_path,sep='\t')
    sce_uniprot = sce_uniprot.fillna('')
    sce_uniprot['Gene Names'] = sce_uniprot['Gene Names'].apply(lambda x:x.split(' '))
    sce_gene2ec_dict = {}
    for index,row in sce_uniprot.iterrows():
        for gene in row['Gene Names']:
            if gene not in sce_gene2ec_dict.keys():
                sce_gene2ec_dict[gene] = []
                sce_gene2ec_dict[gene].append(row['EC number'])
            else:
                sce_gene2ec_dict[gene].append(row['EC number'])

    print(len(sce_gene2ec_dict))
    sce_gene2ec_dict = {k:v for k,v in sce_gene2ec_dict.items() if k.startswith('Q') or k.startswith('Y')}
    sce_gene2ec_dict = {k:[x for x in v if x!=''] for k,v in sce_gene2ec_dict.items()}
    print(len(sce_gene2ec_dict))
    return sce_gene2ec_dict

def filter_gene_with_ec(gene_lst,sce_gene2ec_dict):
    tmp = []
    for gene in gene_lst:
        if len(sce_gene2ec_dict[gene])>0:
            # print(sce_gene2ec_dict[gene])
            tmp.append(gene)
    return tmp

def get_reaction_recovery_ratio(yeast8_reaction_in_rxndb):
    num = 0
    for k,v in yeast8_reaction_in_rxndb.items():
        if len(v)>0:
            num+=1
    print(num/len(yeast8_reaction_in_rxndb))

def check_gpr_accuracy(yeast8_reaction_in_rxndb,model,sce_gene2ec_dict,rxndb_to_model_total_info):
    nogpr = 0
    truegpr = 0
    errorgpr = 0
    true_gene = []
    all_gene = []
    for reaID,rxndb_ID_lst in yeast8_reaction_in_rxndb.items():
        rea_gpr = str(model.reactions.get_by_id(reaID).gpr)
        rea_gpr = rea_gpr.replace('(', '').replace(')', '').replace(' and ', ' ').replace(' or ', ' ').split()
        rea_gpr = filter_gene_with_ec(rea_gpr,sce_gene2ec_dict)

        # print('rea_gpr',rea_gpr)
        rxndb_gpr_lst = []
        for rxndb_ID in rxndb_ID_lst:
            rxndb_gpr = rxndb_to_model_total_info[rxndb_to_model_total_info['NO']==rxndb_ID]['GPR'].to_list()[0]
            rxndb_gpr = rxndb_gpr.replace('(', '').replace(')', '').replace(' and ', ' ').replace(' or ', ' ').split()
            rxndb_gpr_lst += rxndb_gpr
        rxndb_gpr_lst = list(set(rxndb_gpr_lst))
        # print('rxndb_gpr_lst',rxndb_gpr_lst)

        if len(rea_gpr)==0:
            # print('nogpr')
            nogpr +=1
        elif any(gene in rxndb_gpr_lst for gene in rea_gpr):
            # print('truegpr')
            all_gene +=rxndb_gpr_lst
            true_gene += [x for x in rxndb_gpr_lst if x in rea_gpr]
            truegpr+=1
        else:
            errorgpr+=1

    # print('nogpr: ' + str(nogpr),'truegpr: ' + str(truegpr),'errorgpr: ' + str(errorgpr))
    print('ACC',truegpr/(truegpr+errorgpr))
    print('FP',1-len(list(set(true_gene)))/len(list(set(all_gene))))


# def check_gpr_accuracy(yeast8_reaction_in_rxndb,model,sce_gene2ec_dict,rxndb_to_model_total_info):
#     nogpr = 0
#     truegpr = 0
#     errorgpr = 0
#     for reaID,rxndb_ID_lst in yeast8_reaction_in_rxndb.items():
#         rea_gpr = str(model.reactions.get_by_id(reaID).gpr)
#         rea_gpr = rea_gpr.replace('(', '').replace(')', '').replace(' and ', ' ').replace(' or ', ' ').split()
#         rea_gpr = filter_gene_with_ec(rea_gpr,sce_gene2ec_dict)

#         # print('rea_gpr',rea_gpr)
#         rxndb_gpr_lst = []
#         for rxndb_ID in rxndb_ID_lst:
#             rxndb_gpr = rxndb_to_model_total_info[rxndb_to_model_total_info['NO']==rxndb_ID]['GPR'].to_list()[0]
#             rxndb_gpr = rxndb_gpr.replace('(', '').replace(')', '').replace(' and ', ' ').replace(' or ', ' ').split()
#             rxndb_gpr_lst += rxndb_gpr
#         rxndb_gpr_lst = list(set(rxndb_gpr_lst))
#         # print('rxndb_gpr_lst',rxndb_gpr_lst)

#         if len(rea_gpr)==0:
#             # print('nogpr')
#             nogpr +=1
#         elif any(gene in rxndb_gpr_lst for gene in rea_gpr):
#             # print('truegpr')
#             truegpr+=1
#         else:
            
#             errorgpr+=1
#     print('nogpr: ' + str(nogpr),'truegpr: ' + str(truegpr),'errorgpr: ' + str(errorgpr))
#     print(truegpr/(truegpr+errorgpr))

def compare_inchikey0_mnxmeta(inchikey, compare_total_inchikey):
    matched_inchikey = []
    for compare_inchikey in compare_total_inchikey:
        if inchikey== compare_inchikey:
            matched_inchikey.append(compare_inchikey)
    return inchikey, matched_inchikey


def get_smiles2metnetx(mnxmeta_smile_filtered, yeast8_total_inchikey_carbon,yeast8_smiles2metnetx_path,num_processes=60):
    
    res = {}
    compare_smiles_partial = partial(compare_inchikey0_mnxmeta, compare_total_inchikey=yeast8_total_inchikey_carbon)
    with multiprocessing.Pool(num_processes) as p:
    
        for s, t in tqdm(p.imap_unordered(compare_smiles_partial,mnxmeta_smile_filtered,chunksize=100), total=len(mnxmeta_smile_filtered)):
            res[s] = t
    
    # for s, t in tqdm(map(compare_smiles_partial,mnxmeta_smile_filtered), total=len(mnxmeta_smile_filtered)):
    #     res[s] = t
 

    # print(res)

    dump_file(res,yeast8_smiles2metnetx_path)


def smiles_metnetx_reverse(smiles2metnetx_path,total_inchikey0,mnxmeta_smile_inchikey_dict):
    smiles2met = load_pickle(smiles2metnetx_path)
    smiles_lst = []
    for k,v in tqdm(smiles2met.items(),total=len(smiles2met)):
        # if len(v)>0:
        #     print(k,len(v))
        smiles_lst += v

    smiles = {}
    for i in tqdm(total_inchikey0,total=len(total_inchikey0)):
        smiles[i] = []
        for k,v in smiles2met.items():
            if i in v:
                smiles[i].append(k)
    smiles2metnetx = {}
    for k,v in tqdm(smiles.items(),total=len(smiles)):
        smiles_lst = []
        for i in v:
            smiles_lst += mnxmeta_smile_inchikey_dict[i]
        smiles2metnetx[k] = smiles_lst
    dump_file(smiles2metnetx,smiles2metnetx_path)

def merge_metnetx_smiles(yeast_total_smiles,smiles2metnetx_path,yeast_met_file):
    print(len(yeast_total_smiles))
    smiles2metnetx = load_pickle(smiles2metnetx_path)
    for k,v in tqdm(smiles2metnetx.items(),total=len(smiles2metnetx)):
        yeast_total_smiles.extend(v)
    yeast_total_smiles = list(set(yeast_total_smiles))
    dump_file(yeast_total_smiles,yeast_met_file)
    print(len(yeast_total_smiles))

def neutralize_charge(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
 
    # Remove explicit Hs so that they don't interfere with charge neutralization
    mol = Chem.RemoveHs(mol)
    
    # Neutralize charges
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        num_hydrogens = atom.GetTotalNumHs()
        if charge > 0:
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(num_hydrogens - charge if num_hydrogens > charge else 0)
        elif charge < 0:
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(num_hydrogens + abs(charge))
   
    # Convert back to SMILES
    neutral_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
   
    return neutral_smiles

def process_reactions_to_ec_frequency(reactions, gene2ec_dict):
    ec_lst = []
    for rea in reactions:
        tmp = str(rea.gpr).replace(')', '').replace('(', '').replace(' and ', ' ').replace(' or ', ' ').split(' ')
        ec_lst += gene_list2_ec_list(tmp, gene2ec_dict)
    ec_frequency_dict = list_to_frequency_dict(ec_lst)
    return rank_ec_number_dict(ec_frequency_dict)

def process_reactions(reactions):
    gene_lst = []
    for rea in reactions:
        tmp = str(rea.gpr).replace(')', '').replace('(', '').replace(' and ', ' ').replace(' or ', ' ').split(' ')
        gene_lst += tmp
    gene_lst = [x for x in gene_lst if x != '']
    print(len(set(gene_lst)))
    return list_to_frequency_dict(gene_lst)

# def get_cid_from_smiles(smiles):
#     try:
#         # Search for the compound using its SMILES
#         compound = pcp.get_compounds(smiles, 'smiles')
#         time.sleep(1)
#         if compound:
#             return compound[0].cid
#         else:
#             return None
#     except Exception as e:
#         print(f"An error occurred: {e}")
#         return None

def get_cid_from_smiles(smiles):
    try:
        # Search for the compound using its SMILES
        compound = pcp.get_compounds(smiles, 'smiles') or pcp.get_compounds(neutralize_charge(smiles), 'smiles')
        time.sleep(1)
        if compound:
            return compound[0].cid
        else:
            return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None