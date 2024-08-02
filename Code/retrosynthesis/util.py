from rdkit import Chem
import pandas as pd
import re
import cobra
from tqdm import tqdm
import pubchempy as pcp
import json
from common import *

def get_GEM_all_metabolite_with_smiles(model_path):
    if 'mat' in model_path:
        model = cobra.io.load_matlab_model(model_path)
        metabolite_smiles_dict = {}
        for i in model.metabolites:
            if 'SMILES' in i.annotation:
                try: 
                    smiles_tmp = i.annotation['SMILES'][0]
                    metabolite_smiles_dict[i.name] = normalize_smiles(smiles_tmp)
                except:pass
        return metabolite_smiles_dict
    else: 
        print('error model file')
        return {}


def standardize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except:
        return None


def process_ymdb(ymdb_path, output):
    df = pd.read_csv(ymdb_path)
    for i in range(len(df)):
        if df['super_class'][i] != 'Lipids and lipid-like molecules':
            continue
        else:
            df['NAME'][i] = re.sub(r'\(\d+z\)','',df['NAME'][i])
    df['standard_smiles'] = df['SMILES'].apply(standardize_smiles)
    df.to_csv(output, index=False)

def check_metabolite_in_model(yeast_gem_final_path, ymdb_path, output_excel):
    yeast_gem_final = pd.read_csv(yeast_gem_final_path)
    ymdb = pd.read_csv(ymdb_path)
    for i in range(len(ymdb)):
        if ymdb['super_class'][i] != 'Lipids and lipid-like molecules':
            continue
        else:
            ymdb['NAME'][i] = re.sub(r'\(\d+z\)','',ymdb['NAME'][i])
    ymdb['standard_smiles'] = ymdb['SMILES'].apply(standardize_smiles)
    for i in range(len(yeast_gem_final)):
        if '1-mlcl' in yeast_gem_final['NAME'][i]:
            yeast_gem_final['NAME'][i] = yeast_gem_final['NAME'][i].split('(')[0] + '(0:0/' + yeast_gem_final['NAME'][i].split('(')[1]
        elif '2-mlcl' in yeast_gem_final['NAME'][i]: 
            yeast_gem_final['NAME'][i] = yeast_gem_final['NAME'][i].split(')')[0] + '/0:0)'
    yeast_gem_final['inchikey0'] = yeast_gem_final['SMILES'].apply(smiles2inchikey0)
    ymdb['inchikey0'] = ymdb['SMILES'].apply(smiles2inchikey0)
    ymdb['in_model'] = 0
    for i in tqdm(range(len(ymdb))):
        
        if  ymdb.loc[i,'NAME'] in yeast_gem_final['NAME'].tolist():
            ymdb.loc[i,'in_model'] = 1
    for i in tqdm(range(len(ymdb))):

        if  ymdb.loc[i,'inchikey0'] in yeast_gem_final['inchikey0'].tolist():
            ymdb.loc[i,'in_model'] = 1
    ymdb.to_excel(output_excel,index=False)
    yeast_gem_final.to_csv(yeast_gem_final_path,index=False)
    ymdb = ymdb[ymdb['in_model'] == 1]
    print('in model: ',len(ymdb))

def extract_mnxm_id(s):
    
    if re.search(r'(MNXM\d+)',s):
        return re.search(r'(MNXM\d+)',s).group(1)
    else:
        return None
def extract_chebi_id(s):
    if re.search(r'(CHEBI:\d+)',s):
        return re.search(r'(CHEBI:\d+)',s).group(1)
    else:
        return None
def extract_kegg_id(s):
    if re.search(r'(C\d+)',s):
        return re.search(r'(C\d+)',s).group(1)
    else:
        return None
    
def add_smiles_pubchem(smiles,compound_name):
    if pd.isnull(smiles) :
        try:
            compound = pcp.get_compounds(compound_name,'name')[0]
            return compound.isomeric_smiles
            
        except:
            print(f"Error processing {compound_name}")

            return None
    else:
        return smiles
    

def extract_smiles_and_ids_secondary(suppl):
    # Initialize dataframes
    

    # Read the SDF file
    secondary_smiles_list = []
    secondary_id_list = []
   

    # Iterate over the molecules and extract properties
    for mol in tqdm(suppl,total=len(suppl)):
        try:
            if mol is not None:
                smiles = mol.GetProp('SMILES')
                secondary_id = mol.GetProp('Secondary ChEBI ID')
                secondary_smiles_list.append(smiles)
                secondary_id_list.append(secondary_id)
        except:
            continue

    return secondary_smiles_list, secondary_id_list
def extract_smiles_and_ids(suppl):
    # Initialize dataframes
    

    # Read the SDF file
    
    smiles_list = []
    id_list = []

    # Iterate over the molecules and extract properties
    

    for mol in tqdm(suppl,total=len(suppl)):
        try:
            if mol is not None:
                smiles = mol.GetProp('SMILES')
                id = mol.GetProp('ChEBI ID')
                smiles_list.append(smiles)
                id_list.append(id)
        except:
            continue
    return smiles_list, id_list

def update_metMetaNetXID(id,id_list_old,id_list_new):
    
    if id in id_list_old:
        index = id_list_old.index(id)
        return id_list_new[index]
    else:
        return None
    
def update_smiles( smiles,id_column_list,smiles_list,smiles_column):
    if pd.isnull(smiles) :
        try:
            if pd.notnull(smiles_column) :
              
                index = id_column_list.index(smiles_column)
                return smiles_list[index]
            
        except:
            return None
    else:
        return smiles



def process_lipid_names(df):
    dict1 = {'triglyceride':'tg','CDP-diacylglycerol':'cdp-dg','phosphatidyl-N,N-dimethylethanolamine':'pe-nme2'
             ,'phosphatidyl-N-methylethanolamine':'pe-nme','phosphatidylglycerol':'pg','monolysocardiolipin (2':'2-mlcl (2',
             'monolysocardiolipin (1':'1-mlcl (1','1-phosphatidyl-1D-myo-inositol (1':'pi (1','cardiolipin':'cl','phosphatidylethanolamine':'pe'
             ,'phosphatidylcholine':'pc','3-(3-sn-phosphatidyl)-sn-glycerol1-phosphate':'pgp'}
    for i in range(len(df)):
        if df['SMILES'].isna()[i] == False and '*' not in df['SMILES'][i]:
            continue
        else:
            for key in dict1:
                if key in df['NAME'][i]:
                    df['NAME'][i] = df['NAME'][i].replace(key,dict1[key])
    for i in range(len(df)):
        if df['SMILES'].isna()[i] == False and '*' not in df['SMILES'][i]:
            continue
        else:
            if re.search(r'\(1-\d+:\d+',df['NAME'][i]):
                df['NAME'][i] = re.sub(r'\(1-','(',df['NAME'][i])
                df['NAME'][i] = re.sub(r',\s+2-','/',df['NAME'][i])
                if re.search(r',\s+3-',df['NAME'][i]):
                    df['NAME'][i] = re.sub(r',\s+3-','/',df['NAME'][i])
                if re.search(r',\s+4-',df['NAME'][i]):
                    df['NAME'][i] = re.sub(r',\s+4-','/',df['NAME'][i])
            elif re.search(r'\(2-\d+:\d+',df['NAME'][i]):
                df['NAME'][i] = re.sub(r'\(2-','(',df['NAME'][i])
                df['NAME'][i] = re.sub(r',\s+3-','/',df['NAME'][i])
                df['NAME'][i] = re.sub(r',\s+4-','/',df['NAME'][i])
            df['NAME'][i] = re.sub(r'\s+','',df['NAME'][i])
    return df

def sort_fatty_acid_chains(df, output_csv):
    list1 = ['tg','cdp-dg','pi','cl','pe-nme2','pe-nme','pg','2-mlcl','1-mlcl','pc','pe','pgp']
    for i in range(len(df)):
        if df['SMILES'].isna()[i] == False:
            continue
        else:
            if df['NAME'][i].split('(')[0] in list1:
                string_ = df['NAME'][i].split('(')[0]
                # Extract the fatty acid chains
                fatty_acid_chains = re.findall(r'\d+:\d+', df['NAME'][i])

                # Sort the fatty acid chains
                sorted_fatty_acid_chains = sorted(fatty_acid_chains)

                # Reconstruct the TG string with the sorted fatty acid chains
                sorted_tg_string = string_ + '(' + '/'.join(sorted_fatty_acid_chains) + ')'
                df['NAME'][i] = sorted_tg_string
    df.to_csv(output_csv, index=False)

def get_smiles_from_model(yeast_gem_path, yeast_gem_smiles):
    yeast_gem = pd.read_excel(yeast_gem_path, sheet_name='METS')
    with open(yeast_gem_smiles, 'r') as f:
        metabolite_smiles_dict = json.load(f)
    yeast_gem['smiles'] = yeast_gem['NAME'].apply(lambda x: metabolite_smiles_dict.get(x, None))
    return yeast_gem

def drop_uncorrect_smiles(yeast_gem):
    for index, row in tqdm(yeast_gem.iterrows()):
        if not pd.isnull(row['smiles']):
            if not compare_atoms_between_formula_smiles(row['smiles'], row['COMPOSITION']):
                yeast_gem.at[index, 'smiles'] = None
    return yeast_gem