from pathlib import Path

import xarray as xr
from rdkit import Chem
from rdkit.Chem import Descriptors

from pygecko.reaction.layout import Combinatorial_Layout, Product_Layout
from pygecko.reaction.transformation import Transformation

# FBS Implementation for Fragmentation Search
# Only import BDE fragments if the script and model is available
# bde_module = importlib.util.find_spec("pygecko.BDE_fragments")
# if bde_module:
#     import pygecko.BDE_fragments as BDE
from pygecko.fragments_bde import fragments_bde
import subprocess
import tempfile
import json
import pandas as pd
import numpy as np

class Reaction_Array(Combinatorial_Layout):

    '''
    A class for storing information about the combinatorial layout of a reaction array.

    Attributes:
        design (xr.DataArray): xarray DataArray containing the combinatorial design of the reaction array.
        rows (dict): Dictionary containing the mapping of the reaction array rows to the combinatorial dimensions of the
        reaction array.
        columns (dict): Dictionary containing the mapping of the reaction array columns to the combinatorial dimensions
        ofthe reaction array.
    '''

    def __init__(self, layout_file:Path|str, transformation:Transformation, meta_data_file:str|None=None, calculate_fragments:bool=False):
        super().__init__(layout_file, meta_data_file, transformation)
        self.design = xr.DataArray(self.array, dims=['x', 'y'], coords={'x': list(map(chr, range(65, 65+self.array.shape[0]))), 'y': list(range(1, self.array.shape[1]+1))}, attrs=self.meta_data)
        self.rows = {self.design.x.values[i]: self.layout.x.dropna().to_numpy(dtype=str)[i] for i in range(len(self.design.x.values))}
        self.columns = {self.design.y.values[i]: self.layout.y.dropna().to_numpy(dtype=str)[i] for i in range(len(self.design.y.values))}
        if meta_data_file:
            self.__extend_metadata()
        if calculate_fragments:
            self.fragments = self.calculate_fragments()


    def __getitem__(self, pos:str) -> list|str:

        '''
        Takes a well position and returns the substrates in this position.
        '''

        item = self.design.loc[pos[0], int(pos[1:])].item()
        if '.' in item:
            item = item.split('.')
        return item

    def get_product_mw(self, pos:str):

        '''
        Returns the molecular weight of the product in the given position.

        Args:
            pos (str): Well position.

        Returns:
            float: Molecular weight of the product.

        '''

        product = self.get_product(pos)
        mol = Chem.MolFromSmiles(product)
        mw = Descriptors.ExactMolWt(mol)
        return mw

    def get_product_mz(self, pos:str):

        '''
        Returns the mass to charge ratio of the single charged product in the given position.
        Args:
            pos (str): Well position.

        Returns:
            float: Mass to charge ratio of the single charged product.
        '''

        mw = self.get_product_mw(pos)
        return round(mw, 0)


    def get_product(self, pos):

        '''
        Returns the product of the reaction in the given position.

        Args:
            pos (str): Well position.

        Returns:
            str: SMILES string of the product.

        '''

        subst = self[pos]
        product = self.transformation(subst)
        return product

    def get_substrate(self, pos, index=0):
        subst = self[pos]
        if isinstance(subst, list):
            subst = subst[index]
        return subst


    def __extend_metadata(self):

        '''
        Extends the metadata of the layout by replacing the '$row$' and '$column$' entries with the actual compounds
        from the layout.
        '''

        new_entries = []
        del_entries = []
        for i, stock in enumerate(self.meta_data['stock_solutions']):
            if stock['compound'] == '$row$':
                for well in stock['wells']:
                    new_entry = stock.copy()
                    new_entry['compound'] = self.rows[well]
                    new_entry['wells'] = [well]
                    new_entries.append(new_entry)
                    del_entries.append(i)
            if stock['compound'] == '$column$':
                for well in stock['wells']:
                    new_entry = stock.copy()
                    new_entry['compound'] = self.columns[int(well)]
                    new_entry['wells'] = [well]
                    new_entries.append(new_entry)
                    del_entries.append(i)
        for i in sorted(set(del_entries), reverse=True):
            del self.meta_data['stock_solutions'][i]
        self.meta_data['stock_solutions'].extend(new_entries)

    def get_fragments(self, pos: str) -> dict:  # CHANGE: Added method to get fragments
        '''
        Returns the fragments of the product in the given position.
        Args:
            pos (str): Well position.
        Returns:
            dict: Dictionary of fragment SMILES strings and their respective BDEs.
        '''
        product_smiles = self.get_product(pos)
        can_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(product_smiles))
        return self.fragments.get(can_smiles, {})

    def calculate_fragments(self) -> dict:  # CHANGE: Added fragment calculation method
        '''
        Calculate fragments for all products and substrates.
        Returns:
            dict: Fragments for all molecules.
        '''
        #all_smiles = [self.get_product(pos) for pos in self.design.coords["x"].values + self.design.coords["y"].values]
        all_smiles = [self.get_product(f'{x}{y}') for x in self.design.coords["x"].values for y in self.design.coords["y"].values]  # CHANGE: Corrected to loop through all well positions
        all_smiles = list(set(all_smiles))  # Remove duplicates
        result_df = import_bde_fragments(all_smiles)
        fragments_dict = self.get_fragments_dict(result_df)
        return fragments_dict

    def get_fragments_dict(self,
                           fragments_df: pd.DataFrame) -> dict:  # CHANGE: Added method to convert fragments dataframe to dict
        '''
        Converts the fragments dataframe to a dictionary, considering the fragments with the lowest pred_bde first.
        Args:
            fragments_df (pd.DataFrame): Dataframe containing the fragments.
        Returns:
            dict: Dictionary with fragment information.
        '''
        fragments_dict = {}

        # Group by molecule to handle each molecule separately
        for mol, group in fragments_df.groupby('molecule'):
            sorted_group = group.sort_values(by='pred_bde')  # Sort by pred_bde to consider the lowest values first
            fragments_dict[mol] = {}

            for _, row in sorted_group.iterrows():
                fragment1 = row['fragment1']  # CHANGE: Use fragment1
                fragment2 = row['fragment2']  # CHANGE: Use fragment2
                bde = row['pred_bde']

                # Add fragment1 if not already present
                if fragment1 not in fragments_dict[mol]:
                    fragments_dict[mol][fragment1] = bde

                # Add fragment2 if not already present
                if fragment2 not in fragments_dict[mol]:
                    fragments_dict[mol][fragment2] = bde

        return fragments_dict

    def get_product_fragments(self, pos: str) -> dict:  # CHANGE: Added method to get fragments
        '''
        Returns the fragments of the product in the given position.
        '''
        return self.fragments.get(self.get_product(pos), {})


class Product_Array(Product_Layout):

    '''
    A class for storing information about the products of a reaction array.

    Attributes:
        design (xr.DataArray): xarray DataArray containing the product design of the reaction array.
    '''

    def __init__(self, layout_file:Path|str, calculate_fragments:bool=False):
        super().__init__(layout_file)
        self.design = xr.DataArray(self.array, dims=['x', 'y'],
                                   coords={'x': list(map(chr, range(65, 65 + self.array.shape[0]))),
                                           'y': list(range(1, self.array.shape[1] + 1))})
        # FBS addition of fragmentation calculation
        if calculate_fragments:
            self.fragments = self.calculate_fragments()

        print('Product array created.')

    def __getitem__(self, pos:str) -> list|str:

        '''
        Takes a well position and returns the substrates in this position.
        '''

        item = self.design.loc[pos[0], int(pos[1:])].item()
        if '.' in item:
            item = item.split('.')
        return item

    def get_product(self, pos):

        '''
        Returns the product of the reaction in the given position.

        Args:
            pos (str): Well position.

        Returns:
            str: SMILES string of the product.

        '''
        return self[pos]

    def get_product_mw(self, pos: str):
        '''
        Returns the molecular weight of the product in the given position.

        Args:
            pos (str): Well position.

        Returns:
            float: Molecular weight of the product.

        '''

        product = self[pos]
        mol = Chem.MolFromSmiles(product)
        mw = Descriptors.ExactMolWt(mol)
        return mw

    def get_product_mz(self, pos: str):
        '''
        Returns the mass to charge ratio of the single charged product in the given position.
        Args:
            pos (str): Well position.

        Returns:
            float: Mass to charge ratio of the single charged product.
        '''

        mw = self.get_product_mw(pos)
        return round(mw, 0)
    
    def get_fragments(self, pos: str) -> dict:  # CHANGE: Added method to get fragments
        '''
        Returns the fragments of the product in the given position.
        Args:
            pos (str): Well position.
        Returns:
            dict: Dictionary of fragment SMILES strings and their respective BDEs.
        '''
        product_smiles = self.get_product(pos)
        can_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(product_smiles))
        return self.fragments.get(can_smiles, {})
    
    

    def calculate_fragments(self) -> dict:  # CHANGE: Added fragment calculation method
        '''
        Calculate fragments for all products and substrates.
        Returns:
            dict: Fragments for all molecules.
        '''
        #all_smiles = [self.get_product(pos) for pos in self.design.coords["x"].values + self.design.coords["y"].values]
        all_smiles = [self.get_product(f'{x}{y}') for x in self.design.coords["x"].values for y in self.design.coords["y"].values]  # CHANGE: Corrected to loop through all well positions
        all_smiles = list(set(all_smiles))  # Remove duplicates
        result_df = import_bde_fragments(all_smiles)
        fragments_dict = self.get_fragments_dict(result_df)
        return fragments_dict

    def get_fragments_dict(self,
                           fragments_df: pd.DataFrame) -> dict:  # CHANGE: Added method to convert fragments dataframe to dict
        '''
        Converts the fragments dataframe to a dictionary, considering the fragments with the lowest pred_bde first.
        Args:
            fragments_df (pd.DataFrame): Dataframe containing the fragments.
        Returns:
            dict: Dictionary with fragment information.
        '''
        fragments_dict = {}

        # Group by molecule to handle each molecule separately
        for mol, group in fragments_df.groupby('molecule'):
            sorted_group = group.sort_values(by='pred_bde')  # Sort by pred_bde to consider the lowest values first
            fragments_dict[mol] = {}

            for _, row in sorted_group.iterrows():
                fragment1 = row['fragment1']  # CHANGE: Use fragment1
                fragment2 = row['fragment2']  # CHANGE: Use fragment2
                bde = row['pred_bde']

                # Add fragment1 if not already present
                if fragment1 not in fragments_dict[mol]:
                    fragments_dict[mol][fragment1] = bde

                # Add fragment2 if not already present
                if fragment2 not in fragments_dict[mol]:
                    fragments_dict[mol][fragment2] = bde

        return fragments_dict

    def get_product_fragments(self, pos: str) -> dict:  # CHANGE: Added method to get fragments
        '''
        Returns the fragments of the product in the given position.
        '''
        return self.fragments.get(self.get_product(pos), {})



def import_bde_fragments(smiles_list: list) -> pd.DataFrame | None:
    '''
    Integrates the BDE fragmentation model to match mass spectrometry peaks based on predicted fragments.

    Args:
        smiles_list (list): List of SMILES strings of the analytes.

    Returns:
        pd.DataFrame | None: DataFrame with the BDE fragmentation results or None if an error occurs.
    '''
    try:
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".json") as temp_file:
            json.dump(smiles_list, temp_file)
            temp_file.flush()

            # Define the shell script to run
            script_path = 'C:/Users/flori/anaconda3/envs/pyGecko/Lib/site-packages/pygecko/fragments_bde/run_fragmentation_bde.bat'  # Use the .bat script if on Windows
            result = subprocess.run(
                [script_path, temp_file.name],
                capture_output=True,
                text=True,
                shell=True  # Required to run the script through the shell
            )
            if result.returncode != 0:
                print("Error in BDE prediction subprocess:", result.stderr)
                return None

            try:
                # Extract clean JSON data from the stdout
                json_data_start = result.stdout.find("[{")
                json_data_end = result.stdout.rfind("}]") + 2
                if json_data_start == -1 or json_data_end == -1:
                    raise json.JSONDecodeError("Could not find JSON data in the subprocess output", result.stdout,
                                               0)

                json_data = result.stdout[json_data_start:json_data_end]
                data_frame = pd.read_json(json_data)
                return data_frame

            except json.JSONDecodeError as e:
                print(f"Failed to decode JSON from BDE prediction output: {str(e)}")
                return None

    except Exception as e:
        print(f"Error in BDE prediction: {str(e)}")
        return None


# ### FBS ###
# ### Fragmentation match with BDE fragments ###
# def match_fragments(smiles_list: list) -> pd.DataFrame | None:
#     '''
#     Integrates the BDE fragmentation model to match mass spectrometry peaks based on predicted fragments.
# 
#     Args:
#         smiles_list (list): List of SMILES strings of the analytes.
# 
#     Returns:
#         pd.DataFrame | None: DataFrame with the BDE fragmentation results or None if an error occurs.
#     '''
#     try:
#         with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".json") as temp_file:
#             json.dump(smiles_list, temp_file)
#             temp_file.flush()
# 
#             # Define the shell script to run
#             script_path = 'C:/Users/flori/anaconda3/envs/pyGecko/Lib/site-packages/pygecko/fragments_bde/run_fragmentation_bde.bat'  # Use the .bat script if on Windows
#             result = subprocess.run(
#                 [script_path, temp_file.name],
#                 capture_output=True,
#                 text=True,
#                 shell=True  # Required to run the script through the shell
#             )
#             if result.returncode != 0:
#                 print("Error in BDE prediction subprocess:", result.stderr)
#                 return None
# 
#             try:
#                 # Extract clean JSON data from the stdout
#                 json_data_start = result.stdout.find("[{")
#                 json_data_end = result.stdout.rfind("}]") + 2
#                 if json_data_start == -1 or json_data_end == -1:
#                     raise json.JSONDecodeError("Could not find JSON data in the subprocess output", result.stdout,
#                                                0)
# 
#                 json_data = result.stdout[json_data_start:json_data_end]
#                 data_frame = pd.read_json(json_data)
#                 return data_frame
# 
#             except json.JSONDecodeError as e:
#                 print(f"Failed to decode JSON from BDE prediction output: {str(e)}")
#                 return None
# 
#     except Exception as e:
#         print(f"Error in BDE prediction: {str(e)}")
#         return None
# 
# def analyze_gc_ms_spectra(smiles_list):
#     # Get BDE predictions and fragment data
#     bde_fragments = BDE_fragments.main(smiles_list)
#     print("BDE Fragments:", bde_fragments)
#     # Proceed with GC/MS analysis using the received BDE predictions and fragments
# 
# # Example usage
# if __name__ == '__main__':
#     smiles_list = ["C", "CCCCC", "C#C/C(C)=C/CNCC"]
#     # Provide minimal example metadata
#     # example_metadata = {
#     #     'AcqMethodName': 'ExampleMethod',
#     #     'another_metadata_field': 'example_value'
#     # }
#     # example_chromatogram = np.array([[0.0, 1.0, 2.0]])
#     # Create an instance of MS_Injection
#     #injection = MS_Injection(metadata=example_metadata, chromatogram=example_chromatogram, peaks=None, scans=pd.DataFrame())
# 
#     result = match_fragments(smiles_list) #injection.
#     if result is not None:
#         print(result)
# 
#     # match_fragments(self=self, smiles_list=smiles_list)
# 
# # TODO: bde_fragments should only be a list of fragments
# # TODO: Run it in the pyGecko environment with the BDE staying in the BDE env