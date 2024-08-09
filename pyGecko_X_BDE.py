import os
import sys
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from typing import Dict, Iterator
import tensorflow as tf
import tensorflow_addons as tfa
from tensorflow.keras import layers

# Set path to BDEdb2 example predictions
bde_prediction_path = os.path.join(os.getcwd(), 'BDEdb2', 'Example-BDE-prediction')

# Add BDEdb2 prediction path to the Python path if not already present
if bde_prediction_path not in sys.path:
    sys.path.append(bde_prediction_path)

# Import NFP and the CFC preprocessor
import nfp
from preprocess_inputs_cfc import preprocessor

# Defines a custom TensorFlow Keras layer called Slice
class Slice(layers.Layer):
    # Defines the computation at every call of the layer
    def call(self, inputs):
        input_shape = tf.shape(inputs)  # Gets the shape of the input tensor
        num_bonds = input_shape[1] / 2  # Calculates the number of bonds, results in float
        output = tf.slice(inputs, [0, 0, 0], [-1, num_bonds, -1]) # Slices the input tensor to include half the bonds from the second dimension onwards
        output.set_shape(self.compute_output_shape(inputs.shape))   # Sets the shape of the output tensor
        return output

    def compute_output_shape(self, input_shape):    # First dimension (batch size) and the last dimension (features) remain the same
        return [input_shape[0], None, input_shape[2]]   # Second dimension (bonds) is flexible and dynamic

custom_objects = {**nfp.custom_objects, 'Slice': Slice}

# Load the pre-trained model from the specified path
model_path = os.path.join(bde_prediction_path, 'model_3_multi_halo_cfc/best_model.hdf5')
model = tf.keras.models.load_model(model_path, custom_objects=custom_objects)

print("Model loaded")

# Compounds to predict BDEs for
test = np.array(['C', 'C#C/C(C)=C/CNCC', 'CCCCC']) # 'C','CCCCC', 'C#C/C(C)=C/CNCC', 'O=C1NC(=O)NC=C1C', , 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' # 'C#C/C(C)=C/CNCC', ,'CCCCC',

#make the test data graphs

def get_data(smiles):
    ''' Takes a SMILES string and returns a dictionary of input features for the model
    Uses the preprocessor function from preprocess_inputs_cfc.py to generate the input features
    This contains the atom and bond features for the model, as well as the number of atoms and bonds in the molecule'''
    input_dict = preprocessor(smiles)
    input_dict['n_atom'] = len(input_dict['atom'] )
    input_dict['n_bond'] = len(input_dict['bond'] )
    # print(input_dict)
    return input_dict
# Example:
# get_data converts SMILES strings to a dictionary of features
# input_dict = {
#     "atom": [...],  # List of atom features
#     "bond": [...],  # List of bond features
#     "bond_indices": [...],  # Indices for bonds
#     ...
#     "n_atom": len([...]),  # Number of atoms
#     "n_bond": len([...])   # Number of bonds
# }
# These dictionaries serve as input entries for the TensorFlow dataset



# Function to create a TensorFlow dataset from a list of SMILES strings
def create_dataset(smiles_list, batch_size=1000):
    dataset = tf.data.Dataset.from_generator(
        lambda: iter(get_data(smiles) for smiles in smiles_list),
        output_signature={
            **preprocessor.output_signature,
            'n_atom': tf.TensorSpec(shape=(), dtype=tf.int32),
            'n_bond': tf.TensorSpec(shape=(), dtype=tf.int32)
        }
    ).padded_batch(
        batch_size=batch_size,
        padding_values={
            **preprocessor.padding_values,
            'n_atom': tf.constant(0, dtype=tf.int32),
            'n_bond': tf.constant(0, dtype=tf.int32)
        }
    )
    return dataset
# When molecules have varying numbers of features, padding ensures they form a cohesive batch by padding smaller lists in the features with zeros


# Create the dataset for the test SMILES
test_dataset = create_dataset(test)

# Check the structure of the dataset, for debugging purposes
# for data in test_dataset.take(1):
#     print(data)
########



def get_predictions(model, dataset):
    """ Generate predictions using the model and a given dataset """
    return model.predict(dataset, verbose=True)


def organize_predictions(predicted_bdes, test_smiles):
    """ Organize predictions into a DataFrame with molecule and bond indices """
    # Reshape predictions and create a DataFrame
    pred_bdes_flat = predicted_bdes.reshape(-1, 2)
    pred_df = pd.DataFrame(pred_bdes_flat, columns=['pred_bde', 'pred_bdfe'])

    # Assign indices to match test SMILES strings
    num_molecules = predicted_bdes.shape[0]
    num_bonds = predicted_bdes.shape[1]
    molecules_repeated = np.repeat(test_smiles, num_bonds) # Repeat each SMILES string for each bond for which a BDE was predicted
    pred_df['molecule'] = molecules_repeated

    # Function to add bond indices
    def add_bond_indices(df):
        df['bond_index'] = range(num_bonds)
        return df

    # Restructure DataFrame
    pred_df = pred_df.groupby('molecule', group_keys=False).apply(add_bond_indices)

    # Filter non-zero predictions
    pred_df = pred_df[(pred_df['pred_bde'] != 0.000000) & (pred_df['pred_bdfe'] != 0.000000)]

    return pred_df


# Model to predict
predicted_bdes = get_predictions(model, test_dataset)

# Organize predictions
pred_bdes = organize_predictions(predicted_bdes, test)

# Display organized predictions DataFrame (for debugging/verification)
# print(pred_bdes)


# New Helper Class and Functions for Fragmentation

class Molecule:
    ''' Provides methods to generate molceules with explicit hydrogens from Smiles strings '''

    def __init__(self, mol: rdkit.Chem.Mol = None, smiles: str = None) -> None:
        assert (mol is not None) or (smiles is not None), "mol or smiles must be provided"
        self._mol = mol
        self._smiles = smiles
        self._molH = None
        self._is_canon = False

    @property
    def mol(self) -> rdkit.Chem.Mol:
        if self._mol is None:
            self._mol = rdkit.Chem.MolFromSmiles(self._smiles)
        return self._mol

    @property
    def molH(self) -> rdkit.Chem.Mol:
        if self._molH is None:
            self._molH = rdkit.Chem.AddHs(self.mol)
        return self._molH

    @property
    def smiles(self) -> str:
        if (self._smiles is None) or not self._is_canon:
            self._smiles = rdkit.Chem.MolToSmiles(self.mol)
        return self._smiles

def get_fragments(input_molecule: Molecule) -> pd.DataFrame:
    ''' Converts molecules into fragments by iterating over bonds'''
    return pd.DataFrame(fragment_iterator(input_molecule))

def fragment_iterator(input_molecule: Molecule) -> Iterator[Dict]:
    ''' Iterates over bonds in a molecule and fragments it '''
    rdkit.Chem.Kekulize(input_molecule.molH, clearAromaticFlags=True)
    for bond in input_molecule.molH.GetBonds():
        # Usually have high BDE but we want all fragments
        # if bond.IsInRing() or bond.GetBondTypeAsDouble() > 1.9999:
        #     continue
        try:
            mh = rdkit.Chem.RWMol(input_molecule.molH)
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            mh.RemoveBond(a1, a2)
            mh.GetAtomWithIdx(a1).SetNoImplicit(True)
            mh.GetAtomWithIdx(a2).SetNoImplicit(True)
            rdkit.Chem.SanitizeMol(mh)
            fragmented_smiles = rdkit.Chem.MolToSmiles(mh)
            frag1, frag2 = sorted(fragmented_smiles.split("."))
            frag1 = Molecule(smiles=frag1)
            frag2 = Molecule(smiles=frag2)
            yield {
                "molecule": input_molecule.smiles,
                "bond_index": bond.GetIdx(),
                "bond_type": f"{bond.GetBeginAtom().GetSymbol()}-{bond.GetEndAtom().GetSymbol()}",
                "fragment1": frag1.smiles,
                "fragment2": frag2.smiles
            }
        except ValueError:
            continue

# # Corrected Code to Generate Fragment Information and Predictions
# def generate_fragments_and_predictions(test_smiles, predicted_bdes):
#     ''' Iterates over the test SMILES and generates fragments for each molecule
#     Aligns these fragments with the predictions using bond indices
#     Combines all fragments and predicitons into a single DataFrame '''
#     all_fragments = []
#     all_predictions_and_fragments = []
#
#     for idx, smiles in enumerate(test_smiles):
#         molecule = Molecule(smiles=smiles)
#         fragments = get_fragments(molecule)
#
#         # Ensure that the fragments only merge with BDE predictions if the SMILES match
#         matching_fragments = fragments[fragments['molecule'] == smiles]
#
#         # Ensure that the predicted BDEs align with the fragments
#         if len(predicted_bdes[idx][:, 0]) == len(matching_fragments):
#             matching_fragments['pred_bde'] = predicted_bdes[idx][:, 0]
#         else:
#             # Handle the mismatch case appropriately
#             if len(predicted_bdes[idx][:, 0]) > len(matching_fragments):
#                 print(
#                     f"Warning: More predicted BDEs ({len(predicted_bdes[idx][:, 0])}) than fragments ({len(matching_fragments)}) for SMILES {smiles}.")
#                 matching_fragments['pred_bde'] = predicted_bdes[idx][:len(matching_fragments), 0]
#             else:
#                 raise ValueError(
#                     f"Length of predicted BDEs ({len(predicted_bdes[idx][:, 0])}) does not match length of fragments ({len(matching_fragments)})."
#                 )
#
#         all_fragments.append(matching_fragments)
#         all_predictions_and_fragments.append(matching_fragments)
#
#     # Combine all fragments into one DataFrame
#     all_fragments_df = pd.concat(all_fragments, ignore_index=True)
#     return all_fragments_df

#
# def generate_fragments_and_predictions(test_smiles, predicted_bdes):
#     all_fragments = []
#
#     for smiles in test_smiles:
#         molecule = Molecule(smiles=smiles)
#         fragments_df = get_fragments(molecule)
#
#         # Filter predictions to retain only those corresponding to the current molecule
#         molecule_pred_df = predicted_bdes[predicted_bdes['molecule'] == smiles]
#
#         # Ensure `bond_index` alignment between predictions and fragments
#         num_pred_bonds = len(molecule_pred_df)
#         num_frag_bonds = len(fragments_df)
#
#         if num_pred_bonds != num_frag_bonds:
#             print(
#                 f"Warning: Mismatch in number of bonds and predicted BDEs for molecule {smiles}: {num_pred_bonds} predictions vs {num_frag_bonds} fragments.")
#
#         if num_pred_bonds > num_frag_bonds:
#             molecule_pred_df = molecule_pred_df.head(num_frag_bonds)
#         elif num_pred_bonds < num_frag_bonds:
#             fragments_df = fragments_df.head(num_pred_bonds)
#
#         fragments_df['pred_bde'] = molecule_pred_df['pred_bde'].values
#         fragments_df['pred_bdfe'] = molecule_pred_df['pred_bdfe'].values
#
#         all_fragments.append(fragments_df)
#
#     all_fragments_df = pd.concat(all_fragments, ignore_index=True)
#     return all_fragments_df


def reshape_predictions(predicted_bdes, test_smiles):
    num_molecules, num_bonds, num_features = predicted_bdes.shape

    # Flatten predictions
    pred_bdes_flat = predicted_bdes.reshape(-1, num_features)

    # Repeat the molecule names for each bond in predictions
    molecules_repeated = np.repeat(test_smiles, num_bonds)

    # Create a DataFrame
    pred_df = pd.DataFrame(pred_bdes_flat, columns=['pred_bde', 'pred_bdfe'])
    pred_df['molecule'] = molecules_repeated

    # Add bond_index for identification
    bond_indices = np.tile(np.arange(num_bonds), num_molecules)
    pred_df['bond_index'] = bond_indices

    return pred_df


# Reshape predictions into DataFrame
pred_df = reshape_predictions(predicted_bdes, test)


def generate_fragments_and_predictions(test_smiles, predicted_bdes_df):
    all_fragments = []

    # Iterate over each molecule SMILES
    for smiles in test_smiles:
        molecule = Molecule(smiles=smiles)
        fragments_df = get_fragments(molecule)

        # Filter to only current molecule predictions
        # molecule_pred_df = predicted_bdes_df[predicted_bdes_df['molecule'] == smiles]
        # Filter to retain only current molecule predictions and exclude zero predictions
        molecule_pred_df = predicted_bdes_df[(predicted_bdes_df['molecule'] == smiles) & (predicted_bdes_df['pred_bde'] != 0)]
        # Validate counts
        num_pred_bonds = len(molecule_pred_df)
        num_frag_bonds = len(fragments_df)

        if num_pred_bonds != num_frag_bonds:
            print(
                f"Warning: Mismatch in number of bonds and predicted BDEs for molecule {smiles}: {num_pred_bonds} predictions vs {num_frag_bonds} fragments.")

        if num_pred_bonds > num_frag_bonds:
            molecule_pred_df = molecule_pred_df.head(num_frag_bonds)
        elif num_pred_bonds < num_frag_bonds:
            fragments_df = fragments_df.head(num_pred_bonds)

        fragments_df['pred_bde'] = molecule_pred_df['pred_bde'].values
        fragments_df['pred_bdfe'] = molecule_pred_df['pred_bdfe'].values

        all_fragments.append(fragments_df)

    all_fragments_df = pd.concat(all_fragments, ignore_index=True)
    return all_fragments_df



# Set display options to show all columns and full column content
# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_colwidth', None)

# Combine Fragments and Predictions
all_predictions_and_fragments = generate_fragments_and_predictions(test, pred_df)
print(all_predictions_and_fragments)



