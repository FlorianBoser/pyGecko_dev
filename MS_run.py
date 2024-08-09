from BDE_fragments import BDE_fragments


def analyze_gc_ms_spectra(smiles_list):
    # Get BDE predictions and fragment data
    bde_fragments = BDE_fragments.main(smiles_list)
    print("BDE Fragments:", bde_fragments)
    # Proceed with GC/MS analysis using the received BDE predictions and fragments

# Example usage
if __name__ == '__main__':
    smiles_list = ["C", "CCCCC", "C#C/C(C)=C/CNCC"]
    analyze_gc_ms_spectra(smiles_list)

# TODO: bde_fragments should only be a list of fragments
# TODO: Run it in the pyGecko environment with the BDE staying in the BDE env



# from pygecko.parsers import Agilent_FID_Parser, MS_Base_Parser
# from pygecko.parsers.msconvert_wraper import msconvert
# from pygecko.gc_tools import load_sequence
# from pygecko.reaction import Transformation, Product_Array, Reaction_Array
# from pygecko.visualization.visuals import Visualization
# from pygecko.analysis.analysis import Analysis
# from pygecko.reaction import Reaction_Parser
# from pathlib import Path
#
#
# from rdkit import Chem
# from rdkit.Chem import AllChem
# def main():
#
#
#
#     # # APT Bufunctional Plate
#     # # TODO: Raw data file should be tested afterwards
#     ms_raw = 'C:/Users/flori/Doktorarbeit/08_HTE_aaS/PyGecko/data/JSN_JLT_Plate_GC_Data/Rr_MS'
#     mzml_path = 'C:/Users/flori/Doktorarbeit/08_HTE_aaS/PyGecko/data/JSN_JLT_Plate_GC_Data/Rr_mzML'
#
#     layout_path = 'C:/Users/flori/Doktorarbeit/08_HTE_aaS/PyGecko/data/JSN_JLT_Plate_GC_Data/products_H-HCl.csv'
#
#     # Reaction SMARTS used to map substrates to products.
#     # rxn = Transformation(
#     #     '[c:1]1([a:2][a:3][a:4][a:5]1)[Cl,Br:6].[SH1:7][#6:8]>>[c:1]1([a:2][a:3][a:4][a:5]1)[S:7][#6:8]')
#
#     # Create a Well_Plate object from the reaction layout and reaction SMARTS.
#     layout = Product_Array(layout_path)
#
#     # Convert the mzML files to mzXML files.
#     # convert_to_mzml(ms_raw, mzml_path)
#     ms_sequence = MS_Base_Parser.load_sequence(mzml_path, pos=True)
#
#     # Pick peaks in the GC-MS
#     ms_sequence.pick_peaks() # height=28000, prominence_ms = 1
#     # TODO: Standard peak is currently not labeled, but in the MS only mode probably no standard will be used
#     ms_sequence.set_internal_standard(4.02, name='Dodecane', smiles='CCCCCCCCCCCC')
#
#
#     # Calculate the yields of the reactions.
#     yield_array = Analysis.calc_plate_ms_only_yield(ms_sequence, layout)    # TODO: Ignore the Isotope Detection for Regioisomeres?
#
#     import pandas as pd
#     pd.DataFrame(yield_array['quantity']).to_excel('C:/Users/flori/Doktorarbeit/08_HTE_aaS/PyGecko/data/results/v4_JSN_JLT_Plate_GC_Data.xlsx', index=False, header=False)
#
#     # for i in fid_sequence:
#     #     Visualization.view_chromatogram(i)
#     # Generate plate heatmap.
#     Visualization.visualize_plate_qualitative(yield_array['quantity'], well_labels=True,
#                                   row_labels=["A", "B", "C", "D", "E", "F", "G", "H"],
#                                   col_labels=["1", "2", "3", "4", "5", "6" , "7", "8", "9", "10", "11", "12"]
#                                               )#, path='C:/Users/flori/Doktorarbeit/08_HTE_aaS/PyGecko/data/results/v4_FBS_CSN_JLT_reaction_MS_quant.png')    # , path='C:/Users/flori/Doktorarbeit/08_HTE_aaS/PyGecko/data/yieldmap_APT_FBS_v1.png'
#
#
#
#
#
#
#
#
#     print("Yuchee!")
#
#
#
#
#
#
#
# if __name__ == '__main__':
#     def convert_to_mzml(ms_raw: str, output_dir: str):
#         raw_directory = Path(ms_raw)
#         supported_formats = ['.D']
#         raw_files = []
#         for file_format in supported_formats:
#             raw_files.extend(raw_directory.glob(f'*{file_format}'))
#         msconvert(input_files = raw_files, output_dir=output_dir, format='mzML')
#
#     def test_smarts():
#         column_smarts = "c1cc(!@[S](=O))ccc1"
#         row_smarts = "c1cccc(c1)S(=O)(=O)ONCC(Cl)[c,C]"
#         # Define the SMARTS pattern for alkene and bifunctional compound
#         column_pattern = Chem.MolFromSmarts(column_smarts)
#         row_pattern = Chem.MolFromSmarts(row_smarts)
#
#         column = { "1" : "CC1=CC(C)=CC(C)=C1",
#                     "2" : "C1(C=CC=C2)=C2C=CC=C1",
#                     "3" : "COC1=CC2=C(SC(C)=N2)C=C1",
#                     "4" : "CC(N1C)=CC2=C1C=CC=C2",
#                     "5" : "C1(C(C=CC=C2)=C2O3)=C3C=CC=C1",
#                     "6" : "C1(N=C(C=CC=C2)C2=C3)=C3C=CC=C1",
#                     "7" : "O=C1N(C)C2=C(C=CC=C2)C1",
#                     "8" : "CC(SC1=CC=CC=C1)C",
#                     "9" : "COC1=CC=CC=C1Cl",
#                     "10" : "COC1=CC=CC(Br)=C1",
#                     "11" : "OC1=C(C)C=CC=C1C",
#                     "12" : "NC1=C(C)C=CC=C1C"
#                     }
#
#         rows = {"A" : "ClC(CN(C(OC(C)(C)C)=O)OS(C1=CC=C(C)C=C1)(=O)=O)C2=CC=CC=C2",
#                 "B" : "ClC(CN(C(OC(C)(C)C)=O)OS(C1=CC=C(C)C=C1)(=O)=O)C2=CC=C(OC)C=C2",
#                 "C" : "ClC(CN(C(OC(C)(C)C)=O)OS(C1=CC=C(C)C=C1)(=O)=O)C2=CC=C(C(F)(F)F)C=C2",
#                 "D" : "ClC(CN(C(OC(C)(C)C)=O)OS(C1=CC=C(C)C=C1)(=O)=O)(C)C2=CC=CC=C2",
#                 "E" : "ClC(C(C)N(C(OC(C)(C)C)=O)OS(C1=CC=C(C)C=C1)(=O)=O)C2=CC=CC=C2",
#                 "F" : "ClC(CN(C(OC(C)(C)C)=O)OS(C1=CC=C(C)C=C1)(=O)=O)C2=CC=CS2",
#                 "G" : "ClC(CN(C(OC(C)(C)C)=O)OS(C1=CC=C(C)C=C1)(=O)=O)C2=NC=CC=C2",
#                 "H" : "ClC(CN(C(OC(C)(C)C)=O)OS(C1=CC=C(C)C=C1)(=O)=O)CCC2=CC=CC=C2"
#                 }
#
#         for key, value in rows.items():
#
#             row_mol = Chem.MolFromSmiles(value)
#             if row_mol.HasSubstructMatch(row_pattern) == False:
#                 print(f"Not a compound of the rows substance class: {key}")
#
#                 #break
#
#
#             for k, v in column.items():
#                 column_mol = Chem.MolFromSmiles(v)
#                 if column_mol.HasSubstructMatch(column_pattern) == False:
#                     print(f"Not a compound of the row substance class: {k}")
#
#                     #break
#
#                 educts = [row_mol, column_mol]
#                 try:
#                     products = rxn.transform.RunReactants(educts)
#                     p = Chem.MolToSmiles(products[0][0])
#                     well = key + k
#                     print(f"{well} : {p}")
#                 except:
#                     well = key + k
#                     print(f"{well} : X")
#
#     # TODO: Implement chlorinate_products
#     # def chlorinate_products():
#     #     import pandas as pd
#     #     products_df = pd.read_csv('C:/Users/flori/Doktorarbeit/08_HTE_aaS/PyGecko/data/JSN_JLT_Plate_GC_Data/products.csv')
#     #     p =
#
#
#
#     # test_smarts()
#
#     # chlorinate_products()
#
#     main()