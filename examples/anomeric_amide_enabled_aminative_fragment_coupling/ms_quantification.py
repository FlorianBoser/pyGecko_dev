from pygecko.parsers import MS_Base_Parser
from pygecko.parsers.msconvert_wraper import msconvert
from pygecko.reaction import Transformation, Product_Array, Reaction_Array
from pygecko.visualization.visuals import Visualization
from pygecko.analysis.analysis import Analysis
from pathlib import Path
from rdkit import Chem

def main():

    # Define input and output paths.
    ms_raw = '/raw_MS'
    mzml_path = '/mzML'

    layout_path = '/product_smiles.csv'
    layout = Product_Array(layout_path)

    convert_to_mzml(ms_raw, mzml_path)
    ms_sequence = MS_Base_Parser.load_sequence(mzml_path, pos=True)

    # Pick peaks in the GC-MS
    ms_sequence.pick_peaks()
    ms_sequence.set_internal_standard(4.02, name='Dodecane', smiles='CCCCCCCCCCCC')

    # Estimate the yields of the reactions.
    yield_array = Analysis.calc_plate_ms_only_yield(ms_sequence, layout)

    # Generate plate heatmap.
    Visualization.visualize_plate_qualitative(yield_array['quantity'], well_labels=True,
                                  row_labels=["A", "B", "C", "D", "E", "F", "G", "H"],
                                  col_labels=["1", "2", "3", "4", "5", "6" , "7", "8", "9", "10", "11", "12"]
                                              )#, path='/results/Heatmap_Anomeric_Amide_Enabled_Aminative_Fragment_Coupling.png')

    print("Done! Best wishes from the Glorius Group!")

if __name__ == '__main__':
    def convert_to_mzml(ms_raw: str, output_dir: str):
        raw_directory = Path(ms_raw)
        supported_formats = ['.D']
        raw_files = []
        for file_format in supported_formats:
            raw_files.extend(raw_directory.glob(f'*{file_format}'))
        msconvert(input_files = raw_files, output_dir=output_dir, format='mzML')

    main()