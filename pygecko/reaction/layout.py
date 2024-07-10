import pandas as pd
import numpy as np
from pathlib import Path

from pygecko.reaction.transformation import Transformation
from pygecko.reaction.utilities import read_json

class Combinatorial_Layout:

    '''
    A class for storing information about a combinatorial reaction layout.

    Attributes:
        layout (pd.DataFrame): DataFrame containing the combinatorial dimensions of the layout.
        transformation (Transformation): Transformation object.
        meta_data (dict): Dictionary containing the metadata of the layout.
        dimensions (int): Number of dimensions of the layout.
        array (np.ndarray): A numpy array containing the layout.
    '''

    layout: pd.DataFrame
    transformation: Transformation
    meta_data: dict
    dimensions:int
    array: np.ndarray

    __slots__ = 'layout', 'transformation', 'meta_data', 'dimensions', 'array'

    def __init__(self, layout_file: Path|str, meta_data_file:Path|str, transformation:Transformation):
        self.layout = pd.read_csv(layout_file)
        self.transformation = transformation
        if meta_data_file:
            self.meta_data = read_json(meta_data_file)
        else:
            self.meta_data = {}
        self.dimensions = self.layout.shape[1]
        arr1 = self.layout.x.dropna().to_numpy(dtype=str)
        arr2 = self.layout.y.dropna().to_numpy(dtype=str)
        self.array = np.array([[x + '.' + y for y in arr2] for x in arr1])


class Product_Layout:

    '''
    A class for storing information about a product layout.

    Attributes:
        layout (pd.DataFrame): DataFrame containing the product layout.
        dimensions (int): Number of dimensions of the layout.
        array (np.ndarray): A numpy array containing the layout.
    '''

    layout: pd.DataFrame
    dimensions:int
    array: np.ndarray

    __slots__ = 'layout', 'dimensions', 'array'

    def __init__(self, layout_file: Path|str):
        # Checking if the file exists
        if not Path(layout_file).is_file():
            raise FileNotFoundError(f"No such file: {layout_file}")
        # Attempt to read CSV with standard parameters and additional diagnostics
        try:
            self.layout = pd.read_csv(layout_file, header=None, delimiter=';')  # Assuming no header
        except Exception as e:
            # If reading fails, raise an informative error
            raise IOError(f"Failed to read the file: {e}")

        #self.layout = pd.read_csv(layout_file)
        self.dimensions = self.layout.shape[1]
        self.array = self.layout.to_numpy(dtype=str)
        # print('done')