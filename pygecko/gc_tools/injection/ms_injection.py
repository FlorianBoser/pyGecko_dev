from collections import defaultdict

import numpy as np
import pandas as pd
from brainpy import isotopic_variants
from rdkit import Chem
from rdkit.Chem import Descriptors

from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.analysis import Analysis_Settings
from pygecko.gc_tools.injection import Injection
from pygecko.gc_tools.peak import Peak_Detection_MS, MS_Peak

# Only import BDE fragments if the script and model is available
# bde_module = importlib.util.find_spec("pygecko.BDE_fragments")
# if bde_module:
#     import pygecko.BDE_fragments as BDE
from pygecko.fragments_bde import fragments_bde
import subprocess
import tempfile
import json




class MS_Injection(Injection):

    '''
    Class to represent MS injections.

    Attributes:
        chromatogram (np.ndarray): Chromatogram of the injection.
        peaks (dict[float, MS_Peak]): Peaks of the injection.
        scans (pd.DataFrame): Scans of the injection.
        detector (str): Detector used for the injection.
        analysis_settings (Analysis_Settings): Data method used for the injection.
        solvent_delay (float): Solvent delay applied to the injection.
    '''

    chromatogram: np.ndarray
    peaks: dict[float, MS_Peak]|None
    scans: pd.DataFrame
    detector: str
    analysis_settings: Analysis_Settings
    solvent_delay: float

    __slots__ = 'chromatogram', 'peaks', 'scans', 'detector', 'analysis_settings', 'solvent_delay'

    def __init__(self, metadata:dict|None, chromatogram:np.ndarray, peaks:dict|None, scans:pd.DataFrame, pos:bool=False):
        super().__init__(metadata, pos=pos)
        self.chromatogram = chromatogram
        self.peaks = peaks
        self.scans = scans
        self.detector = 'MS'
        self.analysis_settings = Analysis_Settings(chromatogram)
        self.solvent_delay = chromatogram[0][0]

    # def match_mz(self, mz:float) -> MS_Peak|list[MS_Peak]|None:
    #
    #     '''
    #     Returns the peak or a list of peaks whose mass spectra contain the given m/z. Returns None if no peak was found
    #     matching the criteria.
    #
    #     Args:
    #         mz (float): m/z to match.
    #
    #     Returns:
    #         MS_Peak|list[MS_Peak]|None: The peak or a list of peaks whose mass spectra contain the given m/z or None if
    #         no peak was found matching the criteria.
    #
    #     '''
    #
    #     candidates = []
    #     for rt, peak in self.peaks.items():
    #         if mz in peak.mass_spectrum['mz']:
    #             index = np.where(peak.mass_spectrum['mz'] == mz)[0]
    #             if peak.mass_spectrum['rel_intensity'][index][0] > 2 and mz > peak.mass_spectrum['mz'].max()*(2/3):
    #                 candidates.append(peak)
    #     if candidates:
    #         if len(candidates) > 1:
    #             print(f'Multiple peaks with m/z {mz} fitting the calculated isotope pattern were found for {self.sample_name}.')
    #             return candidates
    #         else:
    #             return candidates[0]
    #     return None

    def match_mol(self, smiles:str) -> MS_Peak|None:

        '''
        Returns the peak with the lowest isotope error for the m/z corresponding to the given molecule's parent
        peak. Returns None if no peak was found matching the criteria.

        Args:
            smiles (str): SMILES string of the analyte.

        Returns:
            MS_Peak|None: The peak with the lowest isotope error for the given molecule or None if no peak was
            found matching the criteria.
        '''

        mol = Chem.MolFromSmiles(smiles)
        mz = round(Descriptors.ExactMolWt(mol), 0)
        peak = self.__match_mz_mol(mz, smiles=smiles)
        if peak and peak.flag != 'standard':  # FBS for MS quantification the flag standard is not yet set at this point - however, might be used without standard anyways
        #if peak:
            analyte = Analyte(peak.rt, smiles=smiles)
            peak.analyte = analyte
        return peak

    def __match_mz_mol(self, mz:float, smiles:str) -> MS_Peak|None:

        '''
        Returns the peak with the lowest isotope error for a given m/z and molecule. Returns None if no peak was found
        matching the criteria.

        Args:
            smiles (str): SMILES string of the analyte.
            mz (float): m/z to match.

        Returns:
            MS_Peak|None: The peak with the lowest isotope error for the given m/z and molecule or None if no peak was
            found matching the criteria.
        '''

        candidates = {}
        for rt, peak in self.peaks.items():
            if mz in peak.mass_spectrum['mz']:
                index = np.where(peak.mass_spectrum['mz'] == mz)[0]
                # FBS peak.mass_spectrum['rel_intensity'][index][0] > 4 seems to be way to high: I checked for Standard Dodecane and it was only 1.64
                # Changed to 2.0
                # print(str(self.plate_pos) + " : " + str(rt)+ " : " + str(peak.mass_spectrum['rel_intensity'][index][0]))
                if peak.mass_spectrum['rel_intensity'][index][0] > 2.0 and mz > peak.mass_spectrum['mz'].max() * (
                        2 / 3):  # 2/3
                    # GC/MS 6 has relative low molecular/parent peaks, 2.0 should be fine
                    isotope_error = self.__isotope_check(smiles, peak, mz)
                    if isotope_error:
                        candidates[isotope_error] = peak
        if candidates:
            if len(candidates) > 1:
                # FBS Print information for manual inspection of the peaks with multiple peaks wich could be the analyte
                # Extract the retention times of the peaks in candidates
                retention_times = [peak.rt for peak in candidates.values()]
                print(f'{len(candidates)} peaks with m/z {mz} fitting the calculated isotope pattern were found for {self.sample_name:<20}. Retention times: {str(retention_times):<50}')
            peak = candidates[min(candidates)]  # selects the peak with the lowest isotope error
            # FBS: Due to Regioisomers and Isomers in the project with CSN and JLT, we select the peak with the highest area, within 5% of the isotope error
            # candidates_by_area = {}
            # for p in candidates.values():
            #     candidates_by_area[p.area] = p
            # peak = candidates_by_area[max(candidates_by_area)]

            peak.analyte = Analyte(peak.rt, smiles=smiles)
            return peak
        return None

    def pick_peaks(self, inplace: bool = True,  **kwargs: dict) -> None|dict[float, MS_Peak]:        # FBS Addition of a discovery mode

        '''
        Picks peaks from the injection's chromatogram.

        Args:
            inplace (bool): If True, the peaks are assigned to the injection's peaks attribute. Default is True.
            ms_quantification_mode (str): Methode of MS quantification ('height' or 'area'). Default is None.
            **kwargs: Keyword arguments for the peak picking.

        Returns:
            None|dict[float, MS_Peak]: The peaks of the injection if the inplace argument is False, None
            otherwise.
        '''

        self.analysis_settings.update(**kwargs)
        peaks = Peak_Detection_MS.pick_peaks(self.chromatogram, self.scans, self.analysis_settings)
        if inplace:
            self.peaks = peaks
        else:
            return peaks

    def __isotope_check(self, smiles:str, peak:MS_Peak, mz:float, error_margin:float = 0.055) -> float|None:

        '''
        Returns True if the isotope peak of a given m/z is present in the peak's mass spectrum and False otherwise.

        Args:
            smiles (str): SMILES string of the analyte.
            peak (MS_Peak): Peak to check.
            mz (float): m/z to check.

        Returns:
            float|None: The difference between the theoretical and the measured ratio of the isotope peak of a given m/z
            if the difference is smaller than 5% and None otherwise.
        '''

        if 'Cl' in smiles or 'Br' in smiles:
            diff = 2
        else:
            diff = 1
        return self.__isotopic_ratio_check(smiles, peak, mz, diff, error_margin)


    def __isotopic_ratio_check(self,smiles:str, peak:MS_Peak, mz:float, diff:int, error_margin:float = 0.055) -> float|None:

        '''
        Returns True if the ratio of the isotope peak of a given m/z is within 6% of the theoretical ratio and False if
        otherwise.

        Args:
            smiles (str): SMILES string of the analyte.
            peak (MS_Peak): Peak to check.
            mz (float): m/z to check.
            diff (int): Difference between the m/z of the isotope peak and the m/z of the peak to check.

        Returns:
            float|None: The difference between the theoretical and the measured ratio of the isotope peak of a given m/z
            if the difference is smaller than 5% and None otherwise.
        '''

        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        i = np.where(peak.mass_spectrum['mz'] == mz)[0]
        j = np.where(peak.mass_spectrum['mz'] == mz+diff)[0]
        if not i or not j:
            return False
        ratio = peak.mass_spectrum['rel_intensity'][j][0]/peak.mass_spectrum['rel_intensity'][i][0]
        mol_formula = self.__get_mol_formula_dict(mol)
        isotopic_dist = isotopic_variants(mol_formula, npeaks=3, charge=0)
        theo_ratio = isotopic_dist[diff].intensity/isotopic_dist[0].intensity
        if (theo_ratio - error_margin) < ratio < (theo_ratio + error_margin):   # FBS error_margin was 0.055
            return abs(theo_ratio - ratio)
        else:
            return None

    def __get_mol_formula_dict(self, mol) -> dict:

        '''
        Returns a dictionary with the molecule's elements as keys and the number of atoms of each element as values.

        Args:
            mol (rdkit Mol Object): Molecule to get the formula of.

        Returns:
            dict: Dictionary with the molecule's elements as keys and the number of atoms of each element as values.

        '''

        mol_formula = defaultdict(lambda : 0)
        for atom in mol.GetAtoms():
            mol_formula[atom.GetSymbol()] += 1
        return mol_formula

### FBS ###
### Fragmentation based on BDE ###

    def __calculate_isotope_mz(self, mol):
        '''
        Calculate the isotopic m/z for a given molecule.
        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object.
        Returns:
            float: m/z for the isotopic peak.
        '''
        isotope_masses = {}
        mol_formula = self.__get_mol_formula_dict(Chem.AddHs(mol))

        # if 'Cl' in mol_formula or 'Br' in mol_formula:
        #     diff = 2
        # else:
        #     diff = 1
        
        isotopic_dist = isotopic_variants(mol_formula, npeaks=3, charge=0)
        # isotope_mz = round(isotopic_dist[diff].mz, 0)
        for i in isotopic_dist:
            isotope_masses[round(i.mz, 0)] = Chem.MolToSmiles(mol)
                           
                           
        # fragment_masses = {mz: smiles}  # Initialize with parent peak
        # isotope_mzs = self.__calculate_isotope_mz(mol)  # Using calculated isotopic peak
        # fragment_masses[isotope_mz] = f"{smiles}" #(isotope)
        
        return isotope_masses

    ### Fragmentation based on BDE ###
    def match_mol_and_frags(self, smiles: str, fragments: dict = None,
                            fragment_match_percentage: float = 0.49) -> MS_Peak | None:  # CHANGE: Added fragment_match_percentage
        '''
        Returns the peak with the lowest isotope error for the m/z corresponding to the given molecule's parent
        peak and its fragments. Returns None if no peak was found matching the criteria.
        Args:
            smiles (str): SMILES string of the analyte.
            fragments (dict): Dictionary of fragments for the analyte. Default is None.
            fragment_match_percentage (float): Minimum percentage of fragments that need to be matched. Default is 0.5 (50%).
        Returns:
            MS_Peak|None: The peak with the lowest isotope error for the given molecule or None if no peak was
            found matching the criteria.
        '''
        mol = Chem.MolFromSmiles(smiles)
        mz = round(Descriptors.ExactMolWt(mol), 0)  # Parent peak

        # Generate the list of fragment masses (including the parent and isotope peak)
        fragment_masses = {mz: smiles}  # Initialize with parent peak
        isotope_mzs = self.__calculate_isotope_mz(mol)  # Using 3 calculated isotopic peaks
        fragment_masses = {**fragment_masses, **isotope_mzs}  # Add isotope peaks to the list
        # fragment_masses[isotope_mz] = f"{smiles}" #(isotope)
        if fragments:
            for fragment_smiles in fragments.keys():
                fragment_mol = Chem.MolFromSmiles(fragment_smiles)
                fragment_mz = round(Descriptors.ExactMolWt(fragment_mol), 0)
                # fragment_masses[fragment_mz] = fragment_smiles
                # FBS: Filter to exclude masses below 37 as they were not measured in the GC-MS method
                if fragment_mz >= 37:  # Filter to exclude masses below 37
                    fragment_masses[fragment_mz] = fragment_smiles

        peak = self.__match_mz_mol_and_frags(smiles = smiles, fragment_masses=fragment_masses, fragment_match_percentage = fragment_match_percentage)  # CHANGE
        return peak

    def __match_mz_mol_and_frags(self, smiles: str, fragment_masses: dict,
                                 fragment_match_percentage: float) -> MS_Peak | None:  # CHANGE
        '''
        Returns the peak with the lowest isotope error for given m/z values and molecule. Returns None if no peak was found
        matching the criteria.
        Args:
            fragment_masses (dict): Dictionary of m/z values and corresponding SMILES strings to match.
            fragment_match_percentage (float): Minimum percentage of fragments that need to be matched.
        Returns:
            MS_Peak|None: The peak with the lowest isotope error or None if no peak was found matching the criteria.
        '''
        candidates = []

        for rt, peak in self.peaks.items():
            matched_fragments = []
            for mz in fragment_masses.keys():
                if mz in peak.mass_spectrum['mz']:
                    index = np.where(peak.mass_spectrum['mz'] == mz)[0]
                    matched_fragments.append(mz)    # FBS indent deleted

            if len(matched_fragments) / len(fragment_masses) >= fragment_match_percentage:
                # TODO: Do we still need to check for isotopic error?
                isotope_error = self.__isotope_check(smiles= fragment_masses[matched_fragments[0]], peak = peak,
                                                     mz = matched_fragments[0], error_margin = 0.195)
                candidates.append((peak, len(matched_fragments), isotope_error))

                # if isotope_error:
                #     candidates.append((peak, len(matched_fragments), isotope_error))
                # 
        if candidates:
            # # First, sort by the number of matched fragments in descending order
            # candidates.sort(key=lambda x: x[1], reverse=True)
            # # Among the candidates with the highest number of matched fragments, select the one with the lowest isotope error
            # best_candidate = min(candidates[:len([x for x in candidates if x[1] == candidates[0][1]])], key=lambda x: (x[2] is not None, x[2] if x[2] is not None else float('inf')))
            # peak = best_candidate[0]

            # : Due to Regioisomers and Isomers in the project with CSN and JLT, we select the peak with the highest area, within 5% of the isotope error
            candidates_by_area = {}
            for i in candidates:
                p = i[0]
                candidates_by_area[p.area] = p
            peak = candidates_by_area[max(candidates_by_area)]



            # Fix assignment of analyte SMILES
            peak.analyte = Analyte(peak.rt, smiles=smiles)
            return peak

        return None
            
            
        #     if len(candidates) > 1:
        #         # FBS Print information for manual inspection of the peaks with multiple peaks wich could be the analyte
        #         # Extract the retention times of the peaks in candidates
        #         retention_times = [peak.rt for peak in candidates.values()]
        #         print(
        #             f'{len(candidates)} peaks with m/z {mz} fitting the calculated isotope pattern were found for {self.sample_name:<20}. Retention times: {str(retention_times):<50}')
        #     peak = candidates[min(candidates)]  # selects the peak with the lowest isotope error
        #     # FBS: Due to Regioisomers and Isomers in the project with CSN and JLT, we select the peak with the highest area, within 5% of the isotope error
        #     # candidates_by_area = {}
        #     # for p in candidates.values():
        #     #     candidates_by_area[p.area] = p
        #     # peak = candidates_by_area[max(candidates_by_area)]
        # 
        #     peak.analyte = Analyte(peak.rt, smiles=smiles)
        #     return peak
        # return None
