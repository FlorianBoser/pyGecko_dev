# pyLurch: pyGecko MS Quantification
<img src="docs/pyGecko_icon.png" alt="pyGecko_Logo" width="300" height="300"/>

> pyGecko an open-source Python library for the parsing, processing and analysis of GC-MS and GC-FID raw data.

## Acknowledgement

This modified version of pyGecko builds upon the original [pyGecko](https://github.com/FelixKatz77/pyGecko) library and [publication](https://chemrxiv.org/engage/chemrxiv/article-details/66adfc465101a2ffa8001761), created for the parsing, processing, and analysis of GC-MS and GC-FID raw data.
This fork pyLurch implements functionalities for a qualitative yield estimation soley based on MS data.


## About

With increasing amounts of analytical and metadata generated in HTE, data processing and analysis quickly become a 
workflow's limiting step if conducted manually. The automated processing of analytical data opens up time for chemists 
to focus on relevant outcomes, enables the standardized storage of reaction data, and facilitates the integration of 
analytical methods into closed-loop systems. Herein we present pyGecko, an open-source Python library for the parsing,
processing and analysis of GC-MS and GC-FID raw data. pyGecko offers a variety of analysis tools for the automated or 
semi-automated handling of GC measurements and sequences. This includes the interpretation of measurements in the context 
of the experiment, the automatic identification of internal standards and compound identifications based on retention 
times, the mass of a molecular ion or fragment and spectral comparison. Quantification relative to an internal standard 
can be performed for GC-FID measurements. Results of an analysis as well as chromatograms and spectra can be visualized 
and reported in standardized formats like the Open Reaction Database (ORD) schema. pyGecko is designed to be easily 
integrated into automated workflows and can be used as a stand-alone tool or as a python library.

## Installation

> [!IMPORTANT]
> To read vendor files you need to install the msConvert tool from ProteoWizard. You can download it from [here](http://proteowizard.sourceforge.net/download.html).
> You need to specify the path to the msConvert.exe before the first run of pyGecko.

pyGecko can be installed via pip:

```bash 
git clone https://github.com/FlorianBoser/pyLurch.git
cd pyLurch
pip install ./
```
Afterward the path to the msConvert.exe needs to be specified. This can be done by running the following command:

```bash
cd pygecko
python __init__.py
```
This will prompt you to specify the path to the msConvert.exe file:

```bash
Please provide the path to the msConvert executable or specify it in the config.ini:
```
After that pyLurch is ready to use.


## Documentation
The documentation for pyGecko can be found [here](https://pygecko.readthedocs.io/en/latest/).

## Usage
For non-automated workflows pyGecko is best used with jupyter notebooks. The notebooks folder of the repository contains
examples for the usage of pyGecko for the quantitative analysis of reaction outcomes and spectral matching. The Python 
scripts used to perform the data processing for the publication can be found in the examples folder. 

## Supported File Formats
pyGecko supports the following file formats:

| GC-MS         | GC-FID    |
|---------------|-----------|
| .mzML         | .xy       |
| .mzXML        | .CSV      |
| .D (Agilent)  ||
| .RAW (Thermo) ||

> [!NOTE]
> To achieve the best performance, we recommend using the .mzML file format for GC-MS data.

## Citation & Reference Paper

[Calibration-Free Quantification and Automated Data Analysis for High-Throughput Reaction Screening](https://chemrxiv.org/engage/chemrxiv/article-details/66adfc465101a2ffa8001761)
Felix Katzenburg, Florian Boser, Felix Schäfer, Philipp Pflüger, Calibration-Free Quantification and Automated Data Analysis for High-Throughput Reaction Screening, ChemRxiv, 2024; doi:10.26434/chemrxiv-2024-1ctkh
