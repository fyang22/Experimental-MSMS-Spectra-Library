# Import libraries
import pandas as pd
import glob
import os
from pyteomics import mgf

# Read mgf files
input_mgf_files = glob.glob("Data/pos/MGF_10042025/*/*.mgf")
print(input_mgf_files)
print(len(input_mgf_files))

# Check the sanity of mgf files
for filename in input_mgf_files:
    file= mgf.MGF(source=filename , use_header=True, convert_arrays=2, read_charges=True, read_ions=False, dtype=None, encoding=None)
    print(file)

# Read metadata
df = pd.read_csv(os.path.join("Data", "pos","All_lib_10042025.csv"),encoding='cp1252')

# Combine as a mgf file
with open('Lib_Positive_10042025.mgf','w', encoding="utf-8") as out:
    for filename in input_mgf_files:
        file = os.path.basename(filename).split('.')[0]
        with mgf.read(filename,encoding='latin-1') as spectra:
            spec = next(spectra)
            for index, row in df.iterrows():
                file = row['compound_id']
                compound_name = row['name']
                Precursor_mz = row['mz']
                rt = row['rtime']

            #title = spec['params']['title']
                spec['params']['NAME'] = compound_name
                #spec['params']['CHARGE'] = -1
                spec['params']['Spectrum_type'] = 'MS2'
                spec['params']['PRECURSORTYPE'] = 'M+H'
                spec['params']['PrecursorMZ'] = Precursor_mz
                #spec['params']['C'] = rt
                spec['params']['Comment'] = rt
                #spec['params']['IONMODE'] = 'Negative'

                mgf.write([spec],out)

# Combine to msp
def multiple_mgf_to_msp(mgf_files, msp_file, df):
    # creat an empty writable msp file
    with open(msp_file, 'w', encoding="utf-8") as msp:  
        for mgf_file in mgf_files: # parse mgf files in the directory
            file_base_name = os.path.basename(mgf_file).split('.')[0]
            # match the mgf file name with metadate using 'compound_id'
            matching_row = df[df['compound_id'] == file_base_name]
            

            if not matching_row.empty:
                compound_name = matching_row.iloc[0]['name']
                Precursor_mz = matching_row.iloc[0]['mz']
                
                
                # Read the mgf file and write to MSP
                with mgf.read(mgf_file, encoding='latin-1') as spectra:
                    for spectrum in spectra:
                        msp.write(f"Name: {compound_name}\n")
                        msp.write(f"PrecursorMZ: {Precursor_mz}\n")
                        msp.write(f"PRECURSORTYPE: M+H\n")
                        #msp.write(f"PRECURSORTYPE: M+H\n")
                        msp.write(f"Spectrum_type: MS2\n")
                        msp.write(f"Num peaks: {len(spectrum['m/z array'])}\n")
                        
                        # Write Mass spectral
                        for mz, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):
                            msp.write(f"{mz} {intensity}\n")
                        
                        # Separate spectra with a blank line
                        msp.write("\n")
            else:
                print(f"Warning: No matching compound data found for file {mgf_file}")

multiple_mgf_to_msp(input_mgf_files, 'Lib_Positive_10042025.msp',df)