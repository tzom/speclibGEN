from pyteomics import auxiliary as aux
from pyteomics import mgf
import pandas as pd
from tqdm import tqdm

ONLY_TAKE_FIRST_N = None

import argparse
parser = argparse.ArgumentParser(description='Generate spectra library from msfragger pin file and corresponding mgf file.')
parser.add_argument('--pin', help='msfragger pin file', required=True)
parser.add_argument('--mgf_in', help='input mgf file', required=True)
parser.add_argument('--mgf_out', help='output mgf file', required=True)

#PIN_FILE="/hpi/fs00/home/tom.altenburg/scratch/tmp/yHydra_benchmark/example/qe2_03132014_1WT-1.pin"
#MGF_FILE="/hpi/fs00/home/tom.altenburg/scratch/tmp/yHydra_benchmark/example/qe2_03132014_1WT-1_uncalibrated.mgf"

PIN_FILE=parser.parse_args().pin
MGF_IN=parser.parse_args().mgf_in
MGF_OUT=parser.parse_args().mgf_out


msfragger_pin_cols = ['SpecId', 'Label', 'ScanNr', 'ExpMass', 'retentiontime', 'rank', 'mass_diff_score', 'log10_evalue', 'hyperscore', 'delta_hyperscore', 'matched_ion_num',
                      'peptide_length', 'ntt', 'nmc', 'charge_1', 'charge_2', 'charge_3', 'charge_4', 'charge_5', 'charge_6', 'charge_7', 'Peptide', 'Proteins']


def _read_pin(pin_file):
    pin = pd.read_csv(pin_file, sep='\t', nrows=ONLY_TAKE_FIRST_N, usecols=msfragger_pin_cols)
    return pin

search_results = _read_pin(PIN_FILE)

search_results['is_decoy'] = search_results['Proteins'].str.contains('rev_')
search_results['score'] = -search_results['hyperscore']

search_results_filtered = aux.filter(search_results, key='score', is_decoy='is_decoy', fdr=0.01)
search_results_filtered = search_results_filtered[~search_results_filtered.is_decoy]
#search_results_filtered = search_results_filtered.drop_duplicates(subset=['Peptide'])

print("Identified %s PSMs"%len(search_results_filtered))

search_results_filtered["title"] = search_results_filtered.SpecId.str.split('_').str[:-1].str.join('_')

IDENTIFIED_SPECTRA = set(search_results_filtered["title"].values)

SPECTRA_TO_WRITE = []

print("collecting spectra...")
with mgf.read(MGF_IN) as mgf_reader:
    for i,entry in enumerate(tqdm(mgf_reader)):   
        title = entry['params']['title']     
        if title in IDENTIFIED_SPECTRA:
            row = search_results_filtered.index[search_results_filtered["title"] == title].tolist()
            assert len(row) == 1
            PSM = search_results_filtered.loc[row]
            peptide = PSM.Peptide.str.slice(2,-2)
            entry['params']['SEQ'] = peptide.values[0]
            SPECTRA_TO_WRITE.append(entry)
print("done.")

print("writing mgf...")
with mgf.write(SPECTRA_TO_WRITE, MGF_OUT) as mgf_writer:
    pass
print("done.")