
import numpy as np
from scipy import sparse
from anndata import AnnData
import pandas as pd
import argparse
import sys

def main(fragments):

    data, row, column = [],[],[]

    barcodes = {}
    peaks = {}
    
    for fragment in fragments:

        fragment = fragment.strip().split('\t')
        barcode = fragment[3]
        peak = tuple(fragment[4:7])

        if not peak[0] == '.':

            if not barcode in barcodes:
                barcodes[barcode] = len(barcodes)
            
            if not peak in peaks:
                peaks[peak] = len(peaks)

            data.append(1.0)
            row.append(barcodes[barcode])
            column.append(peaks[peak])

    counts = sparse.coo_matrix(
        (data, (row, column)),
        shape = (len(barcodes), len(peaks))
    ).tocsc()

    obs = pd.DataFrame([(str(i), j) for j,i in barcodes.items()], columns = ['index','barcode']).set_index('index')
    var = pd.DataFrame([(str(i),*j) for j,i in peaks.items()], columns = ['index','chr','start','end']).set_index('index')

    return AnnData(X = counts, obs = obs, var = var)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fragments', type = argparse.FileType('r'), required=True)
    parser.add_argument('-o','--output',type = str, required=True)

    args = parser.parse_args()

    count_matrix = main(args.fragments)

    print('Created count matrix:\n' + repr(count_matrix), file = sys.stderr)
    
    #sparse.save_npz(args.output, count_matrix)
    count_matrix.write(args.output)


