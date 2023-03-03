#!/usr/bin/env python3

import argparse
import pandas

DESCRIPTION = "convert_ids2table - identifier mapping conversion"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('infile', help="input file")
    parser.add_argument('outfile', help="output file")
    args = parser.parse_args()
    ens2symbol = {}
    ens2entrez = {}
    df = pandas.read_csv(args.infile, sep='\t', header=0)
    entrez_lines = df[df['Source'] == 'Entrez Gene ID']
    symbol_lines = df[df['Source'] == 'Gene Name']

    for index, row in entrez_lines.iterrows():
        ens, entrez = row['Preferred_Name'], row['Name']
        ens2entrez[ens] = entrez

    for index, row in symbol_lines.iterrows():
        ens, symbol = row['Preferred_Name'], row['Name']
        ens2symbol[ens] = symbol
    #print(ens2symbol)
    #print(ens2entrez)

    with open(args.outfile, 'w') as outfile:
        outfile.write("preferred\tensembl\tentrez\n")
        for ens, symbol in ens2symbol.items():
            try:
                entrez = ens2entrez[ens]
            except KeyError:
                entrez = ens
            row = "%s\t%s\t%s\n" % (ens, symbol, str(entrez))
            outfile.write(row)
