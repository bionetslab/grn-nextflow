#!/usr/bin/env python
"""Usage: gene_set_enrichment.py -i FILE -o FOLDER

Wrapper script taking as input the network diff file (as produced by boostdiff).

-h --help               show this
-i --inputfile FILE    specify input file
-o --output FOLDER      specify output folder
--verbose    print more text

"""
import pandas as pd
import argparse
import gseapy as gp
import os

def argument_parser():

    parser = argparse.ArgumentParser(description="Gene set enrichment argument parser")

    parser.add_argument('-o', '--output', type=str, help='Specifies the output folder')
    parser.add_argument('-i', '--input_file', type=str, help='Specifies the input file')

    args = parser.parse_args()
    return args


def gene_set_enrichment(input_file, output_dir):

    # Read the TSV file using pandas
    df = pd.read_csv(input_file, sep='\t')

    # Removing singletons
    value_counts = df['target'].value_counts()
    filtered_values = value_counts[value_counts > 1]
    target_names_to_keep = filtered_values.index.tolist()
    mask = df['target'].isin(target_names_to_keep)
    df = df[mask]
    targets = df['target'].unique().tolist()

    data_bases = [('KEGG', ['KEGG_2019_Mouse']),
                  ('MSigDB', ['MSigDB_Hallmark_2020']),
                  ('GO', ['GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021', 'SynGO_2022'])
                 ]

    res = pd.DataFrame()

    for (name, data_base) in data_bases:

        enrichment = gp.enrichr(gene_list=targets,
                            gene_sets=data_base,
                            organism='Mouse',
                            outdir=None,
                            )

        res = pd.concat([res, enrichment.results], axis=0)
        # plotting does not really work well with nf here (because there can be no enriched terms to plot -> no guarenteed plot)
        # ax = gp.barplot(enrichment.results,
        #         column="Adjusted P-value",
        #         group='Gene_set',
        #         size=10,
        #         top_term=5,
        #         figsize=(3,5),
        #         color=['darkred', 'darkblue', 'darkgreen', 'orange'],
        #         ofname=ofname
        #     )
    os.makedirs(output_dir, exist_ok=True)
    res.to_csv(os.path.join(output_dir, 'gse_analysis.tsv'), sep="\t")


if __name__ == '__main__':
    args = argument_parser()
    print(args.input_file)
    print(args.output)
    gene_set_enrichment(input_file=args.input_file, output_dir=args.output)
