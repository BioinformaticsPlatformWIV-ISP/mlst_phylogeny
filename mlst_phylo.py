#!/usr/bin/env python
import argparse
import itertools
import logging
from pathlib import Path
from typing import Optional
from collections.abc import Sequence

import pandas as pd


class MLSTPhylogeny:
    """
    Class to generate phylogenies from sequence typing outputs.

    The expected input is a list of TSV with the following columns:
    Locus (string), Allele (string / int), % identity (float), HSP/Locus length (string).

    For example:
    Locus	Allele	% Identity	HSP/Locus length	Type
    SC0831	1	100.00	129/129	DNA
    SEN0401	10	100.00	978/978	DNA
    SNSL254_A1382	1	100.00	117/117	DNA
    SPAB_04503	4	100.00	108/108	DNA
    STM0834	1	100.00	1581/1581	DNA

    Missing alleles should be indicated with allele '-', multi-hits as '?'.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MLSTPhylogeny.parse_args(args)

    @staticmethod
    def parse_args(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: (optional) arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('--input-tsv', action='append', type=Path)
        parser.add_argument('--output-matrix', type=Path, help='Filtered allele matrix (TSV)', required=True)
        parser.add_argument(
            '--output-dists', type=Path, help='Pairwise distance matrix (TSV)', required=True)
        parser.add_argument(
            '--min-perc-loci', type=int, default=90,
            help='Minimum percentage of loci that should be present in a dataset')
        parser.add_argument(
            '--min-perc-samples', type=int, default=90,
            help='Minimum percentage of datasets where loci should be present')
        return parser.parse_args(args)

    @staticmethod
    def is_perfect(record: pd.Series) -> bool:
        """
        Determines if the given hit is perfect.
        :param record: Input record
        :return: True if perfect, False otherwise
        """
        if str(record['Allele']) in ('-', '?'):
            return False
        try:
            identity = float(record['% Identity'])
        except ValueError:
            return False
        if identity != 100.0:
            return False
        len_hsp, len_locus = record['HSP/Locus length'].split('/')
        if not len_hsp == len_locus:
            return False
        return True

    @staticmethod
    def calculate_distance(alleles_a: pd.Series, alleles_b: pd.Series) -> int:
        """
        Calculates the number of different alleles.
        :param alleles_a: Alleles a
        :param alleles_b: Alleles b
        :return: Number of allelic differences
        """
        total_dist = 0
        for allele_a, allele_b in zip(alleles_a, alleles_b):
            if allele_a == '-' and allele_b == '-':
                continue
            else:
                total_dist += 0 if allele_a == allele_b else 1
        return total_dist

    def parse_tsv_typing_list(self, tsv_in: list[Path]) -> pd.DataFrame:
        """
        Parses a list of tabular input files.
        :param tsv_in: List of input TSV files
        :return: DataFrame with detected alleles
        """
        sample_names = []
        allele_data = []
        for path_tsv in tsv_in:
            logging.debug(f'Parsing file: {path_tsv}')
            isolate_name = path_tsv.name
            data_tsv = pd.read_table(path_tsv)
            data_tsv['is_perfect_hit'] = data_tsv.apply(lambda x: MLSTPhylogeny.is_perfect(x), axis=1)
            alleles_parsed = {r['Locus']: r['Allele'] if r['is_perfect_hit'] else '-' for _, r in data_tsv.iterrows()}
            logging.debug(f"{sum(v != '-' for _, v in alleles_parsed.items()):,}/{len(alleles_parsed):,} perfect hits")
            allele_data.append(alleles_parsed)
            sample_names.append(isolate_name)
        return pd.DataFrame(allele_data, index=sample_names, dtype=str)

    def calculate_distance_matrix(self, allele_data_filtered: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates the pairwise distance matrix
        :param allele_data_filtered: Filtered allele matrix
        :return: Distance matrix
        """
        # Calculate pair-wise distances
        distance_by_dataset_pair = {}
        for dataset_a, dataset_b in itertools.combinations(allele_data_filtered.index, r=2):
            key = tuple(sorted([dataset_a, dataset_b]))
            dist = MLSTPhylogeny.calculate_distance(allele_data_filtered.loc[dataset_a], allele_data_filtered.loc[dataset_b])
            distance_by_dataset_pair[key] = dist

        # Create data frame with pairwise distances
        records_out = []
        for dataset_a in allele_data_filtered.index:
            records_out.append({
                dataset_b: distance_by_dataset_pair.get(tuple(sorted([dataset_a, dataset_b])), 0) for
                dataset_b in allele_data_filtered.index
            })
        return pd.DataFrame(records_out, index=allele_data_filtered.index)

    def filter_allele_matrix(self, allele_data: pd.DataFrame) -> pd.DataFrame:
        """
        Filters the allele matrix by removing:
        - Datasets with less than x% of loci detected
        - Loci present in less than x% of datasets
        :param allele_data: Allele data
        :return: Filtered allele data, loci cutoff, datasets cutoff
        """
        # Filter allele matrix (nb. of loci detected per dataset)
        nb_loci_detected = allele_data.apply(lambda x: len(x) - list(x).count('-'), axis=1)
        cutoff_loci = int(self._args.min_perc_loci * len(allele_data.columns) / 100)
        logging.info(f"Removing datasets with < {cutoff_loci} ({self._args.min_perc_loci}%) loci detected")
        allele_data_filt = allele_data[nb_loci_detected > cutoff_loci]
        logging.info(f"{len(allele_data_filt)} datasets passed filtering")

        # Filter allele matrix (loci detected in nb. of datasets)
        cutoff_datasets = int(self._args.min_perc_samples * len(allele_data_filt) / 100)
        logging.info(f"Removing loci detected in < {cutoff_datasets} ({self._args.min_perc_samples}%) samples")
        locus_present_in_datasets = allele_data_filt.apply(lambda x: len(x) - list(x).count('-'))
        allele_data_filt = allele_data_filt.iloc[:, list(locus_present_in_datasets > cutoff_datasets)]
        logging.info(f"{len(allele_data_filt.columns)} loci passed filtering")
        return allele_data_filt

    def run(self) -> None:
        """
        The main method to construct the MLST phylogeny.
        :return: None
        """
        # Parse the input files
        df_alleles = self.parse_tsv_typing_list(self._args.input_tsv)

        # Filter the allele matrix
        df_alleles_filt = self.filter_allele_matrix(df_alleles)
        df_alleles_filt.to_csv(self._args.output_matrix, index_label='ID', sep='\t')
        logging.info(f'Allele matrix exported to: {self._args.output_matrix}')

        # Calculate the distance matrix
        df_dists = self.calculate_distance_matrix(df_alleles_filt)

        # Check the distance matrix
        if all(max(values) == 0 for _, values in df_dists.iterrows()):
            raise ValueError('Empty distance matrix')

        df_dists.to_csv(self._args.output_dists, index_label='ID', sep='\t')
        logging.info(f'Distance matrix exported to: {self._args.output_dists}')

        # Log the command for GrapeTree
        logging.info(f'You can construct a phylogeny using: grapetree --profile {self._args.output_matrix.name} --method MSTreeV2')


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(module)15s - %(levelname)7s - %(message)s')
    mlst_phylo = MLSTPhylogeny()
    mlst_phylo.run()
