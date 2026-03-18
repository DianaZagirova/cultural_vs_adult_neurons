#! /usr/bin/env python
import pandas as pd
import numpy as np
import scipy
import cooler
import multiprocess
from multiprocess import Pool
from functools import partial
from scipy import stats
import warnings
import  sys
import os
import argparse
import matplotlib.pyplot as plt
import gc
import psutil
from contact_hunter.utils import (
     preproc_data,norm,create_profiles,pvalue_calculation,sign_contact,postprocess_prepare,
     qvalue_calculation,bad_distance,plot_heatmap,significant_contacts_average_heatmap,
     post_processing_for_sign_contacts,create_final_table,cpu_number,get_contacts_cli)

def main():
    parser = argparse.ArgumentParser(description='The program identifies significant contacts for provided genomic regions across predefined distance. The method was borrowed from the paper  DOI:10.1038/nature19847 (Won et al 2016) with a small difference - significant contacts are detected for normalized Hi-C maps (obs/exp). Parameters of the tool were tuned to generate list of significant contacts for regions containing SNP in human genome. SNPs to be tested are credible SNP (after CAVIAR/FINEMAP), background SNPs are all SNPs obtained by GWAS. The tool is suitable to obtain significant contacts for TSSs too, but another background should be provided. Performance on other species or features (e.g. midpoints of ATAC-seq or ChIP-seq peaks) should be tested at first.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("COOL_PATH", help="path to cool-file (or mcool file in format 'file.mcool::resolutions/res' where res is the resolution selected)")
    parser.add_argument("BACKGROUND_LOCI", help="path to tab-delimited file with coordinates of background regions (points), should not contain header, 2 filds are expected: chrom, start")
    parser.add_argument("LOCI_TO_TEST", help="path to a file with loci of interes, same format as for background regions file")
    parser.add_argument("RES", type=int,help="resolution, dependes on the data quality (e.g. 5000-10000 for mammalian Hi-C maps of good quality, 20000-25000 for worse quality);we do not recommend use human Hi-C data with resolution higher than 5000 because it can be very computationally demanding")
    parser.add_argument('DIST',type=int,help="distance from loci of interes (in bp), constraining the field of significant contacts search (e.g. 5000000 for mammalian genome)")
    parser.add_argument('RESULTS_FILE', help="path to output file with significant contacts")
    parser.add_argument("--chr",default='all',type=str,help="comma-separated list of chromosomes to be included, e.g. chr1,chr2; if not specified, inspects all the chromosomes")
    parser.add_argument("--fdr",default=0.01,type=float,help="false discovery rate")
    parser.add_argument("--heatmap_generate",action="store_true",default=False, help='plot an average heatmap for significant contacts and save it to the directory with the tested loci file')
    args = parser.parse_args()
    config = vars(args)
    if "chr" in config.keys():
        config['chr']=[s.strip() for s in config["chr"].split(",")]
    print('Hello. The calculations may take a while.')
    df=get_contacts_cli(config["COOL_PATH"],config["BACKGROUND_LOCI"],config["LOCI_TO_TEST"],config["RES"],config['DIST'],config['chr'],config['fdr'],config['heatmap_generate'])
    df[0].to_csv(config['RESULTS_FILE'],index=None,sep='\t')
    if config['heatmap_generate']:
        path=os.path.abspath(config["LOCI_TO_TEST"])
        path=path.split('/')
        path='/'.join(path[:-1])
        plot_heatmap(df[1],"%s/heatmap_for_significant_contacts.png"%path)

