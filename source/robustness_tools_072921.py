#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:14:51 2021

@author: bryanandrews1
"""

import sys
import numpy
import random
import os

translation_table = {    
    'ACC': "T", 'ATG': "M", 'ACA': "T", 
    'ACG': "T", 'ATC': "I", 'AAC': "N", 
    'ATA': "I", 'AGG': "R", 'CCT': "P", 
    'CTC': "L", 'AGC': "S", 'AAG': "K", 
    'AGA': "R", 'CAT': "H", 'AAT': "N", 
    'ATT': "I", 'CTG': "L", 'CTA': "L", 
    'ACT': "T", 'CAC': "H", 'AAA': "K", 
    'CCG': "P", 'AGT': "S", 'CCA': "P", 
    'CAA': "Q", 'CCC': "P", 'TAT': "Y", 
    'GGT': "G", 'TGT': "C", 'CGA': "R", 
    'CAG': "Q", 'CGC': "R", 'GAT': "D", 
    'CGG': "R", 'CTT': "L", 'TGC': "C", 
    'GGG': "G", 'TAG': "*", 'GGA': "G", 
    'TAA': "*", 'GGC': "G", 'TAC': "Y", 
    'GAG': "E", 'TCG': "S", 'TTT': "F", 
    'GAC': "D", 'CGT': "R", 'GAA': "E", 
    'TCA': "S", 'GCA': "A", 'GTA': "V", 
    'GCC': "A", 'GTC': "V", 'GCG': "A", 
    'GTG': "V", 'TTC': "F", 'GTT': "V", 
    'GCT': "A", 'TTA': "L", 'TGA': "*", 
    'TTG': "L", 'TCC': "S", 'TGG': "W", 
    'TCT': "S"}
    
codon_neighbors = {    
    'ACC': ['CCC', 'GCC', 'TCC', 'AAC', 'AGC', 'ATC', 'ACA', 'ACG', 'ACT'], 
    'ATG': ['CTG', 'GTG', 'TTG', 'ACG', 'AGG', 'AAG', 'ATC', 'ATA', 'ATT'], 
    'ACA': ['CCA', 'GCA', 'TCA', 'AAA', 'AGA', 'ATA', 'ACC', 'ACG', 'ACT'], 
    'ACG': ['CCG', 'GCG', 'TCG', 'AAG', 'AGG', 'ATG', 'ACC', 'ACA', 'ACT'], 
    'ATC': ['CTC', 'GTC', 'TTC', 'ACC', 'AGC', 'AAC', 'ATA', 'ATG', 'ATT'], 
    'AAC': ['CAC', 'GAC', 'TAC', 'ACC', 'AGC', 'ATC', 'AAA', 'AAG', 'AAT'], 
    'ATA': ['CTA', 'GTA', 'TTA', 'ACA', 'AGA', 'AAA', 'ATC', 'ATG', 'ATT'], 
    'AGG': ['CGG', 'GGG', 'TGG', 'ACG', 'AAG', 'ATG', 'AGC', 'AGA', 'AGT'], 
    'CCT': ['ACT', 'GCT', 'TCT', 'CAT', 'CGT', 'CTT', 'CCC', 'CCG', 'CCA'], 
    'ACT': ['CCT', 'GCT', 'TCT', 'AAT', 'AGT', 'ATT', 'ACC', 'ACG', 'ACA'], 
    'AGC': ['CGC', 'GGC', 'TGC', 'ACC', 'AAC', 'ATC', 'AGA', 'AGG', 'AGT'], 
    'AAG': ['CAG', 'GAG', 'TAG', 'ACG', 'AGG', 'ATG', 'AAC', 'AAA', 'AAT'], 
    'AGA': ['CGA', 'GGA', 'TGA', 'ACA', 'AAA', 'ATA', 'AGC', 'AGG', 'AGT'], 
    'CAT': ['AAT', 'GAT', 'TAT', 'CCT', 'CGT', 'CTT', 'CAC', 'CAG', 'CAA'], 
    'AAT': ['CAT', 'GAT', 'TAT', 'ACT', 'AGT', 'ATT', 'AAC', 'AAG', 'AAA'], 
    'ATT': ['CTT', 'GTT', 'TTT', 'ACT', 'AGT', 'AAT', 'ATC', 'ATG', 'ATA'], 
    'CTG': ['ATG', 'GTG', 'TTG', 'CCG', 'CGG', 'CAG', 'CTC', 'CTA', 'CTT'], 
    'CTA': ['ATA', 'GTA', 'TTA', 'CCA', 'CGA', 'CAA', 'CTC', 'CTG', 'CTT'], 
    'CTC': ['ATC', 'GTC', 'TTC', 'CCC', 'CGC', 'CAC', 'CTA', 'CTG', 'CTT'], 
    'CAC': ['AAC', 'GAC', 'TAC', 'CCC', 'CGC', 'CTC', 'CAA', 'CAG', 'CAT'], 
    'AAA': ['CAA', 'GAA', 'TAA', 'ACA', 'AGA', 'ATA', 'AAC', 'AAG', 'AAT'], 
    'CCG': ['ACG', 'GCG', 'TCG', 'CAG', 'CGG', 'CTG', 'CCC', 'CCA', 'CCT'], 
    'AGT': ['CGT', 'GGT', 'TGT', 'ACT', 'AAT', 'ATT', 'AGC', 'AGG', 'AGA'], 
    'CCA': ['ACA', 'GCA', 'TCA', 'CAA', 'CGA', 'CTA', 'CCC', 'CCG', 'CCT'], 
    'CAA': ['AAA', 'GAA', 'TAA', 'CCA', 'CGA', 'CTA', 'CAC', 'CAG', 'CAT'], 
    'CCC': ['ACC', 'GCC', 'TCC', 'CAC', 'CGC', 'CTC', 'CCA', 'CCG', 'CCT'], 
    'TAT': ['CAT', 'GAT', 'AAT', 'TCT', 'TGT', 'TTT', 'TAC', 'TAG', 'TAA'], 
    'GGT': ['CGT', 'AGT', 'TGT', 'GCT', 'GAT', 'GTT', 'GGC', 'GGG', 'GGA'], 
    'TGT': ['CGT', 'GGT', 'AGT', 'TCT', 'TAT', 'TTT', 'TGC', 'TGG', 'TGA'], 
    'CGA': ['AGA', 'GGA', 'TGA', 'CCA', 'CAA', 'CTA', 'CGC', 'CGG', 'CGT'], 
    'CAG': ['AAG', 'GAG', 'TAG', 'CCG', 'CGG', 'CTG', 'CAC', 'CAA', 'CAT'], 
    'CGC': ['AGC', 'GGC', 'TGC', 'CCC', 'CAC', 'CTC', 'CGA', 'CGG', 'CGT'], 
    'GAT': ['CAT', 'AAT', 'TAT', 'GCT', 'GGT', 'GTT', 'GAC', 'GAG', 'GAA'], 
    'CGG': ['AGG', 'GGG', 'TGG', 'CCG', 'CAG', 'CTG', 'CGC', 'CGA', 'CGT'], 
    'CTT': ['ATT', 'GTT', 'TTT', 'CCT', 'CGT', 'CAT', 'CTC', 'CTG', 'CTA'], 
    'TGC': ['CGC', 'GGC', 'AGC', 'TCC', 'TAC', 'TTC', 'TGA', 'TGG', 'TGT'], 
    'GGG': ['CGG', 'AGG', 'TGG', 'GCG', 'GAG', 'GTG', 'GGC', 'GGA', 'GGT'], 
    'TAG': ['CAG', 'GAG', 'AAG', 'TCG', 'TGG', 'TTG', 'TAC', 'TAA', 'TAT'], 
    'GGA': ['CGA', 'AGA', 'TGA', 'GCA', 'GAA', 'GTA', 'GGC', 'GGG', 'GGT'], 
    'TAA': ['CAA', 'GAA', 'AAA', 'TCA', 'TGA', 'TTA', 'TAC', 'TAG', 'TAT'], 
    'GGC': ['CGC', 'AGC', 'TGC', 'GCC', 'GAC', 'GTC', 'GGA', 'GGG', 'GGT'], 
    'TAC': ['CAC', 'GAC', 'AAC', 'TCC', 'TGC', 'TTC', 'TAA', 'TAG', 'TAT'], 
    'GAG': ['CAG', 'AAG', 'TAG', 'GCG', 'GGG', 'GTG', 'GAC', 'GAA', 'GAT'], 
    'TCG': ['CCG', 'GCG', 'ACG', 'TAG', 'TGG', 'TTG', 'TCC', 'TCA', 'TCT'], 
    'TTT': ['CTT', 'GTT', 'ATT', 'TCT', 'TGT', 'TAT', 'TTC', 'TTG', 'TTA'], 
    'GAC': ['CAC', 'AAC', 'TAC', 'GCC', 'GGC', 'GTC', 'GAA', 'GAG', 'GAT'], 
    'CGT': ['AGT', 'GGT', 'TGT', 'CCT', 'CAT', 'CTT', 'CGC', 'CGG', 'CGA'], 
    'GAA': ['CAA', 'AAA', 'TAA', 'GCA', 'GGA', 'GTA', 'GAC', 'GAG', 'GAT'], 
    'TCA': ['CCA', 'GCA', 'ACA', 'TAA', 'TGA', 'TTA', 'TCC', 'TCG', 'TCT'], 
    'GCA': ['CCA', 'ACA', 'TCA', 'GAA', 'GGA', 'GTA', 'GCC', 'GCG', 'GCT'], 
    'GTA': ['CTA', 'ATA', 'TTA', 'GCA', 'GGA', 'GAA', 'GTC', 'GTG', 'GTT'], 
    'GCC': ['CCC', 'ACC', 'TCC', 'GAC', 'GGC', 'GTC', 'GCA', 'GCG', 'GCT'], 
    'GTC': ['CTC', 'ATC', 'TTC', 'GCC', 'GGC', 'GAC', 'GTA', 'GTG', 'GTT'], 
    'GCG': ['CCG', 'ACG', 'TCG', 'GAG', 'GGG', 'GTG', 'GCC', 'GCA', 'GCT'], 
    'GTG': ['CTG', 'ATG', 'TTG', 'GCG', 'GGG', 'GAG', 'GTC', 'GTA', 'GTT'], 
    'TTC': ['CTC', 'GTC', 'ATC', 'TCC', 'TGC', 'TAC', 'TTA', 'TTG', 'TTT'], 
    'GTT': ['CTT', 'ATT', 'TTT', 'GCT', 'GGT', 'GAT', 'GTC', 'GTG', 'GTA'], 
    'GCT': ['CCT', 'ACT', 'TCT', 'GAT', 'GGT', 'GTT', 'GCC', 'GCG', 'GCA'], 
    'TTA': ['CTA', 'GTA', 'ATA', 'TCA', 'TGA', 'TAA', 'TTC', 'TTG', 'TTT'], 
    'TGA': ['CGA', 'GGA', 'AGA', 'TCA', 'TAA', 'TTA', 'TGC', 'TGG', 'TGT'], 
    'TTG': ['CTG', 'GTG', 'ATG', 'TCG', 'TGG', 'TAG', 'TTC', 'TTA', 'TTT'], 
    'TCC': ['CCC', 'GCC', 'ACC', 'TAC', 'TGC', 'TTC', 'TCA', 'TCG', 'TCT'], 
    'TGG': ['CGG', 'GGG', 'AGG', 'TCG', 'TAG', 'TTG', 'TGC', 'TGA', 'TGT'], 
    'TCT': ['CCT', 'GCT', 'ACT', 'TAT', 'TGT', 'TTT', 'TCC', 'TCG', 'TCA']}

codon_synonyms = {
    'CTT': ['CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
    'ATG': [],
    'AAG': ['AAA'],
    'AAA': ['AAG'],
    'ATC': ['ATA', 'ATT'],
    'AAC': ['AAT'],
    'ATA': ['ATC', 'ATT'],
    'AGG': ['AGA', 'CGA', 'CGC', 'CGG', 'CGT'],
    'CCT': ['CCG', 'CCA', 'CCC'],
    'CTC': ['CTG', 'CTA', 'CTT', 'TTA', 'TTG'],
    'AGC': ['AGT', 'TCG', 'TCA', 'TCC', 'TCT'],
    'ACA': ['ACC', 'ACG', 'ACT'],
    'AGA': ['AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'CAT': ['CAC'],
    'AAT': ['AAC'],
    'ATT': ['ATC', 'ATA'],
    'CTG': ['CTA', 'CTC', 'CTT', 'TTA', 'TTG'],
    'CTA': ['CTG', 'CTC', 'CTT', 'TTA', 'TTG'],
    'ACT': ['ACC', 'ACA', 'ACG'],
    'CAC': ['CAT'],
    'ACG': ['ACC', 'ACA', 'ACT'],
    'CAA': ['CAG'],
    'AGT': ['AGC', 'TCG', 'TCA', 'TCC', 'TCT'],
    'CAG': ['CAA'],
    'CCG': ['CCT', 'CCA', 'CCC'],
    'CCC': ['CCT', 'CCG', 'CCA'],
    'TAT': ['TAC'],
    'GGT': ['GGG', 'GGA', 'GGC'],
    'TGT': ['TGC'],
    'CGA': ['AGG', 'AGA', 'CGC', 'CGG', 'CGT'],
    'CCA': ['CCT', 'CCG', 'CCC'],
    'TCT': ['AGC', 'AGT', 'TCG', 'TCA', 'TCC'],
    'GAT': ['GAC'],
    'CGG': ['AGG', 'AGA', 'CGA', 'CGC', 'CGT'],
    'TTT': ['TTC'],
    'TGC': ['TGT'],
    'GGG': ['GGT', 'GGA', 'GGC'],
    'TAG': ['TAA', 'TGA'],
    'GGA': ['GGT', 'GGG', 'GGC'],
    'TAA': ['TAG', 'TGA'],
    'GGC': ['GGT', 'GGG', 'GGA'],
    'TAC': ['TAT'],
    'GAG': ['GAA'],
    'TCG': ['AGC', 'AGT', 'TCA', 'TCC', 'TCT'],
    'TTA': ['CTG', 'CTA', 'CTC', 'CTT', 'TTG'],
    'GAC': ['GAT'],
    'TCC': ['AGC', 'AGT', 'TCG', 'TCA', 'TCT'],
    'GAA': ['GAG'],
    'TCA': ['AGC', 'AGT', 'TCG', 'TCC', 'TCT'],
    'GCA': ['GCC', 'GCG', 'GCT'],
    'GTA': ['GTC', 'GTG', 'GTT'],
    'GCC': ['GCA', 'GCG', 'GCT'],
    'GTC': ['GTA', 'GTG', 'GTT'],
    'GCG': ['GCA', 'GCC', 'GCT'],
    'GTG': ['GTA', 'GTC', 'GTT'],
    'TTC': ['TTT'],
    'GTT': ['GTA', 'GTC', 'GTG'],
    'GCT': ['GCA', 'GCC', 'GCG'],
    'ACC': ['ACA', 'ACG', 'ACT'],
    'TGA': ['TAG', 'TAA'],
    'TTG': ['CTG', 'CTA', 'CTC', 'CTT', 'TTA'],
    'CGT': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG'],
    'TGG': [],
    'CGC': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT']}

neutral_transformation_matrix = {
    'ACC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
    'ATG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'AAG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'AAA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'ATC': {'ACC': 1.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'AAC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 1.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'ATA': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'AGG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
    'CCT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 1.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
    'CTC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'AGC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'ACA': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'AGA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CAT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'AAT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'ATT': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CTG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CTA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'ACT': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
    'CAC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'ACG': {'ACC': 1.0, 'ATG': 1.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CAA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'AGT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CCA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CCG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CCC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
    'TAT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
    'GGT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TGT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 1.0},
    'CGA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CAG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CGC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'GAT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CGG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 1.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
    'CTT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 1.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TGC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 1.0, 'TCT': 0.0},
    'GGG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
    'TAG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 1.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
    'GGA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TAA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'GGC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TAC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
    'GAG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TCG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 1.0, 'TGG': 1.0, 'TCT': 1.0},
    'TTA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TTT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
    'GAC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 1.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'CGT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 1.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'GAA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TCA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 1.0},
    'GCA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'GTA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'GCC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
    'GTC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 1.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'GCG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'GTG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TTC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
    'GTT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'GCT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
    'TGA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
    'TTG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 1.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
    'TCC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
    'TGG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 1.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
    'TCT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0}}
def check_if_position_complete(position, scores_dict, start_pos, wt_DNA_seq):
    check = True
    wt_codon = wt_DNA_seq[(position - start_pos)*3:(position - start_pos)*3+3]
    for j in codon_neighbors[wt_codon]:
        if translation_table[wt_codon] != translation_table[j] and translation_table[j] != '*':
            sub_name = (int(position), translation_table[wt_codon], translation_table[j])
            #print sub_name, scores_dict.get(sub_name)
            if scores_dict.get(sub_name) == None:
                check = False            
    for syn_codon in codon_synonyms[wt_codon]:
        for k in codon_neighbors[syn_codon]:
            if translation_table[syn_codon] != translation_table[k] and translation_table[k] != '*' and translation_table[syn_codon] != '*':
                sub_name = (int(position), translation_table[syn_codon], translation_table[k])
                if scores_dict.get(sub_name) == None:
                    check = False
    if check == False:
        sys.stderr.write("position %s excluded due to missing data\n" % position)
    return check
    
def infer_protein_seq(scores):
    with open(scores,'r') as tsv:
        line_count = 0
        seq_dict = {}
        for line in tsv:
            #print line
            row = line.strip().split('\t')
            #print '\t'.join(row)
            line_count += 1
            if line_count == 1:
                wt_ind = row.index("AAwt")
                pos_ind = row.index("position")
            elif seq_dict.get(row[pos_ind]) == None:
                seq_dict[int(row[pos_ind])] = row[wt_ind]
        #print seq_dict
    tsv.close()
    seq = []
    for k in range(min(seq_dict.keys()),max(seq_dict.keys())+1):
        seq.append(seq_dict.get(k,"-"))
    return min(seq_dict.keys()), ''.join(seq), max(seq_dict.keys()) 

def read_fasta(file):
    seqs = {}
    with open(file,'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                header = line.strip()[1:]
                seqs[header] = ""
            else:
                seqs[header] += line.strip()
    return seqs

def translate(DNA_sequence, frame = 0):
    DNA_sequence = DNA_sequence.upper()
    protein_sequence_list = []
    for i in range(frame,len(DNA_sequence),3):
        protein_sequence_list.append(translation_table[DNA_sequence[i:i+3]])
    return ''.join(protein_sequence_list)

def reverse_translate(protein, CB):
    DNA = []
    for aa in protein:
        DNA.append(numpy.random.choice(CB[aa]['codons'], p = CB[aa]["norm_freqs"]).upper())
    return ''.join(DNA)

def match_protein_sequences(seq1, seq2):
    if len(seq1) != len(seq2):
        sys.stderr.write("Error: Protein sequences don't match length between scores file and fasta\nFrom Fasta:\t")
        sys.stderr.write(seq2+'\nFrom Scores:\t')
        sys.stderr.write(seq1+'\n')
        sys.exit()
    else:
        test = True
        for i in range(0,len(seq1)):
            if seq1[i] != '-' and seq2[i] != '-' and seq1[i] != seq2[i]:
                test = False
    return test

def read_scores_tsv(file_name, lim_high, lim_low):
    All_subs = {}
    line_count = 0
    with open(file_name,'r') as scores:
        for line in scores:
            line_count += 1
            row = line.strip().split('\t')
            if line_count == 1:
                features = row
                pos_ind = features.index("position")
                aa1_ind = features.index("AAwt")
                aa2_ind = features.index("AAmut")
                score_ind = features.index("norm_score")
            else:
                if row[aa1_ind] == row[aa2_ind]:
                    continue
                if row[aa2_ind] == "*":
                    continue
                sub_name = (int(row[pos_ind]), row[aa1_ind], row[aa2_ind])
                if float(row[score_ind]) > lim_high:
                    All_subs[sub_name] = lim_high
                elif float(row[score_ind]) < lim_low:
                    All_subs[sub_name] = lim_low
                else:
                    All_subs[sub_name] = float(row[score_ind])
    return All_subs

def read_transformation_matrix(file_name = None, species = None, process = None):
    if file_name == None:
        if species == None or process == None:
            sys.stderr.write("failed to find transformation matrix\nexiting\n")
            sys.exit()
        else:
            file_name = "./" + species.lower() + "_" + process.lower() + "_transform.txt"
        
    if not os.path.exists(file_name):
        sys.stderr.write("failed to find transformation matrix\nexiting\n")
        sys.exit()
    transform_matrix = {}
    with open(file_name,'r') as mat_file:
        Block = None
        features = {}
        for line in mat_file:
            if line.startswith(">"):
                Block = line.strip()[1:]
            elif Block == "General_info":
                if line.count(':') == 1:
                    features[line.strip().split(':')[0]] = line.strip().split(':')[1]
            
            elif Block == "Array":
                if line.strip().split('\t')[0] == "codon":
                    cols = {}
                    for i in range(1,65):
                        cols[i] = line.strip().split('\t')[i]
                elif line[0] in ['A','C','T','G']:
                    row = line.strip().split('\t')
                    transform_matrix[row[0]] = {}
                    for i in range(1,65):
                        transform_matrix[row[0]][cols[i]] = float(row[i])
        mat_file.close()
    return transform_matrix

def codon_shuffle(DNA_sequence): #Dependencies: module random or method random.shuffle, translation_table from Robustness_scores.py
    codon_library = {}
    for i in range(0,len(DNA_sequence),3):
        if codon_library.get(translation_table[DNA_sequence[i:i+3]]) == None:
            codon_library[translation_table[DNA_sequence[i:i+3]]] = []
        codon_library[translation_table[DNA_sequence[i:i+3]]].append(DNA_sequence[i:i+3])
    for aa in codon_library:
        random.shuffle(codon_library[aa])
     
    new_seq = []
    for i in range(0,len(DNA_sequence),3):
        new_seq.append(codon_library[translation_table[DNA_sequence[i:i+3]]].pop(0))
    return ''.join(new_seq)

def read_codon_freqs(filename):
    codon_freqs = {
            'A': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'C': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'D': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'E': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'F': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'G': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'H': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'I': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'K': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'L': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'M': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'N': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'P': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'Q': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'R': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'S': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'T': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'V': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'W': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            'Y': {'codons':[], 'norm_freqs':[], 'raw_freqs': []},
            }
    with open(filename,'r') as tsv_in:
        lc = 0
        for line in tsv_in:
            lc +=1
            if lc ==1:
                header = line.strip().split('\t')
                codon_ind = header.index('codon')
                aa_ind = header.index('amino_acid')
                rf_ind = header.index('raw_freq')
                nf_ind = header.index('norm_freq')
            else:
                row = line.strip().split('\t')
                codon_freqs[row[aa_ind]]['codons'].append(row[codon_ind])
                codon_freqs[row[aa_ind]]['norm_freqs'].append(float(row[nf_ind]))
                codon_freqs[row[aa_ind]]['raw_freqs'].append(float(row[rf_ind]))
    return codon_freqs

def CRS_score(DNA_seq, position, scores_table, start_pos, tra_mat = None, set_codon = None, soo = "neutral"): # soo being species-of-origin, with name 'species' being used only as a main() variable to avoid confusion, but they mean the same thing
    if tra_mat == None:
        tra_mat = neutral_transformation_matrix
    elif soo == "neutral":
        sys.stderr.write("Error: if a weight matrix is specified, the species must also be specified")
        sys.exit()
    i = 3*(position - start_pos)
    if set_codon == None:
        codon = DNA_seq[i:i+3]
    else:
        codon = set_codon
    sum = 0.0
    denom = 0.0
    for sub in sorted(codon_neighbors.keys()):
        if translation_table[codon] == translation_table[sub]:
            score = 1.0
        elif translation_table[sub] == '*':
            score = 0.0
        else:
            mut = (position, translation_table[codon], translation_table[sub])
            try:
                score = float(scores_table[mut])
            except KeyError:
                score = "missing"
        if tra_mat == None:
            weight = 1
        else:
            weight = tra_mat[codon][sub]    
        if score != "missing":
            sum += (score * weight)
            denom += weight
    CRS = sum / denom
    return CRS

def seq_dist(seq1,seq2,language = "DNA",method = "Hamming"):
    if method == "Hamming":
        dist = 0
        for i in range(0,min([len(seq1),len(seq2)])):
            if seq1[i] != seq2[i]:
                dist +=1
    return dist
