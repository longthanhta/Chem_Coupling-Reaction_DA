# Data-driven studying for chemistry

This repository containing python script for analysis, data mining and machine learning.

## 1 Retrosynthesis

A simple retrosynthesis approach based on template (obtained from ASKCOS project), intergrateracd Fastfoward neuralnetwork help to suggesting a generated reaction is reasonable or not.

## 2 Bond extraction

Based on https://pubs.acs.org/doi/10.1021/ci400442f, this folder containing python script to extract bond information from atom-mapped reaction in SMILES format (SMILES that does not containt reaction map number will not work)


## 3 Feature Selections

Contains different tools used for selecting feature, including LASSO, PCA, pearson.

## 4 Mining tool

Containing graph mining tool to extract information from article, run on Linux and require OSRA and Tesseract OCR pre-installed

## 5 Raw data processing tool

Containing tool to handle xml file extracted from Reaxys

## 6 Distribution

Python script to generating distribution plot based on extracted bond

## 7 Network analysis

Based on https://www.nature.com/articles/nchem.136, conduct graph analysis network for reactions.

## 8 Unsupervised clustering

Contain script for conducting k-means and visualizing interactive plot in html (using plotly).
