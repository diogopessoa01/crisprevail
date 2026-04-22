#!/usr/bin/env python

# Import libraries.
import sys
from shiny import reactive
from shiny.express import input, render, ui
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data.
sample_id = 'hCas9-TRAC-a'
report_file = 'hCas9-TRAC-a.csv'

report = pd.read_csv(report_file)

wild_type = 0.0
indel_rate = 0.0

df = pd.DataFrame()
vals_df = pd.DataFrame()

cigar_dict = {
        'M': 0.0,
        'I': 1.0,
        'D': 2.0,
        'N': 3.0,
        'S': 4.0,
        'H': 5.0,
        'P': 6.0,
        '=': 7.0,
        'X': 8.0,
        'B': 9.0
}

for i in range(report.shape[0]):
    allele = list(report.iloc[i, 0])
    allele_vals = [cigar_dict[x] for x in allele]
    allele_keys = list(map(lambda x: '-' if x == 'D' else x, allele))
    df = pd.concat([df, pd.Series(allele_keys)], axis = 1, ignore_index = True)
    vals_df = pd.concat([vals_df, pd.Series(allele_vals)], axis = 1, ignore_index = True)
    if set(report.iloc[i, 0]) == {'M'}:
        wild_type = report.iloc[i, 2]


indel_rate = 100 - wild_type

vals = [wild_type, indel_rate]
labels = ['Unmodified%', 'Modified%']

y = [round(x, 3) for x in report.percentage]

# app

ui.tags.script(
    src="https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
)
ui.tags.script("if (window.MathJax) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);")


with ui.div(class_="col-md-10 col-lg-8 py-5 mx-auto text-lg-center text-left"):
    ui.h3("CRISPR-based gene editing report:")

with ui.div(class_="col-lg-11 py-5 mx-auto text-center"):
    ui.h2("Gene editing rate:")

    @render.plot()
    def plot_pie():
        fig, ax = plt.subplots()
        ax = plt.pie(vals, labels = labels, autopct = '%1.0f%%')
        return fig

with ui.div(class_="col-lg-11 py-5 mx-auto text-center"):
    ui.h2("Allele fraction:")

    @render.plot()
    def plot_heatmap():
        fig, ax = plt.subplots()
        ax2 = sns.heatmap(vals_df.T.to_numpy(), annot = df.T.to_numpy(), yticklabels = y, linewidth = .5, fmt = '')
        return fig


