#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import streamlit as st

sample_id = sys.argv[1]
report_file = sys.argv[2]

@st.cache_data
def load_data(data):
    df = pd.read_csv(data)
    return df


report = load_data(report_file)

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

st.title('CRISPR-based gene editing report:')

st.subheader('Allele fractions:')
st.write(report)

st.subheader('Indel rate:')
fig, ax = plt.subplots()
ax.pie(vals, labels = labels, autopct = '%1.0f%%')
ax.legend()
st.pyplot(fig)
plt.clf()

y = [round(x, 3) for x in report.percentage]

st.subheader('Allele heatmap:')
heat = sns.heatmap(vals_df.T.to_numpy(), annot = df.T.to_numpy(), yticklabels = y, linewidth = .5, fmt = '')
st.pyplot(heat.get_figure())
