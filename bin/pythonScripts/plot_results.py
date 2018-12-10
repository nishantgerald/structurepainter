import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import pymongo

import collections

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def flatten(d, parent_key='', sep='_'):
    """
    flatten nested Mongodb records
    """
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

client = pymongo.MongoClient()
db = client.lai
results = db.results

df = pd.DataFrame(list(map(flatten, results.find({'batch_size' : 20, 'window_size' : {'$gt' : 2}}))))
fig, ax = plt.subplots()

ax.scatter(df['window_size'], df['overall_acc_avg'], color='red')
ax.set_ylabel('average accuracy\n(across 20\ntest chromosomes)', color='red',
    rotation='horizontal', horizontalalignment='right', multialignment='center')

ax2 = ax.twinx()
ax2.scatter(df['window_size'], df['time'] / 60, color='blue')
ax2.set_ylabel('runtime\n(minutes)', color='blue',
    rotation='horizontal', horizontalalignment='left', multialignment='center')

ax.set_xlabel('window size')
plt.title('window size vs. accuracy and runtime (using sliding windows)')
plt.gcf().set_size_inches([9.4, 4.8])

fig.savefig('results/sliding_windows/window_size_vs_accuracy_and_runtime')
