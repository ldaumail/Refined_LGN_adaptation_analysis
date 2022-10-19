import scipy.io
import pandas as pd

fname = 'concat_origfilenames_completenames_10162022.mat'
mat = scipy.io.loadmat(fname)

dates = list()
p = list()
clusts = list()
monk = list()
for key in mat.keys():
    if key[0] == 'x':
        dates.append(key[1:7])
        monk.append(key[8])
        p.append(key[12:13])
        clust = key[(key.index('uclust') + len('uclust')):(key.index('uclust') + len('uclust') +2)]
        if clust.__contains__('_'):
            clust = int(clust[0])
        else:
            clust = int(clust)
        clusts.append(clust)

df = pd.DataFrame({'date':dates,'penetrationN':p,'clust':clusts,'monkey':monk})

df.head()
with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(df)

n_days = df['date'].nunique()
print('total number of days ' + str(n_days))

#total number of units by day, not accounting for different depths

number_of_units_by_day = df['date'].value_counts()
number_of_units_by_day
df2 = number_of_units_by_day.to_frame().reset_index() 
df2.columns = ['date', 'Number of Units']
df2.sort_values(by=['date'])

number_of_units_per_depth = df.groupby('date')['penetrationN'].value_counts()

# this is telling you the number of units at each depth in the 3 column that you see//all grouped by day
#(p01, p02..correspond to different rounds. i think i almost always moved it 
#or it had moved on its own so we should just call it that)

number_of_units_per_depth
df3 = number_of_units_per_depth.to_frame() #.reset_index() 
df3.columns = ['date', 'Penetration Number','Number of Units']

