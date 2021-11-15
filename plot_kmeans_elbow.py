import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans

human_fn = 'human_zscore.txt'
macaque_fn = 'macaque_zscore.txt'

def get_distortions(fn, n=10):

    df = pd.read_csv(fn,sep=' ',index_col=[0])
    df = df.T
    #df = pd.melt(df, id_vars=['Area'])
    distortions = []
    K = range(1,n)
    for k in K:
        kmeanModel = KMeans(n_clusters=k)
        #kmeanModel.fit(df['value'].values.reshape(-1,1))
        kmeanModel.fit(df.values)
        distortions.append(kmeanModel.inertia_)

    return K, distortions


K, hd = get_distortions(human_fn)
K, md = get_distortions(macaque_fn)

plt.figure(figsize=(12,9))

plt.plot(K, hd, 'bx-', label='Human')
plt.plot(K, md, 'ro-', label='Macaque')
plt.xlabel('k',size=20)
plt.ylabel('Distortion',size=20)
plt.legend( prop={'size': 20} )

plt.title('Elbow plot',size=20)

plt.savefig('kmeans_elbow_plot.png')

