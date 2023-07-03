from scipy import io
figureFolder = os.path.expanduser('~/Data/Matlab-DRIVE/Shared/figures/subspaceAngle/')
plt.close('all')
normalize = True
limitRange = False

for k in range(1,7):

    filename = os.path.join(figureFolder, f'subspaceDist_K={k}.mat')

    objs = io.loadmat(filename, 
                      squeeze_me=True)

    X = pd.DataFrame(data=objs['subspaceDist'], index=objs['rowVar'], columns=objs['rowVar'])
    X = 1-X.abs()
    y = [q.replace('\\theta','$\\theta$').replace('\\delta','$\\delta$').replace('HH','$HH$').replace('HP','$HP$').replace('SPW-Rs','SPW-R') for q in tuple(X.index)]

    if normalize:
        X = (X-X.min())/(X.max()-X.min())

    X.index = X.columns=pd.Index(tuple(y));

    if limitRange:
        M = np.array(X.values)
        for i in range(len(X)):
            M[i,i] = np.nan
        M = np.abs(np.nanmax(M))
        m = np.abs(np.nanmin(M))
    else:
        M, m = np.nanmax(X), np.nanmin(X)

    g = sns.clustermap(X, cmap="Reds", vmin=m, vmax=M)
    g.fig.suptitle(f"K = {k}")
    g.fig.savefig(os.path.join(figureFolder, f"K = {k}.png"))
    g.fig.savefig(os.path.join(figureFolder, f"K = {k}.svg"))
