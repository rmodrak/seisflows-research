
raise NotImplementedError

# compute azimuth
try:
    parts = solver.load(path +'/'+ 'kernels/sum', type='kernel')
    for iproc in range(PAR.NPROC):
        y = parts['Gs'][iproc]
        x = - parts['Gc'][iproc]
        t = 0.5*np.arctan2(y, x)
        dirname = path +'/'+ tag
        filename = 'proc%06d_%s.bin' % (iproc, 'azimuth')
        savebin(t, dirname +'/'+ filename)
except:
    pass

