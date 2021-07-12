

# load the rotation centre and offset (fixed)
rotCent, offset = vis.loadInstParams('aug19')

# load the PFI coordinates of the fiducial and science fibres
fiducials, scienceFibres = fps.getFieldDefinition('')

# get the geometry for a given image
za, inr = fps.getInstConfig(fname)

# retrieve the centroids
centroids = vis.getCentroidsDB(conn, frameIDs)
# image=vis.getImage(fname)
# centroids=vis.getCentroids(image,thresh1,thresh2,fwhmx,fwhmy,boxFind,boxCent,nmin,nmax,maxIt)

# transform fiducial fibres
fidPos, sciPos = fps.getFibrePos(fiducials, scienceFibres, za, inr, rotCent, offset)


# match fiducial fibres
fCentroid = mcs.findHomes(centroids, fidPos, tol)


# get the transformation
trans = fps.getAffine(fidPos[:, 1], fidPos[:, 2], fCentroid[:, 1], fCentroid[:, 2])


# apply the above to the science fibres
sCentroid = mcs.findHomes(centroids, sciPos, tol)
nCentroid = fps.applyAffineFPS(sciPos, trans)


# calculate differences
dx, dy = fps.getDiff(sCentroid, nCentroid)

# plot
visplot.quiverPlot(nCentroid[:, 1], nCentroid[:, 2], dx, dy)
