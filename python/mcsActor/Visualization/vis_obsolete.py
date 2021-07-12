
def approximateCoords(x, y, region):
    """

    Get lower left position and rotation angle from a set of pin-hole mask
    centroids.  This will possibly break if the rotation angle is too high. 

    input

    x,y: 1D numpy arrays of coordinates (pixel)

    region: 4 element array with of the region of interest [x1,x2,y1,y2]

    returns: x1,y1 of lower left corner, angle in radians

    """

    # find points in the ergion

    ind = np.where((x > region[0]) & (x < region[1]) & (y > region[2]) & (y < region[3]))

    # find the minimum and maximum distance from the origin

    dd = x*x+y*y

    ind1 = np.where(dd == dd.max())
    ind2 = np.where(dd == dd.min())

    # calculate the angle of rotation wrt Y axis

    angle = np.arctan2((y[ind1]-y[ind2]), (x[ind1]-x[ind2]))-np.pi/4.

    return x[ind2], y[ind2], angle


def scaleCentroids(x, y, x1, y1, scale):
    """

    scale the centroids to mm at mask, and shift to the origin

    input

    x1,y1: lower left spot
    x,y: coordinates of centroids
    scale: scale factor

    returns: transformed x,y

    """

    xc = (x-x1)/scale
    yc = (y-y1)/scale

    return xc, yc


def scaleMask(xx, yy, angle, flip):
    """

    rotate the mask as needed

    input

    xx,yy: hole positions
    angle: rotation angle
    flip: set to 1 to rotate by 90 degrees

    returns: transformed x,y

    """

    # apply any rotation (for matching purposes only)

    # add 90 degrees to rotation

    if(flip == 1):
        angle = angle+np.pi/2

    # shift to centre, rotate, shift back
    xx = xx-168
    yy = yy-168

    xnew = xx*np.cos(angle)-yy*np.sin(angle)
    ynew = xx*np.sin(angle)+yy*np.cos(angle)

    xx = xnew+168
    yy = ynew+168

    return xx, yy


def maskImage(infile, outfile, x1, y1, x2, y2):

    image = pf.getdata(infile)
    image[0:x1, :] = 0
    image[x2:, :] = 0
    image[:, 0:y1] = 0
    image[:, y2:] = 0
    pf.writeto(outfile, image, clobber=True)
