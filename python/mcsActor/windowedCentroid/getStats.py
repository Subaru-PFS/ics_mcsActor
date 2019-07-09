


def getThresh(image,factor1,factor2):
    mn=image.mean()
    st=image.std()

    thresh1=mn+st*factor1
    thresh2=thresh1*factor2

    return thresh1,thresh2


    
