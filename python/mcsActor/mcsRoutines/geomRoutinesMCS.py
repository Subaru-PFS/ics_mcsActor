
import os
import yaml


def defaultCentroidParams(actorName):
    """
    read the default centroid parameters files. note that the yaml file contains both 
    values and definitions of the parameters.

    Input
       actorName: actor in question. currently agc or mcs

    Returns
      dictionary with the values of the parameters

    """
    if(actorName == "mcs"):
        # currently in the /etc directory of the ics_mcsActor
        path = os.path.join("$MCSACTOR_DIR", "etc", "mcsDefaultCentroidParameters.yaml")
    elif(actorName == "agc"):
        # currently in the /etc directory of the ics_mcsActor
        path = os.path.join("$AGCACTOR_DIR", "etc", "agcDefaultCentroidParameters.yaml")
    else:
        return None

    with open(path, 'r') as inFile:
        defaultParms = yaml.safe_load(inFile)

    # returns just the values dictionary
    return defaultParms['values']


def readCobraGeometry(xmlFile, dotFile):
    """
    read cobra geometry from configuration file/inst_config

    dot positions from CSVfile at the moment, this will change

    The results will be return in whatever unit the input XML file is in
    """

    # geometry XML file

    pfic = pfi.PFI(fpgaHost='localhost', doConnect=False, logDir=None)
    aa = pfic.loadModel([pathlib.Path(xmlFile)])

    # first figure out the good cobras (bad positions are set to 0)
    centersAll = pfic.calibModel.centers
    goodIdx = np.array(np.where(centersAll.real != 0)).astype('int').ravel()

    # then extract the parameters for good fibres only
    centrePos = np.array([goodIdx+1, centersAll[goodIdx].real, centersAll[goodIdx].imag]).T
    armLength = (pfic.calibModel.L1[goodIdx]+pfic.calibModel.L2[goodIdx])

    # number of cobras
    nCobras = len(armLength)

    # at the moment read dots from CSV file
    dotData = pd.read_csv(dotFile, delimiter=",")
    dotPos = np.zeros((len(goodIdx), 4))

    dotPos[:, 0] = goodIdx+1
    dotPos[:, 1] = dotData['x_tran'].values[goodIdx]
    dotPos[:, 2] = dotData['y_tran'].values[goodIdx]
    dotPos[:, 3] = dotData['r_tran'].values[goodIdx]

    return centrePos, armLength, dotPos, goodIdx
