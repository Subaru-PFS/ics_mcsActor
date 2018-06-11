
import numpy as np

def transform(primary,secondary):

    n = primary.shape[0]
    pad = lambda x: np.hstack([x, np.ones((x.shape[0], 1))])
    unpad = lambda x: x[:,:-1]
    X = pad(primary)
    Y = pad(secondary)

    A, res, rank, s = np.linalg.lstsq(X, Y)

    transform = lambda x: unpad(np.dot(pad(x), A))

    B=transform(primary)

    print(B)

primary = np.array([[41., 1158., 0.],
                    [39., 40., 0.],
                    [261., 42., 0.],
                    [256., 1162., 0.]])

secondary = np.array([[610., 560., 0.],
                      [610., -560., 0.],
                      [390., -560., 0.],
                      [390., 560., 0.]])

transform(primary,secondary)
