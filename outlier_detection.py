from scipy.spatial import distance
import numpy as np

DEFAULT_K = 17
DEFAULT_B = 5
DEFAULT_FAR = 0.05
DEFAULT_TOLERANCE = 0.001

def get_mad(data):
    # Calculate median and median absolute deviation (mad) for data

    med = np.median(data)
    mad = np.median(abs(data - med))

    return med, mad

def mads_away(query, data):
    # Obtain how many median absolute deviations the query point is away from the median

    med, mad = get_mad(data)
    mads = abs(query - med) / mad

    return mads

def k_mads_away(query, data, k=DEFAULT_K):
    # Calculations how many median absolute deviation the query point
    # is away from the median of the k nearest neighbors

    distances = distance_to_n(query, data, sort=False)
    knn_idx = np.argsort(distances)[0:k]
    mads = mads_away(query, data[knn_idx])
    
    return mads

def calculate_G(distances, k=DEFAULT_K):
    # Calculate G statistic for query

    lower_bound = int(k - np.floor((k-1)/2) - 1)
    upper_bound = int(k + np.floor(k/2))

    G = np.average(np.sort(distances)[lower_bound:upper_bound])

    return G

def train_aKLPE(data, k=DEFAULT_K, b=DEFAULT_B):
    # TODO: see if it can be vectorized/parallelized
    # calculating distances several times...
    n = len(data)
    gs = np.zeros(n)
    
    distances = distance.cdist(np.reshape(data, (-1,1)), np.reshape(data, (-1,1)))

    for i in range(0, b):
        # randomly split data into two equal parts
        set_1, set_2 = split_set(n)

        # use S2 to calculate G for xi in S1 and viceversa
        # TODO: sort before sending to calculate G?
        for i in set_1:
            gs[i] += calculate_G(distances[i, set_2], k=k)
            
        for i in set_2:
            gs[i] += calculate_G(distances[i, set_1], k=k)

    gs = np.divide(gs, b)

    return gs
    
def test_aKLPE(query, data, data_gs=None, k=DEFAULT_K, b=DEFAULT_B, far=DEFAULT_FAR):
    n = len(data)
    g = 0
    gs = data_gs
    is_outlier = False

    distances = distance_to_n(query, data, sort=False)
    
    if gs is None:
        gs = train_aKLPE(data, k, b)
    
    for i in range(0, b):
        subset = np.random.choice(n, np.floor(n/2), replace=False)
        g += calculate_G(distances[subset], k=k)
        
    g /= b

    p = np.average( g < gs )

    if p < far:
        is_outlier = True

    return is_outlier, p

def get_nr(query, data, k=DEFAULT_K, tol=DEFAULT_TOLERANCE):

    distances = distance_to_n(query, data, sort=False)
    knn_idx = np.argsort(distances)[0][:k]
    knn_distances = distances[0][knn_idx]
    knn = data[knn_idx]
    knn_distances[np.where(knn_distances == 0)[0][0]] = tol

    weight = np.divide(1, np.multiply(np.sum(knn_distances), knn_distances)) 
    
    nr = np.abs(query - np.median(np.multiply(weight, knn)))
    med, mad  = get_mad(knn)

    if nr != 0:
        adj = tol * nr
    else:
        adj = 1 * nr

    nr /= mad + adj

    return nr

def train_nr(data, k=DEFAULT_K, b=DEFAULT_B, tol=DEFAULT_TOLERANCE, far=DEFAULT_FAR):
    # TODO: same idea where I calculate distances several times

    n = len(data)
    rs = np.zeros(n*b)
    rs_idx = 0
    dt = far


    k = np.floor(n/2) if k > np.floor(n/2) else k # according to paper's authors

    for i in range(0, b):
        set_1, set_2 = split_set(n)

        for i in set_1:
            rs[rs_idx] += get_nr(data[i], data[set_2], k=k, tol=tol)
            rs_idx += 1

        for i in set_2:
            rs[rs_idx] += get_nr(data[i], data[set_1], k=k, tol=tol)
            rs_idx += 1

    rs = np.sort(rs)
    m = np.round(far * n * b)
    
    dt *= rs[m - n*b]

    return dt

def test_nr(query, data, h, k=DEFAULT_K, far=DEFAULT_FAR, tol=DEFAULT_TOLERANCE):
    is_outlier = False
    n = len(data)

    set_1, set_2 = split_set(n)

    nr = get_nr(query, data[set_1], k=k, tol=tol)
    
    if nr > h:
        is_outlier = True

    return is_outlier, nr

def distance_to_n(query, data, sort=True):
    distances = distance.cdist([[query]], np.reshape(data, (-1,1)))
    
    if sort:
        distances = distances.sort()
    
    return distances

def split_set(n):
    idx = np.random.choice(n, n, replace=False)
    bound = int(np.floor(n/2))
    set_1 = idx[0:bound]
    set_2 = idx[bound:]

    return set_1, set_2
