#import libraries
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
import scipy
from math import sqrt
from scipy.spatial import distance
from scipy.spatial.distance import cdist
from scipy.stats import chi2_contingency
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import rbf_kernel

def scale_data(df) :
    """Scale a dataframe to get the values between 0 and 1. It returns the scaled dataframe.
    
    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        The dataframe to scale

    Returns
    -------
    pandas.core.frame.DataFrame
        A dataframe with the scaled data
    """

    #initialize and fit the scaler
    scaler = MinMaxScaler()
    scaled = scaler.fit_transform(df)

    #return the scaled dataframe
    return pd.DataFrame(scaled, columns=df.columns.tolist())


def pairwise_euclidean_distance(synthetic_data, real_data) :
    """Compute the pairwise euclidean distances between each pair of real and synthetic records
    
    Parameters
    ----------
    synthetic_data : numpy.ndarray
        Synthetic data records
    real_data : numpy.ndarray
        Real data records

    Returns
    -------
    string
        a string with the mean and std values of the computed pairwise euclidean distances
    """

    #compute the pairwise euclidean distances
    distances = distance.cdist(synthetic_data, real_data, 'euclidean')

    #return the mean and std values of the computed pairwise euclidean distances
    return str(np.round(np.mean(distances),4)) + ' ± ' + str(np.round(np.std(distances),4))


def pairwise_mahalanobis_distance(synthetic_data, real_data) :
    """Compute the pairwise mahalanobis distances between each pair of real and synthetic records

    Parameters
    ----------
    synthetic_data : numpy.ndarray
        Synthetic data records
    real_data : numpy.ndarray
        Real data records

    Returns
    -------
    string
        a string with the mean and std values of the computed pairwise mahalanobis distances
    """

    # Compute the covariance matrix from real data
    cov_matrix = np.cov(real_data, rowvar=False)

    # Add small value to diagonal for numerical stability
    cov_matrix += np.eye(cov_matrix.shape[0]) * 1e-6

    # Compute inverse of covariance matrix
    try:
        inv_cov_matrix = np.linalg.inv(cov_matrix)
    except np.linalg.LinAlgError:
        # If singular, use pseudo-inverse
        inv_cov_matrix = np.linalg.pinv(cov_matrix)

    # Compute pairwise Mahalanobis distances
    distances = distance.cdist(synthetic_data, real_data, 'mahalanobis', VI=inv_cov_matrix)

    # Return the mean and std values
    return str(np.round(np.mean(distances),4)) + ' ± ' + str(np.round(np.std(distances),4))


# def hausdorff_distance(synthetic_data, real_data) :
#     """Compute the hausdorff distance between synthetic dataset and real dataset
    
#     Parameters
#     ----------
#     synthetic_data : numpy.ndarray
#         Synthetic data records
#     real_data : numpy.ndarray
#         Real data records

#     Returns
#     -------
#     float
#         the hausdorf distance between synthetic and real datasets
#     """

#     #compute the hausdorff distance
#     hausdorff_dist = scipy.spatial.distance.directed_hausdorff(synthetic_data, real_data)[0]

#     #return the computed value rounded on 4 decimals
#     return np.round(hausdorff_dist,4)


def hausdorff_distance(synthetic_data, real_data):
    #更新使用马氏距离算
    cov = np.cov(real_data, rowvar=False)
    cov += np.eye(cov.shape[0]) * 1e-6
    try:
        VI = np.linalg.inv(cov)
    except:
        VI = np.linalg.pinv(cov)

    # pairwise mahalanobis distance
    D = cdist(synthetic_data, real_data,
              metric="mahalanobis", VI=VI)

    # directed hausdorff
    d_AB = np.max(np.min(D, axis=1))
    d_BA = np.max(np.min(D, axis=0))
    hausdorff_dist = max(d_AB, d_BA)

    return np.round(hausdorff_dist, 4)


def rts_similarity(synthetic_data, real_data) :
    """Compute the real to synthetic similarity between each pair of synthetic and real records
    
    Parameters
    ----------
    synthetic_data : numpy.ndarray
        Synthetic data records
    real_data : numpy.ndarray
        Real data records

    Returns
    -------
    dict
        a dictionary with the min, mean and max values of the computed similarity values
    """

    #compute the cosine similarity between each pair of synthetic and real records
    str_sim = cosine_similarity(synthetic_data, real_data)

    #return a dictionary with the min, mean and max values of the computed similarity values
    return {'min' : np.round(np.min(str_sim),4), 'mean' : np.round(np.mean(str_sim),4), 'max' : np.round(np.max(str_sim),4)}


def MMD(X, Y, gamma=1.0):
    Kxx = rbf_kernel(X, X, gamma)
    Kyy = rbf_kernel(Y, Y, gamma)
    Kxy = rbf_kernel(X, Y, gamma)

    return Kxx.mean() + Kyy.mean() - 2*Kxy.mean()