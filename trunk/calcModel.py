import rpy2.robjects as robjects
import math

## Calculates the mean using Nadaraya-Watson weights
# @param p The bandwidth variable used in the weight
# @param n The sample size to calculate the mean for
# @param points The data that will be used for calculating the means
# /return mean The calculated mean
def calcMean(p, n, points):
    denom = 0
    numer = 0
    for i in points:
        for j in points[i]:
            denom += math.exp(-1.0/p*pow(i-n,2))
            numer += math.exp(-1.0/p*pow(i-n,2)) * j
##    for i in points:
##        denom += math.exp(-1.0/p * pow(i-n, 2))
##        numer += math.exp(-1.0/p * pow(i-n, 2)) * robjects.r['ave'](robjects.FloatVector(points[i]))[0]
    return numer/denom

## Calculates the variance using Nadaraya-Watson weights
# @param p The bandwidth value used in the weight
# @param n The sample size to calculate the variance for
# @param points The data that will be used for calculating the variance
# @param means Dictionary that contains calculated mean values for the corresponding sample sizes
# /return variance The calculated variance
def calcVar(p, n, points, means):
    denom = 0
    numer = 0
    for i in points:
        for j in points[i]:
            denom += math.exp(-1.0/p*pow(i-n,2))
            numer += math.exp(-1.0/p*pow(i-n,2)) * pow((j-means[i]),2)
##    for i in points:
##        denom += math.exp(-1.0/p * pow(i-n, 2))
##        numer += math.exp(-1.0/p * pow(i-n, 2)) * pow(robjects.r['ave'](robjects.FloatVector(points[i]))[0] - means[i], 2)
    return numer/denom

## Calculates the means and variances for a range of values using the given data
# @param data The data that will be used for calculating the means and variances
# @param minimum The minimum value for which to calculate the mean and variance, default is 0
# @param maximum The maximum value for which to calculate the mean and variance, default is 200
# @param p The bandwidth value used in the Nadaraya-Watson weight
# /return model Dictionary that contains tuples of mean and variance for each sample size
def calcModel(data, minimum = 0, maximum = 200, p=10):
    model = dict()
    means = dict()
    for i in range(minimum, maximum+1):
        means[i] = calcMean(p, i, data)
    for i in range(minimum, maximum+1):
        model[i] = (means[i], calcVar(p, i, data, means))
    return model

## Calculates the probability using the inputed values and calculated means and variances
# @param y The total information loss
# @param n The sample size (number of merged genes, nodes, etc.)
# @param model The model that contains the calculated means and variances
# /return probability The probability that the given amount of information loss for the sample size is due to random merging
def calcProb(y, n, model):
    return robjects.r['pnorm'](y, model[n][0], math.sqrt(model[n][1]))[0]

