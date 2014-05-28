#!/usr/bin/env python

import numpy as np
import sys

__doc__ = "GLM fitting with weighted data points"

def optimize_w ( X, y, q, w, niter=5, stop=1e-5, lm=0.1, verbose=False ):
    """Optimize the w parameters of the model with weighted data

    :Parameters:
        *X*     design matrix
        *y*     responses
        *q*     state probabilities (weights)
        *w*     starting weights

    :Optional parameters:
        *niter* number of iterations
        *stop*  stopping criterion
        *lm*    regularization
        *verbose* obvious
    """
    w,l_,conv = weighted_glm ( X, y, q, w, niter=niter, stop=stop, lm=lm )
    if verbose:
        print "Converged" if conv else "Not converged"
    return np.ravel(np.array(w))

def logistic ( x ):
    """Logistic function"""
    return 1./(1+np.exp(-x))

def weighted_glm ( X, y, q, w, niter=5, stop=1e-5, lm=0.1 ):
    """The actual optimization

    Parameters: see optimize_w
    """
    X = np.matrix ( X )
    w = np.matrix ( w.reshape((-1,1)) )
    y = np.matrix ( y.reshape((-1,1)) )
    q = np.matrix ( q.reshape((-1,1)) )
    eta = logistic ( X*w )
    l = np.sum(q.A*(y.A*np.log(eta.A)+(1-y.A)*np.log(1-eta.A))) - w.T*w

    for iteration in xrange ( niter ):
        z = q.A*(y-eta).A
        grad_l = X.T*z - 2*lm*w
        w_ii = -q.A*eta.A*(1-eta.A)
        H = (X.A*w_ii.reshape ( (-1,1) )).T * X
        H -= 2*lm*np.eye ( len(w) )

        dw = np.linalg.lstsq ( H, grad_l )[0]
        e = np.sqrt ( np.sum ( (dw).A**2 ) )

        w -= 0.9*dw
        eta = logistic ( X*w )

        if e < stop:
            break

    l_ = np.sum(q.A*(y.A*np.log(eta.A)+(1-y.A)*np.log(1-eta.A)))
    wnorm = - w.T*w

    if not l_-wnorm>=l:
        conv = 0.
    else:
        conv = 1.

    return w,l_,conv

