#! /usr/bin/env python
# -*- coding: utf-8 -*-
##  ==========================================================================
##                Project: cfd - POD - Copyright (c) 2009 by CERFACS
##  Type   :
##  File   : matrixalgebra.py
##  Vers   : V1.0
##  Chrono : No  Date       Author                 V   Comments
##           1.0 11/08/2009 Braconnier             0.1 Creation
##  ==========================================================================

from numpy import shape, diagonal, sqrt, float, ones, zeros, transpose, identity, dot
import numpy.linalg as nl
#


def findindex(S, tol, maxp):
    (ns, ms) = shape(S)
    xx = diagonal(S)
    sx = sum(xx, 0)
    for i in range(ms):
        index = i + 1
        sumpa = sum(xx[:i + 1], 0) / sx
        if sumpa >= tol:
            break
    index = min(index, maxp)
    return index


#


def norm(x):
    n = len(x)
    normx = 0.
    for i in range(n):
        normx = normx + x[i] * x[i]
    normx = sqrt(normx)
    return normx


#


def mat_yy(dim):
    n = 2 ** dim
    yy = zeros([n, dim], Float)
    for j in range(dim):
        k = 2 ** (j + 1)
        nk = n / k
        for i in range(k):
            yy[i * nk:(i + 1) * nk, j:j + 1] = (-1) ** (i + 1) * ones([nk, 1],
                    Float)
    return yy


#


def svd(A):
    (U, S, V) = singular_value_decomposition(A)
    n = len(S)
    S = identity(n, Float) * S
    V = transpose(V)
    return (U, S, V)


#


def dtpr(x):
    n = len(x)
    psx = 0.
    for i in range(n):
        psx = psx + x[i, 0] * x[i, 0]
    return psx


#


def eig(A):
    (D, V) = Heigenvectors(A)
    n = len(D)
    Eig = zeros([n, n], Float)
    for i in range(n):
        Eig[i, i] = sqrt(abs(D[i]))
    V = transpose(V)
    return (Eig, V)


#


def house(a):
    (n, m) = shape(a)
    v = zeros([n, 1], Float)
    for i in range(n):
        v[i, 0] = a[i, 0]
    na = norm(a)
    if na > 1e-15:
        if a[0, 0] >= 0.:
            beta = a[0, 0] + na
        else:
            beta = a[0, 0] - na
        v[1:n, 0] = v[1:n, 0] / beta
    v[0, 0] = 1.
    return v


#


def rowhouse(A, v):
    (n, m) = shape(A)
    A1 = zeros([n, m], Float)
    for i in range(n):
        for j in range(m):
            A1[i, j] = A[i, j]
    beta = -2. / dtpr(v)
    w = beta * matrixmultiply(transpose(A), v)
    A1 = A1 + matrixmultiply(v, transpose(w))
    return A1


#


def qr(A):
    (m, n) = shape(A)
    A1 = zeros([m, n], Float)
    for i in range(m):
        for j in range(n):
            A1[i, j] = A[i, j]
    v = zeros([m, 1], Float)
    for j in range(n):
        v[j:m, 0:1] = house(A1[j:m, j:j + 1])
        A1[j:m, j:n] = rowhouse(A1[j:m, j:n], v[j:m, 0:1])
        if j < m - 1:
            A1[j + 1:m, j:j + 1] = v[j + 1:m, 0:1]
    R = zeros([n, n], Float)
    for i in range(n):
        for j in range(i, n):
            R[i, j] = A1[i, j]
    Q = identity(m, Float)
    for j in range(n - 1, -1, -1):
        v = zeros([m, 1], Float)
        v[j, 0] = 1.
        v[j + 1:m, 0] = A1[j + 1:m, j]
        Q[j:m, j:m] = rowhouse(Q[j:m, j:m], v[j:m, 0:1])
    Q = Q[:, :n]
    return (Q, R)


#
