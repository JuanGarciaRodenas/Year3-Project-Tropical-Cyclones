#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 16:03:03 2023

@author: juan
"""
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


xValues = np.array([31, 15.5, 3.6, -15.5, -32, -56, -81, -63])
yValues = np.array([1000, 800, 600, 400, 300, 200, 100, 50])

XValues = deepcopy(xValues)
YValues = deepcopy(yValues)


index = 0

for x in range(len(yValues) - 1):
    X = np.array([])
    Y = np.arange(yValues[x + 1] + 1, yValues[x], step = 1)
    Y = np.flip(Y)
    for i in range(len(Y)):
        Value = ((yValues[x] - Y[i]) * xValues[x + 1] - (yValues[x + 1] - Y[i]) * xValues[x]) * (1 / (yValues[x] - yValues[x + 1]))
        X = np.append(X, np.array([Value]))
    if x == 0:
        index = 1
    else:
        index = index + (yValues[x - 1] - yValues[x])
    XValues = np.insert(XValues, index, X)
    YValues = np.insert(YValues, index, Y)
        


            
            


plt.figure(figsize = (3, 4))
plt.plot(XValues, YValues)
plt.ylim(1000, 0)

