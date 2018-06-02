from __future__ import division
from collections import defaultdict
import os,sys
import StringIO

import random
import copy
import math
import multiprocessing
import time
import sys
import csv
import numpy as np


# Short-time-series distance function
def sts(t,x,y): #t is the time interval (1/120s), x is the benchmark time series, y is the test time series
	inside = np.zeros(len(x))
	for i in range(len(x)-1):
		print x[i+1]-x[i] 
		inside[i] = np.square(((x[i+1]-x[i])/t)-((y[i+1]-y[i])/t))
	middle = sum(inside)
	answer = np.sqrt(middle)
	return answer