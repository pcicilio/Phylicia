
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
from scipy.spatial.distance import euclidean

from fastdtw import fastdtw

x = np.array([[1,1], [2,2], [3,3], [4,4], [5,5]])
y = np.array([[2,2], [3,3], [4,4]])
distance, path = fastdtw(x, y, dist=euclidean)
print(distance)
