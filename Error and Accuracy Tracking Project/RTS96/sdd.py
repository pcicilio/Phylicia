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

def sdd(x,y):
	T = 1/120*len(x)
	d_inner = np.zeros(len(x))
	s_inner = np.zeros(len(x))
	s_outer = np.zeros(len(x))
	for i in range(len(x)-1):
		d_inner[i] = abs(x[i]-y[i])
		s_inner[i] = x[i]-y[i]
	d = 1/T*sum(d_inner)
	s_total = sum(s_inner)
	for i in range(len(x)-1):
		s_outer[i] = (s_inner[i]-s_total)
	s = 1/T*sum(s_outer)
	return 1/2*(d+s)
	