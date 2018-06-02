# File:"C:\Users\psse\Desktop\Phylicia\Load Placement Study\Test Systems\Poland\POLAND2383DYNP\Codes\RawLoadReductionBus10.py", generated on SUN, OCT 08 2017  19:06, PSS(R)E release 34.02.00
# Run Generator Outages on all generators using the 6 dynamic files: ZIP, benchmark (50%) CLM, 1Q, 1Q, 3Q, 4Q

from __future__ import division
from collections import defaultdict
import os,sys
import StringIO

#PSSE_PATH = r"C:\Program Files (x86)\PTI34\PSSBIN"
#if not os.path.exists(PSSE_PATH):
#	print "bad path"

#Change to your PSS/e Location and set up paths
sys.path.append(r"C:\Program Files (x86)\PTI\PSSE34\PSSPY27") #Give the path to PSSBIN to imoport psspy
sys.path.append(r"C:\Program Files (x86)\PTI\PSSE34\PSSBIN")
sys.path.append(r"C:\Program Files (x86)\PTI\PSSE34\PSSLIB")
sys.path.append(r"C:\Program Files (x86)\PTI\PSSE34\EXAMPLE")
os.environ['PATH'] = (r"C:\Program Files (x86)\PTI\PSSE34\PSSPY27;" + r"C:\Program Files (x86)\PTI\PSSE34\PSSBIN;" + r"C:\Program Files (x86)\PTI\PSSE34\PSSE34EXAMPLE;" + os.environ['PATH'])
import psse34 #addition necessary for new version 34
import psspy
import pssarrays
import redirect
import dyntools
import pssplot
from dtw import dtw
#import pssexcel
import random
import copy
import math
import multiprocessing
import time
import sys
import csv
import numpy as np
import redirect
import silence
from silence import silence
from sdd import sdd
from sts import sts 
from random import randint
redirect.psse2py() 
#import matplotlib
_i=psspy.getdefaultint()
_f=psspy.getdefaultreal()
_s=psspy.getdefaultchar()
redirect.psse2py()
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import cosine as cosine
from scipy.spatial.distance import chebyshev as chebyshev
from scipy.spatial.distance import correlation as correlation
from fastdtw import fastdtw


## Calculate the Error at different accuracies
#   Where the benchmark is CLMs placed at a specific number of buses on specific loads, and ZIP everywhere else
#   The Accuracy Tests then randomly chose loads so that an x percentage are on the same loads as in the benchmark


## Initializations and File Loading ============================================================================================
# Time how long this takes
t0 = time.time()
# Load load bus numbers
load_bus_info = np.genfromtxt('RTSLoads.csv',delimiter=',') #first column is bus number, second column is list number to match with pms
total_loads = 13
print "Loads with CLM is %d out of %d total loads" %(total_loads,len(load_bus_info[:,1]))
#===========================================================================================================================


## Function: Create dynamic file ===========================================================================================
def Create_DYR(name, buses): #buses is matrix size(1/2 load buses,2) where the first column is the actual bus # and column 2 is the list bus #(where this bus appears in the 
							 #in the pms list. name is string to call the dynamic file
	# Load parameter list
	pms = np.genfromtxt('load_parameters.csv',delimiter=',') #randomized parameters for each load 
	#Open existing dynamic file`
	f = open('RTS96DYN.dyr','r') #read in file
	ZIP = f.read() #saves contents of file to variable
	filename = "dynamic_%s.dyr" %name
	f = open(filename,'w')
	f.write(ZIP + '\n')
	for i in range(len(buses[:,1])):
     									#type	 ID    type         12   IT 2  133  27  146  48  0, 0,  J (-1.25 or 0.8)    +1 +2 (0.04 or 0)   +3 (0.04 or 0)   +4    +5    +6 +7 +8  +9  +10   +11      +12    +13  +14 +15+16+17          +18          +19          +20           +21         +22         +23
		f.write(str(int(buses[i,0]))+", 'USRLOD', '1', 'CMLDBLU2', 12, 1, 2, 133, 27, 146, 48, 0, 0, -1, 0, 0.04, 0.04, 0.75, 0.08, 1, 1, 0, 0.9, 1.1, 0.00625, 1.025, 1.04, 30, 5, 0, 0," + " \n "
		#+18                +19              +20                 +21               +22         
		+str(pms[0,int(buses[i,1])])+","+str(pms[1,int(buses[i,1])])+","+str(pms[2,int(buses[i,1])])+","+str(pms[3,int(buses[i,1])])+","+str(pms[4,int(buses[i,1])])+ " \n "
		#+23    +24                +25           26   
		+"1,"+str(pms[5,int(buses[i,1])])+","+str(pms[6,int(buses[i,1])])+", 1,  " + " \n"
		#+27   +28    +29  +30  +31  
		+" 2, 0.452, 1, 0.548, 0, " + " \n "
		#+32 +33  +34  +35 +36
		+"2, -0.5, 1, 1.5, -1," + " \n "
		# +37 +38  +39  +40   +41   +42    +43    +44    +45  +46   +47                 +48              +49             +50   +51   +52  +53   +54                   +55           +56   
		+"3, 0.75, 0.04, 1.8, 0.12, 0.104, 0.095, 0.0021, 0.1, 0, "+str(pms[7,int(buses[i,1])])+","+str(pms[8,int(buses[i,1])])+","+str(pms[9,int(buses[i,1])])+", 0.1, 99999, 0.5, 0.02, "+str(pms[10,int(buses[i,1])])+","+str(pms[11,int(buses[i,1])])+", 0.1," + " \n "
		#+57 +58    +59  +60   +61   +62  +63   +64    +65 +66 +67                       +68               +69                +70             +71     +72              +73                  +74                +75         +76   
		+"3, 0.75, 0.03, 1.8, 0.19, 0.14, 0.2, 0.0026, 0.5, 2, "+str(pms[12,int(buses[i,1])])+","+str(pms[13,int(buses[i,1])])+","+str(pms[14,int(buses[i,1])])+","+str(pms[15,int(buses[i,1])])+", 0.05, "+str(pms[16,int(buses[i,1])])+","+str(pms[17,int(buses[i,1])])+","+str(pms[18,int(buses[i,1])])+","+str(pms[19,int(buses[i,1])])+", 0.05, " + " \n "
		#+77 +78    +79   +80  +81  +82   +83   +84    +85  +86        +87           +88                 +89                 +90              +91  +92                  +93                +94               +95      
		+"3, 0.75, 0.03, 1.8, 0.19, 0.14, 0.2, 0.0026, 0.1, 2, "+str(pms[20,int(buses[i,1])])+","+str(pms[21,int(buses[i,1])])+","+str(pms[22,int(buses[i,1])])+","+str(pms[23,int(buses[i,1])])+", 0.5, "+str(pms[24,int(buses[i,1])])+","+str(pms[25,int(buses[i,1])])+","+str(pms[26,int(buses[i,1])])+","+str(pms[27,int(buses[i,1])])+"," + " \n "
		#+96    +97           +98           +99   +100 +101 +102      +103           +104  +105  +106
		+"0.05, 0.0666, "+str(pms[28,int(buses[i,1])])+", 0.025, 0.1, 1, 0.98, "+str(pms[29,int(buses[i,1])])+", 0.1, 0.1, 0.0," + " \n "
		#+107  +108  +109 +110 +111  +112 +113 +114
		+"0.0, 1.0, 6.0, 2.0, 12.0, 3.2, 11.0, 2.5," + " \n "
		#+115     +116                    +117          +118 +119   +120                +121                +122                 +123             +124         +125  +126  +127           +128           +129             +130  +131           +132
		+"0.86, "+str(pms[30,int(buses[i,1])])+","+str(pms[31,int(buses[i,1])])+", 1.0, -3.3, "+str(pms[32,int(buses[i,1])])+","+str(pms[33,int(buses[i,1])])+","+str(pms[34,int(buses[i,1])])+","+str(pms[35,int(buses[i,1])])+","+str(pms[36,int(buses[i,1])])+", 0.7, 1.9, 0.025, "+str(pms[37,int(buses[i,1])])+","+str(pms[38,int(buses[i,1])])+", 0.1, 9999," + " \n "+str(pms[39,int(buses[i,1])])+" /  CMLDBLU2 " + "\n")
		i = i+1
	f.close()
#===========================================================================================================================



## Function: Run_SIM() =====================================================================================================
def Run_SIM(x,dyr_file,out_file): #inputs are strings\
	dyre = r"""C:\Users\psse\Desktop\Phylicia\Error and Accuracy Tracking Project Sp18\RTS96\%s""" %dyr_file
	out = r"""C:\Users\psse\Desktop\Phylicia\Error and Accuracy Tracking Project Sp18\RTS96\Channels\opt_%s.out""" %out_file
	print dyr_file
	ierr = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] #checking for errors
	output = StringIO.StringIO()
	with silence(output):
		ierr[0] = psspy.psseinit(200000) #need to have this high, otherwise there are not enough output channels
		ierr[1] = psspy.case(r"""C:\Users\psse\Desktop\Phylicia\Error and Accuracy Tracking Project Sp18\RTS96\RTS96DYN.sav""")
		ierr[2] = psspy.fdns([0,0,0,1,1,0,99,0])
		ierr[3] = psspy.cong(0)
		ierr[4] = psspy.conl(0,1,1,[0,0],[ 100.0,0.0,0.0, 100.0])
		ierr[5] = psspy.conl(0,1,2,[0,0],[ 100.0,0.0,0.0, 100.0])
		ierr[6] = psspy.conl(0,1,3,[0,0],[ 100.0,0.0,0.0, 100.0])
		ierr[7] = psspy.ordr(0)
		ierr[8] = psspy.fact()
		ierr[9] = psspy.tysl(0)
		ierr[10] = psspy.dyre_new([1,1,1,1],dyre,"","","") 
		ierr[11] = psspy.chsb(0,1,[-1,-1,-1,1,13,0]) #record voltage
		ierr[12] = psspy.chsb(0,1,[-1,-1,-1,1,12,0]) #record frequency
		ierr[13] = psspy.chsb(0,1,[-1,-1,-1,1,1,0]) #angle
		ierr[14] = psspy.chsb(0,1,[-1,-1,-1,1,16,0]) #line P & Q
		ierr[15] = psspy.strt_2([0,0],out)
		ierr[16] = psspy.run(0, 0.1,1,1,0)
		#ierr[17] = psspy.dist_branch_fault(217,218,r"""1""",1, 230.0,[0.0,-0.2E+10]) #Line Fault, NaN (network no good)
		#ierr[17] = psspy.dist_bus_fault(211,1, 230.0,[0.0,-0.2E+10]) #bus Fault, NaN (network no good)
		#a = int(x[0])
		#b = int(x[1])
		#ierr[17] = psspy.branch_chng_3(a,b,r"""1""",[0,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],"") #Line Outage
		x = int(x)
		print "before machine change"
		ierr[17] = psspy.machine_chng_2(x,r"""1""",[0,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f]) #Generator Outage
		print "after machine change"
		ierr[18] = psspy.change_channel_out_file(out)
		ierr[19] = psspy.run(0, 0.5,1,1,0) #this was 10
		psspy.dist_clear_fault(1)
		psspy.change_channel_out_file(out)
		psspy.run(1, 10.0,1,1,0)
		ierr[20] = psspy.delete_all_plot_channels()
		print "completed simulation"
	
	print ierr 
	run_output = output.getvalue()
	
	current_error = 0
	
	if "Network not converged" in run_output:
		print "Network not converged" #need to have something in if statement otherwise you get an indentation error
		result = 0 #this will go out to a if condition to rerun the program with a different selection of buses at this accuracy
		current_error = 1
		#raise SystemExit #this will quit the program
	elif "NaN" in run_output:
		print "NaN, network is no good"
		result = 0 #this will go out to a if condition to rerun the program with a different selection of buses at this accuracy
		current_error = 1
		#raise SystemExit #this will quit the program
	if current_error == 0 and "INITIAL CONDITIONS CHECK O.K." in run_output:
		print "continuing with study..."
	
        #Gather the data and output to excel
		data = dyntools.CHNF(out) #getting data from channel.out file
		d,e,z=data.get_data() #gathering data from data in dictionary format
	
        #Concatenate data so all data from one simulation is in one file
		c = 1 #must start at 1, not zero
        #Save Frequency and Voltage
		while c < 726: 
			if c < 100: #Record Angle
				v=z[c]
				new_v = ", ".join(str(i) for i in v) #this removes the brackets at the beginning and end of the list so can be processed in matlab 
				a = np.matrix(new_v) #make it into a matrix
				if c ==1:
					ang_all = np.copy(a)
				else: 
					ang_all = np.concatenate((ang_all,a),axis=0) #changed to concatenate vertically to test them all individually
			if c > 99 and c < 173: #Record Frequency
				v=z[c]
				new_v = ", ".join(str(i) for i in v) #this removes the brackets at the beginning and end of the list so can be processed in matlab 
				f = np.matrix(new_v) #make it into a matrix
				if c ==100:
					f_all = np.copy(f)
				else: 
					f_all = np.concatenate((f_all,f),axis=0) #changed to concatenate vertically to test them all individually
			if c > 172 and c < 246: #Record voltage magnitude
				v=z[c]
				new_v = ", ".join(str(i) for i in v) #this removes the brackets at the beginning and end of the list so can be processed in matlab 
				f = np.matrix(new_v) #make it into a matrix
				if c == 173:
					all = np.copy(f)
				else:
					all = np.concatenate((all,f),axis=0) #changed to concatenate vertically to test them all individually
			if c > 245 and c < 726: #Record P and Q
				if float(c/2) == int(c/2): #P , even numbers
					v=z[c]
					new_v = ", ".join(str(i) for i in v) #this removes the brackets at the beginning and end of the list so can be processed in matlab 
					f = np.matrix(new_v) #make it into a matrix
					if c == 246:
						P_all = np.copy(f)
					else:
						P_all = np.concatenate((P_all,f),axis=0) #changed to concatenate vertically to test them all individually
				else: #Q, odd numbers
					v=z[c]
					new_v = ", ".join(str(i) for i in v) #this removes the brackets at the beginning and end of the list so can be processed in matlab 
					f = np.matrix(new_v) #make it into a matrix
					if c == 247:
						Q_all = np.copy(f)
					else:
						Q_all = np.concatenate((Q_all,f),axis=0) #changed to concatenate vertically to test them all individually
			c = c+1
		result = [all, f_all, ang_all, P_all, Q_all] #0 is voltage, 1 is frequency
	return result
#===================================================================================================================================================	

## Define Function to check that new load is chosen for solution ==========================================================================================
def check_doubles(x,soln,bus_list):
	for p in range(len(soln[:,0])): #check if x is already in the solution
		if bus_list[x,1] == soln[p,1]: #if in the solution change x
			x = np.random.randint(len(bus_list)) #chose integer within the length range of all benchmark buses, and choose one
			x = check_doubles(x,soln,bus_list) #recursion depth issue!
	return x
##=================================================================================================================================================
	


## Define Check Accuracy =====================================================================================================================================
def check_acc(soln,buses_bench):
	bus_check = np.zeros(len(buses_bench))
	for i in range(len(buses_bench)):
		for j in range(len(soln)):
			if buses_bench[i,0] == soln[j,0]:
				bus_check[i] = 1
	perc_correct = sum(bus_check)/len(buses_bench)*100
	return perc_correct
##==============================================================================================================================================================

## Define Bus Lists at Different Accuracy Levels Randomly =========================================================================================================
def bus_list(accuracy,total_loads,buses_bench,non_buses_bench): #accuracy is the accuracy percent desired
	accuracy = accuracy/100 #turn into a percentage
	soln = np.zeros((total_loads,2)) #set initial solution #first column is bus number, second column is list number to match with pms
	bench_total = int(accuracy*total_loads) #number of buses in solution chosen from the benchmark is accuracy*total_loads (total_loads is the number of CLM in the system)
	non_bench_total = int((1-accuracy)*total_loads)
	if bench_total+non_bench_total < total_loads: #if the total CLMs does not equal 365
		non_bench_total += 1
	for i in range(bench_total):
		option = np.random.randint(len(buses_bench)) #chose integer within the length range of all benchmark buses, and choose one
		option = check_doubles(option,soln,buses_bench) #make sure the option is not already in the solution
		soln[i,:] = buses_bench[option,:] #save both the list number and bus number in the same order
		
	for i in range(bench_total,(bench_total+non_bench_total)):
		option = np.random.randint(len(non_buses_bench)) #chose integer within the length range of all on benchmark buses, and choose one
		option = int(check_doubles(option,soln,non_buses_bench)) #make sure the option is not already in the solution
		soln[i,:] = non_buses_bench[option,:]
	
	return soln		
##===============================================================================================================================================================

## Run Accuracy Bus List Generator ============================================================================================================================
#Create bus lists of specific accuracy
def accuracy_bus_list(total_loads,buses_bench,non_buses_bench):
    accuracy = [0,1/13*100,2/13*100,3/13*100,4/13*100,5/13*100,6/13*100,7/13*100,8/13*100,9/13*100,10/13*100,11/13*100,12/13*100,100] #Put in Benchmark to verify everything is working correctly, and to calculate base noise error
    bus_list_actual_nums = np.zeros((total_loads,(len(accuracy)*2)))
    for i in range(len(accuracy)):
        soln = bus_list(accuracy[i],total_loads,buses_bench,non_buses_bench)
        start = 2*i
        end = (2*i)+2
        bus_list_actual_nums[:,start:end] = soln
    filename = "Bus_Lists.csv" 
    np.savetxt(filename, bus_list_actual_nums, delimiter=",")
    return bus_list_actual_nums
##===========================================================================================================================================================


	
## Run Error and Accuracy Check ================================================================================================================================

## Create CSV Files for Cases
#load benchmark nd non benchmark buses lists
buses_bench = np.genfromtxt('Buses_Bench_RTS96.csv',delimiter=',') #first column is bus number, second column is list number to match with 
non_buses_bench = np.genfromtxt('Non_Bench_Buses_RTS96.csv', delimiter=",")

#Run Accuracy Bus List Generator
bus_list_actual_nums = accuracy_bus_list(total_loads,buses_bench,non_buses_bench)

#Run the different accuracy levels and record the error
branches = np.genfromtxt('RTS96_branches.csv',delimiter=',')
accuracy = [0,1/13*100,2/13*100,3/13*100,4/13*100,5/13*100,6/13*100,7/13*100,8/13*100,9/13*100,10/13*100,11/13*100,12/13*100,100]
V_ED = np.zeros((len(accuracy),100))
V_DTW = np.zeros((len(accuracy),100))
V_COS = np.zeros((len(accuracy),100))
V_CHE = np.zeros((len(accuracy),100))
V_COR = np.zeros((len(accuracy),100))
V_MH = np.zeros((len(accuracy),100))
F_ED = np.zeros((len(accuracy),100))
F_DTW = np.zeros((len(accuracy),100))
F_COS = np.zeros((len(accuracy),100))
F_CHE = np.zeros((len(accuracy),100))
F_COR = np.zeros((len(accuracy),100))
F_MH = np.zeros((len(accuracy),100))
Ang_ED = np.zeros((len(accuracy),100))
Ang_DTW = np.zeros((len(accuracy),100))
Ang_COS = np.zeros((len(accuracy),100))
Ang_CHE = np.zeros((len(accuracy),100))
Ang_COR = np.zeros((len(accuracy),100))
Ang_MH = np.zeros((len(accuracy),100))
P_ED = np.zeros((len(accuracy),100))
P_DTW = np.zeros((len(accuracy),100))
P_COS = np.zeros((len(accuracy),100))
P_CHE = np.zeros((len(accuracy),100))
P_COR = np.zeros((len(accuracy),100))
P_MH = np.zeros((len(accuracy),100))
Q_ED = np.zeros((len(accuracy),100))
Q_DTW = np.zeros((len(accuracy),100))
Q_COS = np.zeros((len(accuracy),100))
Q_CHE = np.zeros((len(accuracy),100))
Q_COR = np.zeros((len(accuracy),100))
Q_MH = np.zeros((len(accuracy),100))
errors = 0
for h in range(100): #repeat test 100 times
    k = 0
    while k < (len(accuracy)):
		start = 2*k
		end = (2*k)+2
		soln = bus_list_actual_nums[:,start:end]
		#Run TEST and BENCH
		#Run Test
		name = "Accuracy_%d" %accuracy[k]
		Create_DYR(name,soln)
		dynamic = "dynamic_%s.dyr" %name
		out = "Accuracy_bus_%s.out" %name
		#Choose random generator to have outage
		gen_buses = [101,102,107,113,114,115,116,118,121,122,123,201,202,207,213,214,215,216,218,221,222,223,301,302,307,313,314,315,316,318,321,322,323]
		gen_bus = gen_buses[randint(0,len(gen_buses)-1)]
		print "Before simulation"
		result_test = Run_SIM(gen_bus,dynamic,out) #array of one row containing all voltage output CHANGED to an array of (rows = # of voltages, columns = # of time instances)
		print "After simulation"
		#branch = branches[randint(0,len(branches)-1)]
		#result_test = Run_SIM(branch,dynamic,out) #array of one row containing all voltage output CHANGED to an array of (rows = # of voltages, columns = # of time instances)
		#Run Bench
		#Run Benchmark
		name = "benchmark"
		Create_DYR(name,buses_bench) 
		dynamic = "dynamic_%s.dyr" %name
		out = "Benchmark_%s.out" %name
		result_bench = Run_SIM(branch,dynamic,out) #array of one row containing all voltage output CHANGED to an array of (rows = # of voltages, columns = # of time instances)
		if result_test != 0 and result_bench != 0:
			[V_TEST,F_TEST,Ang_TEST,P_TEST,Q_TEST] = result_test
			[V_BENCH,F_BENCH,Ang_BENCH,P_BENCH,Q_BENCH] = result_bench
			#Voltage
			V_ED[k,h] = np.sum(np.sqrt(np.sum(np.square(V_BENCH-V_TEST), axis = 1))*1000) #times by 1000 so the objnow isn't so small, sums horizontally
			V_MH[k,h] = np.sum(np.sum(np.absolute(V_BENCH-V_TEST), axis = 1)*1000) #manhattan distance
            #dtw and sts distance
			V_len = len(V_BENCH[:,1])
			#print V_len
			#dist_all = np.zeros(V_len)
			d_all_DTW = np.zeros(V_len)
			d_all_COS = np.zeros(V_len)
			d_all_CHE = np.zeros(V_len)
			d_all_COR = np.zeros(V_len)
			#d_all_SDD = np.zeros(V_len)
			for j in range(V_len):
				x = V_BENCH[j]
				y = V_TEST[j]
				#dist_all[j] = sts((1/120),x,y)
				#print dist_all[j] #still zero
				d_all_DTW[j], path = fastdtw(x, y, dist=euclidean)
				d_all_COS[j] = cosine(x,y)
				d_all_CHE[j] = chebyshev(x,y)
				d_all_COR[j] = correlation(x,y)
				#d_all_SDD[j] = sdd(x,y)
			#V_STS[k,h] = np.sum(dist_all)
			V_DTW[k,h] = np.sum(d_all_DTW)
			V_COS[k,h] = np.sum(d_all_COS)
			V_CHE[k,h] = np.sum(d_all_CHE)
			V_COR[k,h] = np.sum(d_all_COR)
			#V_SDD[k,h] = np.sum(d_all_SDD)
    
            #Frequency
			F_ED[k,h] = np.sum(np.sqrt(np.sum(np.square(F_BENCH-F_TEST), axis = 1))*1000) #times by 1000 so the objnow isn't so small, sums horizontally
			F_MH[k,h] = np.sum(np.sum(np.absolute(F_BENCH-F_TEST), axis = 1)*1000) #manhattan distance
            #dtw and sts distance
			F_len=len(F_BENCH[:,1])
			#dist_all = np.zeros(F_len)
			d_all_DTW = np.zeros(F_len)
			d_all_COS = np.zeros(F_len)
			d_all_CHE = np.zeros(F_len)
			d_all_COR = np.zeros(F_len)
			#d_all_SDD = np.zeros(F_len)
			for j in range(F_len):
				x = F_BENCH[j]
				y = F_TEST[j]
				#dist_all[j] = sts((1/120),x,y)
				d_all_DTW[j], path = fastdtw(x, y, dist=euclidean)
				d_all_COS[j] = cosine(x,y)
				d_all_CHE[j] = chebyshev(x,y)
				d_all_COR[j] = correlation(x,y)
				#d_all_SDD[j] = sdd(x,y)
			#F_STS[k,h] = np.sum(dist_all)
			F_DTW[k,h] = np.sum(d_all_DTW)
			F_COS[k,h] = np.sum(d_all_COS)
			F_CHE[k,h] = np.sum(d_all_CHE)
			F_COR[k,h] = np.sum(d_all_COR)
			#F_SDD[k,h] = np.sum(d_all_SDD)

			#Angle
			Ang_ED[k,h] = np.sum(np.sqrt(np.sum(np.square(Ang_BENCH-Ang_TEST), axis = 1))*1000) #times by 1000 so the objnow isn't so small, sums horizontally
			Ang_MH[k,h] = np.sum(np.sum(np.absolute(Ang_BENCH-Ang_TEST), axis = 1)*1000) #manhattan distance
            #dtw and sts distance
			Ang_len=len(Ang_BENCH[:,1])
			d_all_DTW = np.zeros(Ang_len)
			d_all_COS = np.zeros(Ang_len)
			d_all_CHE = np.zeros(Ang_len)
			d_all_COR = np.zeros(Ang_len)
			for j in range(Ang_len):
				x = Ang_BENCH[j]
				y = Ang_TEST[j]
				d_all_DTW[j], path = fastdtw(x, y, dist=euclidean)
				d_all_COS[j] = cosine(x,y)
				d_all_CHE[j] = chebyshev(x,y)
				d_all_COR[j] = correlation(x,y)
			Ang_DTW[k,h] = np.sum(d_all_DTW)
			Ang_COS[k,h] = np.sum(d_all_COS)
			Ang_CHE[k,h] = np.sum(d_all_CHE)
			Ang_COR[k,h] = np.sum(d_all_COR)
			
			#P
			P_ED[k,h] = np.sum(np.sqrt(np.sum(np.square(P_BENCH-P_TEST), axis = 1))*1000) #times by 1000 so the objnow isn't so small, sums horizontally
			P_MH[k,h] = np.sum(np.sum(np.absolute(P_BENCH-P_TEST), axis = 1)*1000) #manhattan distance
            #dtw and sts distance
			P_len=len(P_BENCH[:,1])
			d_all_DTW = np.zeros(P_len)
			d_all_COS = np.zeros(P_len)
			d_all_CHE = np.zeros(P_len)
			d_all_COR = np.zeros(P_len)
			for j in range(P_len):
				x = P_BENCH[j]
				y = P_TEST[j]
				d_all_DTW[j], path = fastdtw(x, y, dist=euclidean)
				d_all_COS[j] = cosine(x,y)
				d_all_CHE[j] = chebyshev(x,y)
				d_all_COR[j] = correlation(x,y)
			P_DTW[k,h] = np.sum(d_all_DTW)
			P_COS[k,h] = np.sum(d_all_COS)
			P_CHE[k,h] = np.sum(d_all_CHE)
			P_COR[k,h] = np.sum(d_all_COR)
			
			#Q
			Q_ED[k,h] = np.sum(np.sqrt(np.sum(np.square(Q_BENCH-Q_TEST), axis = 1))*1000) #times by 1000 so the objnow isn't so small, sums horizontally
			Q_MH[k,h] = np.sum(np.sum(np.absolute(Q_BENCH-Q_TEST), axis = 1)*1000) #manhattan distance
            #dtw and sts distance
			Q_len=len(Q_BENCH[:,1])
			d_all_DTW = np.zeros(Q_len)
			d_all_COS = np.zeros(Q_len)
			d_all_CHE = np.zeros(Q_len)
			d_all_COR = np.zeros(Q_len)
			for j in range(Q_len):
				x = Q_BENCH[j]
				y = Q_TEST[j]
				d_all_DTW[j], path = fastdtw(x, y, dist=euclidean)
				d_all_COS[j] = cosine(x,y)
				d_all_CHE[j] = chebyshev(x,y)
				d_all_COR[j] = correlation(x,y)
			Q_DTW[k,h] = np.sum(d_all_DTW)
			Q_COS[k,h] = np.sum(d_all_COS)
			Q_CHE[k,h] = np.sum(d_all_CHE)
			Q_COR[k,h] = np.sum(d_all_COR)
			
            #Save solutions as you go
			#volt
			np.savetxt("V_ED_GO_10s_clear_Error.csv",V_ED,delimiter=",")
			np.savetxt("V_MH_GO_10s_clear_Error.csv",V_MH,delimiter=",")
			np.savetxt("V_DTW_GO_10s_clear_Error.csv",V_DTW,delimiter=",")
			np.savetxt("V_COS_GO_10s_clear_Error.csv",V_COS,delimiter=",")
			np.savetxt("V_CHE_GO_10s_clear_Error.csv",V_CHE,delimiter=",")
			np.savetxt("V_COR_GO_10s_clear_Error.csv",V_COR,delimiter=",")
			#F
			np.savetxt("F_ED_GO_10s_clear_Error.csv",F_ED,delimiter=",")
			np.savetxt("F_MH_GO_10s_clear_Error.csv",F_MH,delimiter=",")
			np.savetxt("F_DTW_GO_10s_clear_Error.csv",F_DTW,delimiter=",")
			np.savetxt("F_COS_GO_10s_clear_Error.csv",F_COS,delimiter=",")	
			np.savetxt("F_CHE_GO_10s_clear_Error.csv",F_CHE,delimiter=",")	
			np.savetxt("F_COR_GO_10s_clear_Error.csv",F_COR,delimiter=",")
			#Ang			
			np.savetxt("Ang_ED_GO_10s_clear_Error.csv",Ang_ED,delimiter=",")
			np.savetxt("Ang_MH_GO_10s_clear_Error.csv",Ang_MH,delimiter=",")
			np.savetxt("Ang_DTW_GO_10s_clear_Error.csv",Ang_DTW,delimiter=",")
			np.savetxt("Ang_COS_GO_10s_clear_Error.csv",Ang_COS,delimiter=",")	
			np.savetxt("Ang_CHE_GO_10s_clear_Error.csv",Ang_CHE,delimiter=",")	
			np.savetxt("Ang_COR_GO_10s_clear_Error.csv",Ang_COR,delimiter=",")
			#P Flows
			np.savetxt("P_ED_GO_10s_clear_Error.csv",P_ED,delimiter=",")
			np.savetxt("P_MH_GO_10s_clear_Error.csv",P_MH,delimiter=",")
			np.savetxt("P_DTW_GO_10s_clear_Error.csv",P_DTW,delimiter=",")
			np.savetxt("P_COS_GO_10s_clear_Error.csv",P_COS,delimiter=",")	
			np.savetxt("P_CHE_GO_10s_clear_Error.csv",P_CHE,delimiter=",")	
			np.savetxt("P_COR_GO_10s_clear_Error.csv",P_COR,delimiter=",")	
			#Q Flows
			np.savetxt("Q_ED_GO_10s_clear_Error.csv",Q_ED,delimiter=",")
			np.savetxt("Q_MH_GO_10s_clear_Error.csv",Q_MH,delimiter=",")
			np.savetxt("Q_DTW_GO_10s_clear_Error.csv",Q_DTW,delimiter=",")
			np.savetxt("Q_COS_GO_10s_clear_Error.csv",Q_COS,delimiter=",")	
			np.savetxt("Q_CHE_GO_10s_clear_Error.csv",Q_CHE,delimiter=",")	
			np.savetxt("Q_COR_GO_10s_clear_Error.csv",Q_COR,delimiter=",")	
		else: 
			bus_list_actual_nums = accuracy_bus_list(total_loads,buses_bench,non_buses_bench) #change the bus list is the network didn't converge
			k = k-1
			errors += 1
		k += 1
    bus_list_actual_nums = accuracy_bus_list(total_loads,buses_bench,non_buses_bench) #change the bus list for every new data set
    
#OutputTime
t1 = time.time()
total_time = t1-t0
print "The total time is %f" %total_time #7 minutes to do 13 simulations. 
print "the number of errors is %d" %errors
#np.savetxt("Errors_10s_randomized",errors)   
