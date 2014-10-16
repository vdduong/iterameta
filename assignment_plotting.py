# plotting assignments
# Viet Dung Duong, Ruben Staub, Torsten Herrmann

import matplotlib.pyplot as plt
import numpy as np 
from spectrum import *
from class_metabolite import *
from assignment_clustering import *
from assignment_assessment import *
from iterative_matching import *
import math
import re
import glob
import os
from pylab import *

tol_H = 0.1
threshold = 0.05

def show_result(assignment_pattern, tocsy_pattern, peak_list, dict_assignment):
	fig=plt.figure()
	ax_view=fig.add_subplot(1,1,1)
	
	#ax_view.scatter([float(tocsy_pattern[item[0]][0]) for item in assignment_pattern], \
	#	[float(tocsy_pattern[item[0]][1]) for item in assignment_pattern], color='blue', s=20.)
	ax_view.scatter([float(tocsy_pattern[item][0]) for item in tocsy_pattern.keys()], \
		[float(tocsy_pattern[item][1]) for item in tocsy_pattern.keys()], color='blue', s= 20.)
	
	#ax_view.scatter([float(peak_list[item[1]][0]) for item in assignment_pattern], \
	#	[float(peak_list[item[1]][1]) for item in assignment_pattern], color='red')
	ax_view.scatter([float(peak_list[item][0]) for item in peak_list.keys()], \
		[float(peak_list[item][1]) for item in peak_list.keys()], color='red', s= 60.)
	set_exp = set()
	for peak_tocsy in dict_assignment.keys():
			for peak_exp in dict_assignment[peak_tocsy].keys(): set_exp.add(peak_exp)
	#ax_view.scatter([float(peak_list[item][0]) for item in set_exp], \
	#	[float(peak_list[item][1]) for item in set_exp], color='red', s=20.)
	ax_view.invert_yaxis()
	ax_view.invert_xaxis()
	ax_view.set_xlabel('F1 (ppm)')
	ax_view.set_ylabel('F2 (ppm)')
	ax_view.set_title('Threonine')
	
	#ax_view.annotate('1', xy = (4.121, 3.652), xytext = (4.0, 3.7), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('A', xy = (4.241, 3.573), xytext = (4.0, 3.3), arrowprops=dict(facecolor='black', width = 0.01))
	
	#peak_233 = (float(peak_list[233][0]), float(peak_list[233][1]))
	#peak_227 = (float(peak_list[227][0]), float(peak_list[227][1]))
	
	#ax_view.annotate('1', xy = peak_227, xytext = (4.2, 4.2), arrowprops=dict(facecolor='black', width = 0.01)) # for hippuric
	#ax_view.annotate('2', xy = peak_233, xytext = (3.7, 3.7), arrowprops=dict(facecolor='black', width = 0.01)) # for hippuric
	
	#ax_view.annotate('1#', xy = (1.258,1.298), xytext = (1.2, 1.22), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('2#', xy = (1.258,3.653), xytext = (1.2, 3.6), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('3#', xy = (1.258,4.323), xytext = (1.2, 4.3), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('1', xy = (1.398,1.399), xytext = (1.5, 1.5), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('2', xy = (1.398,3.652), xytext = (1.5, 3.8), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('3', xy = (1.398,4.323), xytext = (1.5, 4.5), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('A', xy = (1.318,1.318), xytext = (1.5, 1.2), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('B', xy = (1.318,3.573), xytext = (1.5, 3.42), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('C', xy = (1.318,4.241), xytext = (1.5, 4.12), arrowprops=dict(facecolor='black', width = 0.01))
	fig.show()

if __name__ == '__main__':
	path_local = os.getcwd()
	#peak_list, volume_list = spectrum(path_local + '/hippuric.peaks')
	peak_list, volume_list = spectrum(path_local + '/threonine_2.peaks')
	#name_file = path_local+'/database/bmse000408.str'
	name_file = path_local+'/database/bmse000049.str'
	#peak_list, volume_list = spectrum(path_local + "/43-sn10.peaks")
	dict_assignment = {}
	c = Metabolite(name_file)
	shift = c.initial_shift(name_file)
	c = Metabolite(name_file, shift)
	c.tocsy_pattern = c.tocsy()
	tocsy_initial = c.tocsy_pattern
	len_tocsy_initial = len(tocsy_initial)
	dict_assignment = assignment(c, peak_list)
	dict_patterns = pattern_production(dict_assignment, peak_list, c, len_tocsy_initial, threshold)
	dict_best_patterns = best_assignment(dict_patterns, dict_assignment)
	#assignment_showing(dict_best_patterns, dict_assignment, peak_list, c, len_tocsy_initial)
	pattern = dict_best_patterns[1]
	print dict_assignment
	print pattern
	tocsy_pattern = c.tocsy_pattern
	show_result(pattern, tocsy_pattern, peak_list, dict_assignment)
	show()
