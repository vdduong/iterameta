### assignment assessment
### routines for assignment assessment, 3 routines:
### number of peaks found in experimental peak list in comparison to theoretical pattern
### number of peaks found in assignment cluster (convergence)
### score of fractional Hausdorff distance
### Viet Dung Duong, Ruben Staub, Torsten Herrmann

import math
from assignment_clustering import *
from iterative_matching import *
from class_metabolite import *
from spectrum import *

dict_assignment = {}
tol_H = 0.1
threshold = 0.05

def peaks_found(dict_assignment):
	'''return the number of peaks found in experimental peak list in comparison 
	to theoretical pattern, duplicate removed'''
	nb_peaks_found = len(dict_assignment.keys())
	set_assignments = set()
	for key_tocsy in dict_assignment.keys():
		if len(dict_assignment[key_tocsy]) == 0: nb_peaks_found-=1
		else:
			for key_exp in dict_assignment[key_tocsy].keys():
				set_assignments.add(key_exp)
	
	return min(nb_peaks_found, len(set_assignments))

def peaks_clustered(len_pattern, dict_patterns):
	for key_pattern in dict_patterns.keys():
		set_assignments = set()
		for assignment_ in dict_patterns[key_pattern]:
			set_assignments.add(assignment_[0])
		if float(len(set_assignments)) < 0.6*float(len_pattern):
			del dict_patterns[key_pattern]
	return dict_patterns
	
def hausdorffDistance(fractionHausdorff, metabolite, dict_assignment, peak_list):
	'''compute the Hausdorff distance between the TOCSY pattern and the experimental
	peak list
	return the score, given the pattern fraction used'''
	list_distance_h = []
	tocsy_pattern = metabolite.tocsy_pattern
	for key_tocsy in dict_assignment.keys():
		min_distance = float('inf')
		cs_tocsy_1 = float(tocsy_pattern[key_tocsy][0])
		cs_tocsy_2 = float(tocsy_pattern[key_tocsy][1])
		if len(dict_assignment[key_tocsy]) > 0:
			for key_exp in dict_assignment[key_tocsy].keys():
				cs_exp_1 = float(peak_list[key_exp][0])
				cs_exp_2 = float(peak_list[key_exp][1])
				distance_h = math.sqrt((cs_tocsy_1-cs_exp_1)**2.0 + (cs_tocsy_2-cs_exp_2)**2.0)
				if distance_h < min_distance: min_distance=distance_h
			list_distance_h.append(min_distance)
		else:
			for key_exp in peak_list.keys():
				cs_exp_1 = float(peak_list[key_exp][0])
				cs_exp_2 = float(peak_list[key_exp][1])
				distance_h = (cs_tocsy_1-cs_exp_1)**2.0+(cs_tocsy_2-cs_exp_2)**2.0
				if distance_h < min_distance: min_distance=distance_h
			list_distance_h.append(math.sqrt(min_distance))
	list_distance_h.sort(reverse=False)
	nb_peaks_hausdorff = int(fractionHausdorff*len(list_distance_h))
	return list_distance_h[nb_peaks_hausdorff-1]


if __name__ == '__main__':
	path_local = str(os.getcwd())
	peak_list, volume_list = spectrum(path_local + '/43-sn10.peaks')
	name_file = path_local + '/database/bmse000408.str'
	c = Metabolite(name_file)
	shift = c.initial_shift(name_file)
	c = Metabolite(name_file, shift)
	c.tocsy_pattern = c.tocsy()
	tocsy_initial = c.tocsy_pattern
	dict_assignment = assignment(c, peak_list)
	dict_patterns = pattern_production(dict_assignment, peak_list, c, tocsy_initial, threshold)
	dict_patterns = peaks_clustered(len(tocsy_initial), dict_patterns)
	c = update_shift(c, peak_list, dict_assignment, dict_patterns)
	c.tocsy_pattern = c.tocsy()
	print 'peaks found : ', peaks_found(dict_assignment)
	fractionHausdorff = 0.7
	print hausdorffDistance(fractionHausdorff, c, dict_assignment, peak_list)
					
