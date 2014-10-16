### ITERATIVE MATCHING
### includes :
### - a function listing all the BMRB database
### - a function computing the assignments between the theoretical pattern of the 
### metabolite and the experimental peak list
### - a function updating the chemical shifts of the metabolite from the previous assignment
### - a function repeating the previous steps for 7 times, returning the final assignment,
### the corresponding Hausdorff distance, the total and local flow.
### - a function writing out the output into a xml file to be read by TopSpin
### Viet Dung Duong, Ruben Staub, Torsten Herrmann

from spectrum import *
from class_metabolite import *
from assignment_clustering import *
from assignment_assessment import *
from assignment_plotting import *
import math
import re
import glob
import os

tol_H = 0.1
threshold= 0.04

def database_bmrb() :
	"""
	return the list of all files contained in the folder database
	"""
	path_local = os.getcwd()
	path_final = str(path_local) + '/database/*.str'
	list_file = glob.glob(path_final)
	return list_file

def assignment(metabolite, peak_list):
	"""
	return the dictionary of available assignments
	a dictionary in Python is the equivalent of a hash table
	we use a dictionary of dictionaries to note the assignments and their corresponding probabilities
	"""
	global tol_H
	dict_assignment = dict()
	for key_tocsy in metabolite.tocsy_pattern.keys():
		item_tocsy = metabolite.tocsy_pattern[key_tocsy]
		cs_theor_1 = float(item_tocsy[0])
		cs_theor_2 = float(item_tocsy[1])
		dict_assignment[key_tocsy] = {}
		for key_peak_list in peak_list.keys():
			item_peak_list = peak_list[key_peak_list]
			cs_exp_1 = float(item_peak_list[0])
			cs_exp_2 = float(item_peak_list[1])
			if abs(cs_exp_1 - cs_theor_1) <= tol_H and \
				abs(cs_exp_2 - cs_theor_2) <= tol_H :
				proba = math.exp(-0.5*(cs_exp_1 - cs_theor_1)**2/ tol_H**2)*\
						math.exp(-0.5*(cs_exp_2 - cs_theor_2)**2/ tol_H**2)
				proba = float("%.4f"%proba)
				dict_assignment[key_tocsy][key_peak_list] = proba
	return dict_assignment

def update_shift(metabolite, peak_list, dict_assignment, dict_patterns):
	'''return the updated shifts after the assignment clustering'''
	global tol_H
	shift_list = list()
	for key in metabolite.shift.keys():
		if 'H' in key:
			item = metabolite.shift[key]
			shift_list.append(item)
	shift_list = remove_duplicate(shift_list)
	shift_dict = {}
	for item_shift_list in shift_list:
		shift_dict[item_shift_list] = 0.0
	for key_shift_dict in shift_dict.keys():
		sum_ponderation = 0.0
		time_ponderation = 0.0
		for key_pattern in dict_patterns.keys():
			pattern = dict_patterns[key_pattern] # pattern is actually a list
			for pair in pattern:
				peak_tocsy = metabolite.tocsy_pattern[pair[0]]
				peak_exp = peak_list[pair[1]]
				if key_shift_dict in peak_tocsy[0]:
					sum_ponderation+=float(peak_exp[0])*dict_assignment[pair[0]][pair[1]]
					time_ponderation+=dict_assignment[pair[0]][pair[1]]
				if key_shift_dict in peak_tocsy[1]:
					sum_ponderation+=float(peak_exp[1])*dict_assignment[pair[0]][pair[1]]
					time_ponderation+=dict_assignment[pair[0]][pair[1]]
		
		if time_ponderation > 0.0: 
			sum_ponderation = sum_ponderation/time_ponderation
			shift_dict[key_shift_dict]='%.4f'%(sum_ponderation)
		else:
			shift_dict[key_shift_dict]='%.4f'%(sum_ponderation)
	
	for key_shift in metabolite.shift.keys():
		item_shift = metabolite.shift[key_shift]
		for key_shift_dict in shift_dict:
			if key_shift_dict == item_shift and float(shift_dict[key_shift_dict]) != 0.000:
				metabolite.shift[key_shift] = shift_dict[key_shift_dict]
	metabolite.tocsy_pattern = metabolite.tocsy()
	return metabolite

def assignment_showing(dict_best_patterns, dict_assignment, peak_list, metabolite, len_tocsy_initial, volume_list):
	tocsy_pattern = metabolite.tocsy_pattern
	max_pattern = 0
	for key_pattern in dict_best_patterns.keys():
		if max_pattern < len(dict_best_patterns[key_pattern]): max_pattern = len(dict_best_patterns[key_pattern])
		else:pass
	nb_peaks_found = min(peaks_found(dict_assignment), max_pattern)
	
	print metabolite.name, 'was found with %i peaks over %i'%(nb_peaks_found, len_tocsy_initial)
	
	print '_THEORETICAL__EXPERIMENTAL__PROBA'
	mean_proba, time_proba = 0.0, 0.0
	if len(dict_best_patterns.keys()) == 1:
		for key_pattern in dict_best_patterns.keys():
			#print 'pattern '%key_pattern
			pattern = dict_best_patterns[key_pattern]
			for item in pattern:
				key_tocsy = item[0]
				key_exp = item[1]
				try:
					print tocsy_pattern[key_tocsy][0], tocsy_pattern[key_tocsy][1],'|',\
					peak_list[key_exp][0], peak_list[key_exp][1], '|', dict_assignment[key_tocsy][key_exp], volume_list[key_exp]
					mean_proba+=float(dict_assignment[key_tocsy][key_exp])
					time_proba+=1
				except KeyError:
					pass
			print 'Mean matching probability: %.4f'%(mean_proba/time_proba)
			print '________'
		#show_result(pattern, tocsy_pattern, peak_list)
	else:
		for key_pattern in dict_best_patterns.keys():
			print 'Assignment possibility number: %i'%key_pattern
			pattern = dict_best_patterns[key_pattern]
			for item in pattern:
				key_tocsy = item[0]
				key_exp = item[1]
				try:
					print tocsy_pattern[key_tocsy][0], tocsy_pattern[key_tocsy][1],'|',\
					peak_list[key_exp][0], peak_list[key_exp][1], '|', dict_assignment[key_tocsy][key_exp], volume_list[key_exp]
					mean_proba+=float(dict_assignment[key_tocsy][key_exp])
					time_proba+=1
				except KeyError:
					pass
			print 'Mean matching probability: %.4f'%(mean_proba/time_proba)
			print '________'
	return None

def best_assignment(dict_patterns, dict_assignment):
	dict_best_patterns = {}
	nb_patterns = 0
	for key_pattern in dict_patterns.keys():
		best_pattern = []
		pattern = dict_patterns[key_pattern]
		list_tocsy_in_pattern = set(item[0] for item in pattern) # the set of TOCSY peaks found in the pattern
		for peak_tocsy in list_tocsy_in_pattern:
			max_proba_matching = 0.0
			best_match = None
			for item in pattern:
				if item[0] == peak_tocsy:
					proba_matching = dict_assignment[peak_tocsy][item[1]]
					if proba_matching > max_proba_matching:
						max_proba_matching = proba_matching
						best_match = item[1]
					else: 
						pass
			best_pattern.append((peak_tocsy, best_match))
		nb_patterns +=1
		dict_best_patterns[nb_patterns]=best_pattern
	return dict_best_patterns

			
if __name__ == '__main__':
	peak_list = {}
	path_local = str(os.getcwd())
	#peak_list, volume_list = spectrum(path_local + "/43-sn10.peaks")
	peak_list, volume_list = spectrum(path_local + '/threonine.peaks')
	#list_file = database_bmrb()
	#list_file = [path_local+'/database/bmse000408.str']
	list_file = [path_local + '/database/bmse000049.str']
	for name_file in list_file :
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
		assignment_showing(dict_best_patterns, dict_assignment, peak_list, c, len_tocsy_initial, volume_list)
		#c = update_shift(c, peak_list, dict_assignment, dict_patterns)
		#for key in c.shift.keys():
		#	print key, c.shift[key]
