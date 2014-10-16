# running file
# Viet Dung Duong, Ruben Staub, Torsten Herrmann

from spectrum import *
from iterative_matching import *
from assignment_assessment import *
from assignment_clustering import *
from assignment_plotting import *
import time

start = time.clock()

tol_H = 0.1
threshold = 0.05

number_compounds_found = 0
print '\n'*40
peak_list = {}
path_local = str(os.getcwd())
peak_list, volume_list = spectrum(path_local + "/43-sn10.peaks")
#peak_list, volume_list = spectrum(path_local + '/PeakList.peaks')
#peak_list, volume_list = spectrum(path_local + '/hippuric.peaks')
#peak_list, volume_list = spectrum(path_local+'/metahunter.peaks')
#peak_list, volume_list = spectrum(path_local + '/synthetic_20.peaks')
list_file = database_bmrb()
list_name_1 = set() # list of metabolite names 
for name_file in list_file :
	io_compound = True
	dict_assignment = {}
	metabolite = Metabolite(name_file)
	shift = metabolite.initial_shift(name_file)
	metabolite = Metabolite(name_file, shift)
	metabolite.tocsy_pattern = metabolite.tocsy()
	tocsy_initial = metabolite.tocsy_pattern
	len_tocsy_initial = len(tocsy_initial)
	dict_assignment = assignment(metabolite, peak_list)
	if peaks_found(dict_assignment) >= 0.6*len(tocsy_initial) and len(tocsy_initial) >=4:
		nb_cycle = 0
		while nb_cycle <= 4 and io_compound==True:
			dict_assignment = assignment(metabolite, peak_list)
			dict_patterns = pattern_production(dict_assignment, peak_list, metabolite, len_tocsy_initial, threshold)
			dict_patterns = peaks_clustered(len(tocsy_initial), dict_patterns)
			if len(dict_patterns) == 0: 
				io_compound = False
				break
			metabolite = update_shift(metabolite, peak_list, dict_assignment, dict_patterns)
			tocsy_pattern = metabolite.tocsy()
			if len(tocsy_pattern) < len_tocsy_initial: 
				io_compound = False
				break
			nb_cycle+=1

		if io_compound:
			fractionHausdorff = 0.8
			hausdorff = hausdorffDistance(fractionHausdorff, metabolite, dict_assignment, peak_list)
			if hausdorff < 0.05 and io_compound:
				dict_best_patterns = best_assignment(dict_patterns, dict_assignment)
				#assignment_showing(dict_best_patterns, dict_assignment, peak_list, metabolite, len_tocsy_initial, volume_list)
				#print 'score of %.2f -Hausdorff distance: %.4f'%(fractionHausdorff*100, hausdorff)
				#print '_'*20
				name = metabolite.name
				if name not in list_name_1:
					number_compounds_found+=1
					list_name_1.add(name)
			else:pass
	else: pass
	
print 'In total, %i compounds are found'%(number_compounds_found)
end = time.clock()
print '%s'%(end-start)
print list_name_1
# plotting

# in the result, the initial shifts, the shifted shifts and the corresponding standard deviation 

peak_list, volume_list = spectrum(path_local + '/PeakList.peaks')
#peak_list, volume_list = spectrum(path_local + '/hippuric.peaks')
#peak_list, volume_list = spectrum(path_local+'/metahunter.peaks')
#peak_list, volume_list = spectrum(path_local + '/synthetic_20.peaks')
#list_file = database_bmrb()
list_name_2 = set() # list of metabolite names 
for name_file in list_file :
	io_compound = True
	dict_assignment = {}
	metabolite = Metabolite(name_file)
	shift = metabolite.initial_shift(name_file)
	metabolite = Metabolite(name_file, shift)
	metabolite.tocsy_pattern = metabolite.tocsy()
	tocsy_initial = metabolite.tocsy_pattern
	len_tocsy_initial = len(tocsy_initial)
	dict_assignment = assignment(metabolite, peak_list)
	if peaks_found(dict_assignment) >= 0.6*len(tocsy_initial) and len(tocsy_initial) >=4:
		nb_cycle = 0
		while nb_cycle <= 4 and io_compound==True:
			dict_assignment = assignment(metabolite, peak_list)
			dict_patterns = pattern_production(dict_assignment, peak_list, metabolite, len_tocsy_initial, threshold)
			dict_patterns = peaks_clustered(len(tocsy_initial), dict_patterns)
			if len(dict_patterns) == 0: 
				io_compound = False
				break
			metabolite = update_shift(metabolite, peak_list, dict_assignment, dict_patterns)
			tocsy_pattern = metabolite.tocsy()
			if len(tocsy_pattern) < len_tocsy_initial: 
				io_compound = False
				break
			nb_cycle+=1

		if io_compound:
			fractionHausdorff = 0.8
			hausdorff = hausdorffDistance(fractionHausdorff, metabolite, dict_assignment, peak_list)
			if hausdorff < 0.05 and io_compound:
				dict_best_patterns = best_assignment(dict_patterns, dict_assignment)
				#assignment_showing(dict_best_patterns, dict_assignment, peak_list, metabolite, len_tocsy_initial, volume_list)
				#print 'score of %.2f -Hausdorff distance: %.4f'%(fractionHausdorff*100, hausdorff)
				#print '_'*20
				name = metabolite.name
				if name not in list_name_2:
					number_compounds_found+=1
					list_name_2.add(name)
			else:pass
	else: pass
	
print 'In total, %i compounds are found'%(number_compounds_found)
end = time.clock()
print '%s'%(end-start)

commun = list_name_1.intersection(list_name_2)
print len(commun)
print commun
