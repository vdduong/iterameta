### ASSIGNMENT CLUSTERING
### assignment clustering using quality threshold clustering
### Viet Dung Duong, Ruben Staub, Torsten Herrmann


from spectrum import *
from class_metabolite import *
from iterative_matching import *
from assignment_assessment import *
from pylab import *
import matplotlib.pyplot as plt
import math
import re
import os

tol_H = 0.09
peak_list = {}

def list_point_append(dict_assignment, peak_list, metabolite) :
	"""
	create the list of secondary chemical shifts
	scattering (X,Y) in order to take the k-mean algo in count"""
	list_point = []
	for key_tocsy in dict_assignment.keys():
		peak_tocsy = metabolite.tocsy_pattern[key_tocsy] # name the peak in tocsy
		cs_tocsy_1 = float(peak_tocsy[0]) # assigned chemical shifts
		cs_tocsy_2 = float(peak_tocsy[1])
		for key_exp in dict_assignment[key_tocsy].keys():
			cs_exp_1 = float(peak_list[key_exp][0])
			cs_exp_2 = float(peak_list[key_exp][1])
			coord_1 = float('%.4f'%(cs_tocsy_1 - cs_exp_1))
			coord_2 = float('%.4f'%(cs_tocsy_2 - cs_exp_2))
			list_point.append((coord_1, coord_2, key_tocsy, key_exp))
	return list_point

###############################
## QUALITY THRESHOLD CLUSTERING


def distanceBetweenPoints_qt(peak_1, peak_2):
	"""
	define the distance between two points in quality threshold clustering
	"""
	distance = math.sqrt( (peak_1[0] - peak_2[0])**2 + \
							(peak_1[1]-peak_2[1])**2)
	return distance

def distanceBetweenPoint_List_qt(peak_1, list_1) :
	list_distance = []
	for item in list_1 :
		distance = distanceBetweenPoints_qt(peak_1, item)
		list_distance.append(distance)
	list_distance.sort(reverse=True)
	return list_distance[0]

def cluster_generation(list_point, threshold):
	"""
	using quality threshold clustering algorithm to cluster good assignments
	"""
	
	class Cluster:
		def __init__(self, diam=0.0, number=0):
			self.diam = diam
			self.number = number
			self.elements = []
			
	dict_distance = {}
	for point_1 in list_point:
		dict_distance[point_1] = {}
		for point_2 in list_point:
			dict_distance[point_1][point_2] = distanceBetweenPoints_qt(point_1, point_2)
	nb_ = 0
	clusters = {}
	while list_point:
		dict_clusters = {}
		for point in list_point:
			dict_clusters[point]=Cluster()
			for point_ in list_point:
				if dict_distance[point][point_] < threshold: 
					dict_clusters[point].number+=1
					dict_clusters[point].elements.append(point_)
					if dict_clusters[point].diam < dict_distance[point][point_]:
						dict_clusters[point].diam = dict_distance[point][point_]
		diameter_min = float('inf')
		nb_max = 0
		key_out = None
		for key_cluster in dict_clusters.keys():
			if dict_clusters[key_cluster].number > nb_max:
				nb_max = dict_clusters[key_cluster].number
				dict_clusters[key_cluster].diam = diameter_min
				key_out = key_cluster
			elif dict_clusters[key_cluster].number < nb_max:
				pass
			else:
				if dict_clusters[key_cluster].diam < diameter_min:
					diameter_min = dict_clusters[key_cluster]
					key_out = key_cluster
				else: pass
		nb_ +=1
		clusters[nb_] = []
		for item in dict_clusters[key_out].elements:
			clusters[nb_].append(item)
			list_point.remove(item)
	return clusters




def pattern_production(dict_assignment, peak_list, metabolite, len_tocsy_initial, threshold):
	"""
	from the results of qt_clustering, we retain only the patterns that contain over 
	50 percent of the total peak numbers
	"""
	dict_patterns = {}
	list_point = list_point_append(dict_assignment, peak_list, metabolite)
	clusters = cluster_generation(list_point, threshold)
	nb_patterns = 0
	for nb_cluster in clusters.keys():
		cluster = clusters[nb_cluster]
		nb_patterns+=1
		dict_patterns[nb_patterns]= [(item[2], item[3]) for item in cluster]
	dict_patterns = peaks_clustered(len_tocsy_initial, dict_patterns) # assignment cluster assessment
	return dict_patterns

def show_cluster(list_point, cluster):
	'''show the cluster'''
	fig=plt.figure()
	ax_view=fig.add_subplot(1,1,1)
	ax_view.scatter([float(point[0]) for point in list_point], \
		[float(point[1]) for point in list_point], color='green',s=60.)
	ax_view.scatter([float(point[0]) for point in cluster], \
		[float(point[1]) for point in cluster], color='orange',s=60.)
	ax_view.set_xlabel('F1 (ppm)')
	ax_view.set_ylabel('F2 (ppm)')
	ax_view.set_title('Hippuric assignment clustering')
	
	#ax_view.annotate('1#', xy = (0.06,0.02), xytext = (0.07,0.03), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('2#', xy = (0.06,-0.08), xytext = (0.07, -0.07), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('3#', xy = (0.06,-0.082), xytext = (0.07, -0.092), arrowprops=dict(facecolor='black', width = 0.01))
	#ax_view.annotate('1, 2, 3', xy = (-0.08,-0.081), xytext = (-0.06, -0.06), arrowprops=dict(facecolor='black', width = 0.01))
	ax_view.annotate('1', xy = (-0.036,-0.03), xytext = (0.0,-0.05), arrowprops=dict(facecolor='black', width = 0.01))
	ax_view.annotate('2', xy = (-0.04, 0.085), xytext = (-0.1,0.1), arrowprops=dict(facecolor='black', width = 0.01))
	fig.show()
	return None
	
#########################
## TEST of the module
if __name__ == '__main__':
	path_local = str(os.getcwd())
	peak_list, volume_list = spectrum(path_local + "/43-sn10.peaks")
	#peak_list, volume_list = spectrum(path_local+'/threonine.peaks')
	#list_file = database_bmrb()
	name_file = path_local + '/database/bmse000408.str'
	#name_file = path_local + '/database/bmse000049.str'
	list_file = [name_file]
	nb_compounds = 0
	for name_file in list_file :
		metabolite = Metabolite(name_file)
		shift = metabolite.initial_shift(name_file)
		metabolite = Metabolite(name_file, shift)
		metabolite.tocsy_pattern = metabolite.tocsy()
		tocsy_initial = metabolite.tocsy_pattern
		len_tocsy_initial = len(tocsy_initial)
		threshold = 0.05
		dict_assignment = assignment(metabolite, peak_list)
		dict_patterns = pattern_production(dict_assignment, peak_list, metabolite, len_tocsy_initial, threshold)
		if len(dict_patterns) > 0 and len(metabolite.tocsy_pattern) >=4.0 :
			nb_compounds +=1
			list_point = list_point_append(dict_assignment, peak_list, metabolite)
			for key in dict_patterns.keys():
				print key, dict_patterns[key]
		list_point = list_point_append(dict_assignment, peak_list, metabolite) 
		print list_point
		list_point_1 = list(list_point)
		clusters = cluster_generation(list_point_1, threshold)
		cluster = clusters[1]
		show_cluster(list_point, cluster)
		show()
	print nb_compounds


