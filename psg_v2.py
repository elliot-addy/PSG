'''
Script containing all functions required to construct Penalised Sparse Grids
'''

## PACKAGES ##

import numpy as np
import itertools as it
import pickle

## DEFINE FUNCTIONS ##

def to_sum_k_rec(dim, level):
	'''
	Returns list of all combinations of dim-many non-negative integers that sum
	to level. [ADAPTED FROM ONLINE FORUM - REFERENCE]
	'''
	if dim == 1:
		yield (level,)
	else:
		for x in range(0, level+1):
			for i in to_sum_k_rec(dim - 1, level - x):
				yield (x,) + i


def index_to_design(index):                                                     
	'''
	For a given multi-index, return additional components of sparse grid not
	included by indices closed-below index.

	Args:
		index: List of dim-length Lists representing a multi-index for a sparse
		grid component.
	Returns:
		List of dim-length Lists, each representing a point in dim-dimensional
		space.
	'''
	dim  = len(index)
	# Initialise list of one-dimensional sequences.
	one_d_arrays = []
	# Iterate over each dimension.
	for j in range(0,dim):
		# Number of points in each dimension governed by muli-index.
		one_d_arrays.append(
			np.arange(-0.5+0.5**(index[j]+1),0.5,0.5**(index[j])))
	# Sparse grid component given by cartesian product of one-dim sequences.
	points = it.product(*one_d_arrays)

	return list(points)


def construct_sg(level, dim):
	'''
	Constructs the set of points in an isotropic sparse grid with chosen
	parameters.
	'''
	# Initialise empty list of points.
	point_list = []
	# Iterate over all levels.
	for level in range(0,level+1):
		# Iterate over all indices.
		for index in list(to_sum_k_rec(dim, level)):
			# Append points from each index to point list.
			point_list += index_to_design(list(index))

	return point_list


def construct_psg(level, penalty, dim):
	'''
	Construct penalised sparse grid, points stored as list.

	Args:
		level: Int, level of construciton of the sparse grid.
		penalty: List of Ints, penalty in each dimension.
		dim: Int, dimension.
	Returns:
		List of Tuples, coordiates of points in penalised sparse grid.
	'''
	# Initialise save bool.
	issaved = False
	ischanged = False

	# Open saved dictionary of isotropic designs.
	try:
		with open('sg_dict.pkl', 'rb') as sg_dict_file:
			sg_dict = pickle.load(sg_dict_file)
		issaved = True
	except:
		print('Unable to access saved sparse grids.')

	# Initialise design with single point at origin.
	design = [np.zeros(dim)]

	## ANISOTROPIC PART ##

	# Iterate over all non-trivial dimensions.
	for k in range(1, dim+1):
		# Check if any compenents of dimension k exist.
		if level >= k:
			# Iterate over all subsapces of dimension k.
			for J in it.combinations(range(1, dim+1), k):
				# Derive penalty of component.
				component_penalty = 0
				for j in J:
					component_penalty += penalty[j-1]
				# Check if component is non-trivial in this direction
				if level >= k + component_penalty:
					# Using saved dictionary?
					if issaved == True:
						# Check if sparse grid component is saved.
						try:
							sg_component = sg_dict[
									(level - component_penalty - k, k)
								]
						# If not, construct and save.
						except:
							ischanged = True
							print(
								f'Saving new component, level \
							{level - component_penalty - k}, dimension = {k}'
						)
							# Create component isotropic sparse grid of 
							# required size.
							sg_component = construct_sg(
									level - component_penalty - k, k
								)
							# Save new component in dictionary.
							sg_dict[(level - component_penalty - k, k)] \
								= sg_component
							print('Saved.')
					else:
						# If not using saved dictionary, construct component.
						sg_component = construct_sg(
								level - component_penalty - k, k
							)

					# Iterate over each sub-region in k-dim subspace J, either
					# side of each axis.
					for region in it.product([-1,1], repeat = k): #?
						# Iterate over each point.
						for point in sg_component:
							# Iterate over each coordiante.
							projected_point = np.zeros(dim)
							for i,j in enumerate(J):
								projected_point[j-1] += 0.25*region[i] \
									+ 0.5*point[i]
							design.append(projected_point)

	# Store sparse grid dictionary if changes have been made.
	if ischanged == True:
		with open('sg_dict.pkl', 'wb') as sg_dict_file:
			pickle.dump(sg_dict, sg_dict_file)
	
	return design


if __name__ == '__main__':

	## PACKAGES ##
	import matplotlib.pyplot as plt

	def main():
		'''
		Construct Penalised Sparse Grid and plot for dimensions <= 3.
		'''
		
		# PSG parameters.
		dim = 100
		level = 5
		penalty = [i-1 for i in range(dim)]

		print_design = True

		# Function.
		design = np.array(construct_psg(level, penalty, dim))
		if print_design == True:
			print(design)
			print(f'Num points = {np.shape(design)}')

		# Plotting 2D.
		if dim == 2:
			fig, ax = plt.subplots()
			ax.grid()
			ax.set_ylim(-0.5,0.5)
			ax.set_xlim(-0.5,0.5)
			ax.scatter(*design.T)
			plt.show()
		# Plotting 3D.
		elif dim == 3:
			fig = plt.figure()
			ax = fig.add_subplot(projection='3d')
			ax.set_ylim(-0.5,0.5)
			ax.set_xlim(-0.5,0.5)
			ax.set_zlim(-0.5,0.5)
			ax.scatter(*design.T)
			plt.show()

	main()
