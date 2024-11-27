'''
Script containing all functions required to construct Penalised Sparse Grids
'''

## PACKAGES ##

import numpy as np
import itertools as it

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


def construct_sg_differences(level, dim):
	'''
	Constructs the set of points in SG(level) / SG(level-1), i.e. the
	additional points added to successive isotropic sparse grids.
	'''
	# Return empty list if level is negative.
	if level < 0:
		return [] # [CHANGED FROM ORIGINAL]
	else:
		# Initialise empty list of points.
		point_list = []
		# Iterate over all additional indices.
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
	# Initialise design as empty list.
	design = []

	## ANISOTROPIC PART ##

	# Iterate over all dimensions.
	for j in range(1, dim+1):
		# Iterate over differences between successive penalties.
		for l in range(1, penalty[j]-penalty[j-1]+1):
			# Find addtional points between levels, done in a lower dimension.
			addition = construct_sg_differences(level-penalty[j]+l,j)
			# Embed new points in higher ambient dimension. [IMPROVE]
			ammended_addition = []
			for point in addition:
				ammended_addition.append(
					list(point)+list(np.zeros(dim-j))
				)
			# Append new points to design.
			design += ammended_addition
	
	## ISOTROPIC PART ##

	# Check if isotropic part is non-empty.
	if level > penalty[dim-1]:
		# Iterate over levels up to level of construciton of isotropic part,
		# given by difference between overall level and size of largest
		# penalty.
		for l in range(level-penalty[dim-1]):
			# Append points from each level of isotropic part to design.
			design += construct_sg_differences(l, dim)

	return design


if __name__ == '__main__':

	## PACKAGES ##
	import matplotlib.pyplot as plt

	def main():
		'''
		Construct Penalised Sparse Grid and plot for dimensions <= 3.
		'''
		
		# PSG parameters.
		dim = 3
		level = 5
		penalty = [0,1,2,3]

		# Plotting parameters.
		print_design = True

		# Function.
		design = np.array(construct_psg(level, penalty, dim))
		if print_design == True:
			print(design)

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
