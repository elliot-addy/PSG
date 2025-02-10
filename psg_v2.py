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


def construct_sg(level, dim):
	'''
	Constructs the set of points in an isotropic sparse grid with chosen
	parameters.
	'''
	# [TO DO] FIRST CHECK IF SPARSE GRID IS ALREADY STORED!!!!!!!!!!!!!!!!!!!!!
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
	# Initialise design with single point at origin.
	design = [np.zeros(dim)]

	## ANISOTROPIC PART ##

	# Iterate over all non-trivial dimensions.
	for k in range(1, dim+1):
		# Check if any compenents of dimension k exist.
		if level >= k:
			# Iterate over all subsapces of dimension k.
			for J in it.combinations(range(1, dim+1), k):
				print('J')
				print(J)
				# Derive penalty of component.
				component_penalty = 0
				for j in J:
					component_penalty += penalty[j-1]
				print('component penalty')
				print(component_penalty)
				# Check if component is non-trivial in this direction
				if level >= k + component_penalty:
					# Create component isotropic sparse grid of required size.
					sg_component = construct_sg(level-component_penalty-k, k)
					print('sg component')
					print(sg_component)
					# Iterate over each sub-region in k-dim subspace J, either
					# side of each axis.
					addition = []
					for region in it.product([-1,1], repeat=k): #?
						print('region')
						print(region)
						# Iterate over each point.
						for point in sg_component:
							print('point')
							print(point)
							# Iterate over each coordiante.
							projected_point = np.zeros(dim)
							for i,j in enumerate(J):
								projected_point[j-1] += 0.25*region[i] \
									+ 0.5*point[i] 
								print('projected point')
								print(projected_point)
								design.append(projected_point)
	return design


if __name__ == '__main__':

	## PACKAGES ##
	import matplotlib.pyplot as plt

	def main():
		'''
		Construct Penalised Sparse Grid and plot for dimensions <= 3.
		'''
		
		# PSG parameters.
		dim = 2
		level = 4
		penalty = [1,2,3,4]

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
