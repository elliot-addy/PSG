'''
Script containing all functions required to construct Penalised Sparse Grids 
and correpsonding posterior Gaussian processes.
'''

# Packages.
import numpy as np
import itertools as it

def index_to_design(index):                                                     
	'''
	For a given multi-index, return additional components of sparse grid
	compared to index -1 in all dimensions [REWORK THIS].

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
	for j in range(0,dim):
		# Number of points in each dimension governed by muli-index.
		one_d_arrays.append(
			np.arange(-0.5+0.5**(index[j]+1),0.5,0.5**(index[j]))
 		)
	# Sparse grid component given by cartesian product of one-dim sequences.
	points = it.product(*one_d_arrays)
	return list(points


if __name__ == '__main__':

	def main():
		'''
		Construct Penalised Sparse Grid and plot for dimensions <= 3.
		'''
		print('hello')
		pass

	main()


