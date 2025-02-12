'''
Implementation of 'fast prediciton' algorithm, given by Algorithm 1 in Plumlee
2014, 'Fast prediction of deterministic functions using sparse grid
experimental designs,' generalised to Penalised Sparse Grids.
'''

## PACKAGES ##

import numpy as np
from scipy.linalg import cholesky, toeplitz, solve_triangular
from math import comb
from itertools import product
from operator import itemgetter

## SCRIPTS ##

from psg_v2 import construct_psg, to_sum_k_rec
from kernels import maternKernel


class MaternPSGEmulator:
	'''
	Object for emulating functions parameterised over [-0.5,0.5]^dim at points
	arranged in a penalised sparse grid using a separable matern kernel.
	'''

	def __init__(
		self,
		dim,
		level,
		penalty,
		func,
		nu_list,
		lengthscale_list,
		sigma_list
	):
		# Store variables in object.
		self.dim = dim
		self.func = func
		self.nu_list = nu_list
		self.lengthscale_list = lengthscale_list
		self.level = level
		self.penalty = penalty
		self.sigma_list = sigma_list
		
		# Initialise kernels.
		self.construct_kernels()

		# Construct penalised sparse grid dictionary.
		self.construct_psg_dict()


	def construct_kernels(self):
		'''
		Create and store list of kernel objects for each dimension.
		'''
		kernelObject_list = []
		for i in range(self.dim):
			kernelObject_list.append(maternKernel(
				self.nu_list[i],
				self.lengthscale_list[i],
				self.sigma_list[i]
			))
		self.kernel_list = kernelObject_list
		# Reset matrix inverse dictionary. [MAKE THIS NICER]
		self.L_dict = {}
		
	def construct_psg_dict(self):
		'''
		Initialise a dictionary with sparse grid poits as keys (Lists) and 
		2-Lists as values.
		Args:
		self.level: Int, level of construction of the sparse grid.
			self.penalty: List of Ints, penalty in each dimension.
			self.dim: Int, dimension.
			self.func: Function, objective function for sampling.
		Returns:
			Dictionary with sparse grid points as keys and 2-Lists as values,
			where the first entry of list is func evaluated at the grid point.
		'''
		keys = construct_psg(self.level, self.penalty, self.dim)
		self.dict = {tuple(key): [self.func(key),0] for key in keys}
		# Reset matrix inverse dictionary. [MAKE THIS NICER]
		self.L_dict = {}

	def contribution_scaling(self, multi_index):
		'''
		Scales the contribution of each sparse grid component for inversion
		calculation. Function 'a' in PLumlee 2014.

		Args:
			self: MaternPSGEmulator object.
			multi_index: List of self.dim many non-negative integers.
		Returns:
			Float, scales contribution from given multi-index's sub-system.
		'''
		# Find 1-norm of multi-index.
		norm = np.sum(multi_index)
		return (-1)**(self.level-norm)*comb(self.dim-1,self.level-norm)
		
	def ordered_index_to_design(self, multi_index):
		'''
		Order points in cartesian products to match order of Kronecker system.

		Args:
			self: MaternPSGEmulator object.
			multi_index: List of self.dim many non-negative integers.
		Returns:
			List of grid-points , ordered with respect to Kronecker linear
			system.
		'''
		# Initialise list of 1D designs.
		design_list = []
		# Iterate through dimensions.
		for j in range(self.dim):

			exponent = max(1, multi_index[j] + 1 - self.penalty[j])
			design_list.append(np.arange(
				-0.5 + 0.5**exponent,
				0.5,
				0.5**exponent
			))

		return list(product(*design_list))
			

	def calculate_weights(self):
		'''
		Calculate weight vector. Adaptation of 'Algorithm 1' in Plumlee 2014. 
		'''
		# Create list of multi-indices required, P(L):
		multi_index_list = []
		for l in range(self.level+1):
			multi_index_list += list(to_sum_k_rec(self.dim, l))

		# Iterate over all multi-indices.
		for multi_index in multi_index_list:
			# Find subset of sparse grid points corresponding to multi-index,
			# arranged in Kronecker order.
			print('multi-index')
			print(multi_index)

			keys = self.ordered_index_to_design(multi_index)
			# Initilaise inverse Kronecker product matrix.
			L = np.ones(1)
			# Iterate over each index corresponding to a dimension.
			for j, index in enumerate(multi_index):
				# Find penalty in dimension j.
				penalty = self.penalty[j]
				# Check if 1D covar matrix Cholesky decomposition is stored in 
				# matrix dictionary.
				if (j, index-penalty) in self.L_dict:
					L_component = self.L_dict[(j, index-penalty)]
				else:
					# If inverse 1D covar matrix is not stored, calculate using
					# Cholesky decomposition.
					# Determine size of matrix.
					n = int(2**(index+1-penalty) - 1)
					# Determine 1D design.
					one_d_array = np.arange(
						-0.5+0.5**(index+1-penalty),
						0.5,
						0.5**(index+1-penalty)
					)
					# Initalise first row of Toeplitz matrix.
					column = [1]
					# Iterate through column.
					for i in range(1,n):
						column.append(self.kernel_list[j](
							one_d_array[i],
							one_d_array[0]
						))
					# Construct matrix.
					covar_mat = toeplitz(column,column)
					# Cholesky decomposition.
					L_component = cholesky(covar_mat, lower = True)
					print('L_component')
					print(L_component)
					# Store matrix in dictionary.
					self.L_dict[(j, index-penalty)] = L_component

				# Successively build Kronecker product of Cholesky matrices.
				L = np.kron(L, L_component)
			
			# Solve triangular systems.
			# [THINK ABOUT KEYS][SPARSE PRODUCT?]
			print('keys')
			print(keys)
			print(itemgetter(*keys)(self.dict))
			print(L)
			data, old_weights = np.array(itemgetter(*keys)(self.dict)).T
			print(data)
			if isinstance(data, float):
				data, old_weights = [data], [old_weights]
			b = solve_triangular(
				L,
				np.array(data),
				lower = True
			)
			weight_update = solve_triangular(
				L.T,
				b
			)
			# Calculate weight update.
			new_weights = old_weights + \
				self.contribution_scaling(multi_index) * weight_update
			# Store new weights in dictionary.
			self.dict.update(list(zip(keys, list(zip(data, new_weights)))))

if __name__ == '__main__':
	
	## PACKAGES ##
	import matplotlib.pyplot as plt

	def main():
		'''
		Create MaternPSGEmulator object and test.
		'''
		# Sparse grid parameters.
		dim = 2
		level = 0
		penalty = [0,0,0,0]
		
		# Emulator parameters.
		nu_list = [0.5,1.5,0.5]
		lengthscale_list = list(2**np.array(penalty[:dim]))
		sigma_list = [1,1,1]

		# Define test funciton.
		def func(x):
			'''
			Test function used in place of simulator.

			Args:
				x: List of Floats, point in [-0.5,0.5]^dim.
			Returns:
				Float.	
			'''
			return 1.0 #x[0]**2 + x[1]**2 
		
		# Construct MaternPSGEmulator object.
		PSGEmulatorObject = MaternPSGEmulator(
			dim,
			level,
			penalty,
			func,
			nu_list,
			lengthscale_list,
			sigma_list
		)

		print(PSGEmulatorObject.dict)
		PSGEmulatorObject.calculate_weights()
		print(PSGEmulatorObject.dict)

		if dim == 1:
			fig, ax = plt.subplots()
			points = list(PSGEmulatorObject.dict.keys())
			print('points')
			print(points)
			weights = [PSGEmulatorObject.dict[point][1] for point in points]
		
			def	kernel_approximation(MaternObject, arg, points, weights):
				'''
				Kernel approximation
				'''
				total = 0
				for i in range(len(points)):
					total +=  weights[i]*MaternObject.kernel_list[0](
									points[i][0],
									arg
									)
				return total
		
			x_list = np.linspace(-0.5,0.5,100)
			y_list = [kernel_approximation(
					PSGEmulatorObject,
					x,
					points, 
					weights
				) for x in x_list]

			ax.plot(x_list,y_list)
			plt.show()
			print(y_list)

		elif dim == 2:
			fig = plt.figure()
			ax = fig.add_subplot(projection='3d')
			points = list(PSGEmulatorObject.dict.keys())
			print('points')
			print(points)
			weights = [PSGEmulatorObject.dict[point][1] for point in points]
			print('weights')
			print(weights)
			def kernel_approximation(
					MaternObject,
					arg1,
					arg2,
					points,
					weights
				):
				'''
				Kernel approximation.
				'''
				total = 0
				for i in range(len(points)):
					total += weights[i] *\
									MaternObject.kernel_list[0](
										points[i][0],
										arg1
									) *\
									MaternObject.kernel_list[1](
										points[i][1],
										arg2
									)
				return total

			x = y = np.linspace(-0.5,0.5,100)
			X, Y = np.meshgrid(x, y)
			z = np.array(kernel_approximation(
							PSGEmulatorObject,
							np.ravel(X),
							np.ravel(Y),
							points,
							weights
						)
					)

			Z = z.reshape(X.shape) 
			
			ax.plot_surface(X, Y, Z)
			plt.show()

	main()
	


