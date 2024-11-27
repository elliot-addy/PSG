'''
Implementation of 'fast prediciton' algorithm, given by Algorithm 1 in Plumlee
2014, 'Fast prediction of deterministic functions using sparse grid
experimental designs,' generalised to Penalised Sparse Grids.
'''

## SCRIPTS ##

from psg.py import construct_psg
from kernels.py import maternKernel


class MaternPSGEmulator:
	'''
	Object for emulating functions parameterised over [-0.5,0.5]^dim at points
	arranged in a penalised sparse grid using a separable matern kernel.
	'''

	def __init__(
		self,
		dim,
		nu_list,
		lengthscale_list,
		level,
		penalty,
		sigma_list
	):
		# Store variables in object.
		self.dim = dim
		self.nu_list = nu_list
		self.lengthscale_list = lengthscale_list
		self.level = level
		self.penalty = penalty
		self.sigma_list = sigma_list
		
		# Initialise kernels.
		self.construct_kernels

	def construct_kernels(
		self
	):
		'''
		Create and store list of kernel objects for each dimension.
		'''
		kernelObject_list = []
		for i in range(dim):
			kernelObject_list.append(maternKernel(
				nu_list[i],
				lengthscale_list[i],
				sigma_list[i]
			))
		self.kernelObject_list = kernelObject_list
		
	def construct_psg_dict(level, penalty, dim, func):
		'''
		Initialise a dictionary with sparse grid poits as keys (Lists) and 
		2-Lists as values.
		Args:
			level: Int, level of construction of the sparse grid.
			penalty: List of Ints, penalty in each dimension.
			dim: Int, dimension.
			func: Function, objective function for sampling.
		Returns:
			Dictionary with sparse grid points as keys and 2-Lists as values,
			where the first entry of list is func evaluated at the grid point.
		'''
		keys = construct_psg(level, penalty, dim)
		psg_dict = {key: [func(keys), 0] for key in keys}
	
		return psg_dict



if __name__ == '__main__':
	
	def main():
		'''
		Create 
		'''
