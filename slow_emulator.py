## PACKAGES ##
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_triangular

## MODULES ##
from psg_v2 import construct_psg
from kernels import maternKernel


def main():
	'''
	Constructs a basic GP emulator of a funciton f in (-0.5,0.5)^d with zero 
	prior mean and calculates L2 
	error.
	'''
	## PARAMETERS ##

	# Intepolation parameters.
	dim = 2
	level = 4
	penalty = np.arange(0, dim, 1, dtype=np.float64)

	# Kernel parameters.
	nu_array = 0.5*np.ones(dim)
	lengthscale_array = 2**penalty
	sigma_array = np.ones(dim)

	def func(x):
		'''
		Test function.
		'''
		return 1.

	# Plotting parameters.
	res = 100

	## DATA ##

	x_data = np.array(construct_psg(level, penalty, dim))
	f_data = np.array([func(x) for x in x_data])

	em = emulator(
			dim,
			x_data,
			f_data,
			nu_array,
			lengthscale_array,
			sigma_array
		)

	print('Plotting...')

	if dim == 1:
		x_array = np.linspace(-0.5, 0.5, res)
		m_array = np.array([em([x]) for x in x_array])
		plt.plot(x_array, m_array)
		plt.show()

	if dim == 2:
		x_array = y_array = np.linspace(-0.5, 0.5, res)
		X, Y = np.meshgrid(x_array, y_array)
		z_array = np.array(
					[em(x) for x in np.array([np.ravel(X), np.ravel(Y)]).T]
				)
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')
		Z = z_array.reshape(X.shape)
		ax.plot_surface(X, Y, Z)
		plt.show()

class emulator:
	
	def __init__(
		self,
		dim,
		x_data,
		f_data,
		nu_array,
		lengthscale_array,
		sigma_array
	):
		'''
		Attributes
		----------

		dim : int
			Dimension of function domain.
		x_data : array-like
			Points in (-0.5,0.5)^d to be interpolated at.
		f_data : array-like
			Function evaluations at x_data.
		nu_array : array-like
			Smoothness parameter of kernel in each axial direction.
		lengthscale_array : array-like
			Lengthscale parameter of kernel in each axial direction.
		sigma_array : array-like
			Standard deviation parameter of kernel in each axial direction.
		kernel_array : array-like of maternKernel objects.
			Represents separable Matern kernel.
		cholesky : array-like
			Lower triangular Cholesky decomposition of covariance matrix.
		weights : array-like
			Weights for kernel interpolant.
		'''

		self.dim = dim
		self.x_data = x_data
		self.f_data = f_data
		self.nu_array = nu_array
		self.lengthscale_array = lengthscale_array
		self.sigma_array = sigma_array
		
		# Constuct and store kernel_array.
		self.construct_kernel()
	
		# Construct and store Cholesky factorisation of covariance matrix.
		self.construct_cholesky()

		# Calculate weights and store.
		self.calculate_weights()
		
	def construct_kernel(self):
		'''
		Construct np.array representing separabel Matern kernel.
		'''
		kernel_array = np.zeros(self.dim, dtype=object)
		# Loop through each dimension.
		for j in range(self.dim):
			kernel_array[j] = maternKernel(
									self.nu_array[j],
									self.lengthscale_array[j],
									self.sigma_array[j]
								)
		self.kernel_array = kernel_array

	def covariance(self, x1, x2):
		'''
		Calcualtes the covariance between x1 and x2 in (-0.5,0.5)^d.
		'''
		log_covar = 0
		for i in range(self.dim):
			log_covar += np.log(self.kernel_array[i](x1[i], x2[i]))
		return np.prod(self.sigma_array) * np.exp(log_covar)
	
	def construct_cholesky(self):
		'''
		Construct covariance kernel and store Cholesky factorisation.
		'''
		print('Building Cholesky matrix... ')
		# Intitialise matrix.
		covar_matrix = np.zeros(
							(self.x_data.shape[0], self.x_data.shape[0]),
							dtype=np.float64
						)
		for i in range(self.x_data.shape[0]):
			for k in range(self.x_data.shape[0]):
				covar_matrix[i,k] = self.covariance(
										self.x_data[i],
										self.x_data[k]
									)
		self.cholesky = np.linalg.cholesky(covar_matrix)
		print('Done.')
	
	def calculate_weights(self):
		'''
		Calcualtes weights given Cholesky factrisation of covariance matrix.
		'''
		print('Solving Cholesky system...')
		y = solve_triangular(self.cholesky, self.f_data, lower=True)
		weights = solve_triangular(self.cholesky.T, y, lower=False)
		
		self.weights = weights
		print('Done.')

	def __call__(self, arg):
		'''
		Returns value of interpolant at given point in (-0.5,0.5)^d.
		'''
		m = 0
		for i, x in enumerate(self.x_data):
			m += self.weights[i] * self.covariance(arg, x)
		return m
		
		
if __name__ == '__main__':
	main()
		

