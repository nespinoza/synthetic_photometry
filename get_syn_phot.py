import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import UnivariateSpline


####################### OPTIONS ########################

input_metallicities = ['m01','m02','m03','m05','m10',\
                       'm15','m20','m25','m30','m35',\
                       'm40','m45','m50','p00','p01',\
                       'p02','p03','p05','p10']

input_vturb = ['0.0']
band = 'i_sdss'

########################################################
fluxes_folder = 'model_fluxes'

# First, take care of the bandpass; generate a univariate 
# spline that spans it:
w,R = np.loadtxt('bandpasses/'+band+'.dat')
S = UnivariateSpline(w,R,s=0,k=1)

# Integrate the response function to get the normalization:
Sintegrand = lambda x,S: S(x)
normalization = (integrate.quad(Sintegrand,np.min(w),np.max(w), \
                 args=(S,),full_output=1))[0]

# Now iterate through each spectrum, get magnitudes, save:
for met in input_metallicities:
    
