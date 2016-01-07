import glob
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import UnivariateSpline


####################### OPTIONS ########################

output_fname = 'mags_i_band.dat'
input_metallicities = ['m01','m02','m03','m05','m10',\
                       'm15','m20','m25','m30','m35',\
                       'm40','m45','m50','p00','p01',\
                       'p02','p03','p05','p10']

input_vturb = '2.0'
band = 'i_sdss'

########################################################
fluxes_folder = 'model_fluxes'

# First, take care of the bandpass; generate a univariate 
# spline that spans it:
w,R = np.loadtxt('bandpasses/'+band+'.dat',unpack=True)
S = UnivariateSpline(w,R,s=0,k=1)

# Integrate the response function to get the normalization:
min_w = np.min(w)
max_w = np.max(w)
Sintegrand = lambda x,S: S(x)
normalization = (integrate.quad(Sintegrand,min_w,max_w, \
                 args=(S,),full_output=1))[0]

# Prepare flux integral:
integrand = lambda x,F,S,norm: S(x)*F(x)/norm

# Prepare output file:
fout = open('output_files/'+output_fname,'w')

fout.write('#Temp (K)\tMet\tGrav\tVturb\tMagnitude\n')
# Now iterate through each spectrum, get magnitudes, save:
for met in input_metallicities:
    if met[0] == 'p':
       current_metallicity = met[1]+'.'+met[2]
    else:
       current_metallicity = '-'+met[1]+'.'+met[2]
    print current_metallicity
    all_folders = glob.glob(fluxes_folder+'/'+met+'/*')
    for folder in all_folders:
        current_temperature = folder.split('/')[-1]
        all_files = glob.glob(folder+'/*vturb_'+input_vturb+'_*')
        for current_file in all_files:
            current_gravity = current_file.split('flux_grav_')[1].split('_')[0]
            w,f = np.loadtxt(current_file,unpack=True) 
            Flux = UnivariateSpline(w,f,s=0,k=1)
            integrated_flux = (integrate.quad(integrand,min_w,max_w,\
                               args=(Flux,S,normalization),full_output=1))[0]
            magnitude = -2.5*np.log10(integrated_flux)
            fout.write(current_temperature+'\t'+current_metallicity+'\t'+current_gravity+'\t'+input_vturb+'\t'+str(magnitude)+'\n')
fout.close()
