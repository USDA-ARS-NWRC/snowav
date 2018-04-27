import netCDF4 as nc
import numpy as np


def output_nc(snow):
    '''
    Outputs attributes of snow identified in config file as netcdf file
    '''   
    for var in snow.ncvars:
        
        print('Saving %s to netcdf in %s%s%s'%(var,snow.figs_path,var,snow.name_append))
        p = getattr(snow,var)
        ny = len(p[:,0])
        nx = len(p[0,:])
        my_data = np.reshape(p,(nx,ny),'F')
        f = nc.Dataset('%s%s%s.nc'%(snow.figs_path,var,snow.name_append), 'w')
        f.createDimension('nx',nx)
        f.createDimension('ny',ny)
        data = f.createVariable(var, 'f4', ('nx','ny'))
        data[:] = np.flipud(my_data)
        f.close()  
    
        