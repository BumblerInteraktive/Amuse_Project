from amuse.lab import read_set_from_file
from amuse.units import units
from scipy.spatial import cKDTree
from scipy import stats
import numpy
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
import os

def read_from_file(data_path, index):
    filename = os.path.join(data_path,"gas_particle_data_i{0:04}.hdf5".format(index.value_in(units.Myr)))
    filename2 = os.path.join(data_path,"star_particle_data_i{0:04}.hdf5".format(index.value_in(units.Myr)))
    gas_particles = read_set_from_file(filename, format="amuse")
    star_particles = read_set_from_file(filename2, format="amuse")
    return gas_particles, star_particles
def create_density_histogram(gas,t_end, filename):
    xyz = numpy.vstack([gas.x.value_in(units.kpc),gas.y.value_in(units.kpc),gas.z.value_in(units.kpc)])
    kde = stats.gaussian_kde(xyz)
    densities = 1e4*kde(xyz)
    densities*=1.989e33
    densities/=(3.086e21)**3
    # Create and save the histogram
    plt.hist(densities, bins=30, color='blue', edgecolor='black', range=(0.0,7.5e-29))
    plt.xlabel('Density (g/cm^3)')
    plt.ylabel('Count')
    plt.title('Gas Particle Densities at t='+str(int(t_end.value_in(units.Myr)))+'Myr')
    plt.savefig(filename)
    plt.close()

def density_scatter( stars ,gas,filename,sort = True, bins = 20,**kwargs )   :

    """
    Scatter plot colored by 2d histogram
    """
    gasx,gasy=gas.x.value_in(units.kpc), gas.y.value_in(units.kpc)
    starsx,starsy=stars.x.value_in(units.kpc), stars.y.value_in(units.kpc)
    fig,ax= plt.subplots()
    
    data , x_e, y_e = numpy.histogram2d(gasx, gasy, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , numpy.vstack([gasx,gasy]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[numpy.where(numpy.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = gasx[idx], gasy[idx], z[idx]
    ax.scatter(starsx,starsy, alpha=0.5, s=1, lw=0)
    ax.scatter( x, y, c=z, **kwargs,alpha=0.5, s=1, lw=0 )
    plt.xlim(-50, 50)
    plt.ylim(-50, 50)

    norm = Normalize(vmin = 0, vmax = 1)
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm),ax=ax)
    cbar.ax.set_ylabel('Density')

    fig.savefig(filename)  
    plt.close()


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--num_fig", dest="num_fig", default = 60,
                      help="Number of data snapshots made [%default]")
    result.add_option("--dir", dest="datapath", default = os.path.join("..", "data", "raw", "data_00velocity"),
                      help="Relative directory path to data")
    result.add_option("--t_end", unit=units.Myr,
                      dest="t_end", default = 3000|units.Myr,
                      help="End of the simulation [%default]")
    return result

if __name__ == '__main__':
    o, arguments  = new_option_parser().parse_args()
    t_ends = numpy.linspace(o.t_end.value_in(units.Myr)/o.num_fig, o.t_end.value_in(units.Myr), o.num_fig) |units.Myr
    plot_path=os.path.join(o.datapath, "plots")
    os.makedirs(plot_path)
    for t_end in t_ends:
        gas, stars = read_from_file(o.datapath, t_end)
        create_density_histogram(gas,t_end, filename=os.path.join(plot_path, "histogram_t"+str(int(t_end.value_in(units.Myr)))+"Myr.png"))
        density_scatter(stars, gas,  filename=os.path.join(plot_path, "densities_scatter_t"+str(int(t_end.value_in(units.Myr)))+"Myr.png"))


