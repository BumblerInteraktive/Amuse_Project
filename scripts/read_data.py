from amuse.lab import read_set_from_file
from amuse.units import units
from making_channel_hydro_stars import make_plot
from scipy.spatial import cKDTree
from scipy import stats
import numpy
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn

def read_from_file(index):
    filename = "data_60velocity/gas_particle_data_i{0:04}.hdf5".format(index.value_in(units.Myr))
    filename2 = "data_60velocity/star_particle_data_i{0:04}.hdf5".format(index.value_in(units.Myr))
    # filename = "/home/lea/Amuse-env/gas_particle_data_i{0:04}.hdf5".format(index.value_in(units.Myr))
    # filename2 = "/home/lea/Amuse-env/star_particle_data_i{0:04}.hdf5".format(index.value_in(units.Myr))
    gas_particles = read_set_from_file(filename, format="amuse")
    star_particles = read_set_from_file(filename2, format="amuse")
    return gas_particles, star_particles
def create_density_histogram(gas,t_end, filename):
    #gas_positions=numpy.array([gas.x.value_in(units.kpc),gas.y.value_in(units.kpc),gas.z.value_in(units.kpc)])
    # particle_mass = 1.0  # Adjust based on your data
    # volume = numpy.prod(numpy.max(gas_positions, axis=0) - numpy.min(gas_positions, axis=0))
    # densities = len(gas_positions) / volume
    # volume_per_particle = 1.0  # Adjust based on your simulation parameters
    # densities = numpy.zeros(gas_positions.shape[0])
    # for i in range(gas_positions.shape[0]):
    #     densities[i] = gas_positions[i].sum() / volume_per_particle

    # # Create a k-d tree for efficient neighbor search
    # kdtree = cKDTree(gas_positions)

    # # Set a smoothing length for the Gaussian kernel
    # smoothing_length = 0.1  # You may need to adjust this based on your simulation parameters

    # # Calculate density for each particle
    # densities = numpy.zeros(gas_positions.shape[0])


    # for i in range(gas_positions.shape[0]):
    #     distances, neighbors = kdtree.query(gas_positions[i], k=32)  # Adjust 'k' based on the number of neighbors
    #     weights = numpy.exp(-0.5 * (distances / smoothing_length)**2)
    #     densities[i] = numpy.sum(weights) / (numpy.pi * smoothing_length**2)
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
def density_scatter( stars ,gas,filename,sort = True, bins = 20,save=False, **kwargs )   :

    """
    Scatter plot colored by 2d histogram
    """
    gasx,gasy=gas.x.value_in(units.kpc), gas.y.value_in(units.kpc)
    starsx,starsy=stars.x.value_in(units.kpc), stars.y.value_in(units.kpc)
    fig,ax= plt.subplots()
    
    data , x_e, y_e = numpy.histogram2d(gasx, gasy, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , numpy.vstack([gasx,gasy]).T , method = "splinef2d", bounds_error = False)
    #print(numpy.nanmin(z),numpy.nanmax(z))
    
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

    #norm = Normalize(vmin = numpy.min(z), vmax = numpy.max(z))
    norm = Normalize(vmin = 0, vmax = 1)
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm),ax=ax)
    cbar.ax.set_ylabel('Density')
    fig.show()
    if save:
        fig.savefig(filename)  

num_fig = 60
t_endpoint = 3000 | units.Myr
t_ends = numpy.linspace(t_endpoint.value_in(units.Myr)/num_fig,t_endpoint.value_in(units.Myr),num_fig) |units.Myr
for t_end in t_ends:
    if t_end> 2200 | units.Myr:
        break
    gas, stars = read_from_file(t_end)
    create_density_histogram(gas,t_end, filename="/home/lea/Amuse-env/pictures/flyby_histogram_t"+str(int(t_end.value_in(units.Myr)))+"Myr.png")
    #   density_scatter(stars, gas,  filename="/home/lea/Amuse-env/pictures/densities_scatter_t"+str(int(t_end.value_in(units.Myr)))+"Myr.png", save=True)
#     make_plot(stars, gas, filename="/home/lea/Amuse-env/pictures/Galaxy_merger_t"+str(int(t_end.value_in(units.Myr)))+"Myr.png")

# gas,stars = read_from_file(2200.0| units.Myr)
# # create_density_histogram(stars, filename="/home/lea/Amuse-env/pictures/densities_final.png")
# # plt.close()
# # gas,stars = read_from_file(100.0| units.Myr)
# create_density_histogram(stars,2200.0| units.Myr, filename="/home/lea/Amuse-env/pictures/densities_fast_vel2.png")

# density_scatter(stars, gas,  filename="/home/lea/Amuse-env/pictures/densities_scatter_fast_vel2.png", save=True)
