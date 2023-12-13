import numpy
import copy
from matplotlib import pyplot
from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model
from prepare_figure import single_frame
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver
from amuse.couple import bridge
from amuse.ext.composition_methods import *
# from amuse.community.gadget2.interface import Gadget2
# from amuse.community.hermite.interface import Hermite
from amuse.community.bhtree.interface import BHTree
from amuse.community.fi.interface import Fi



def make_plot(disk1,  gasparticles, filename,disk2=[], starsparticles=[]):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    figure = single_frame(x_label, y_label, logy=False, xsize=14, ysize=14)
    from distinct_colours import get_distinct
    c = get_distinct(4)
    pyplot.xlim(-50, 50)
    pyplot.ylim(-50, 50)

    pyplot.scatter(disk1.x.value_in(units.kpc), disk1.y.value_in(units.kpc),
                   c=c[0], alpha=0.5, s=1, lw=0)
    if len(disk2)>0:
        pyplot.scatter(disk2.x.value_in(units.kpc), disk2.y.value_in(units.kpc),
                   c=c[1], alpha=0.5, s=1, lw=0)
    
    if len(gasparticles)>0:
        pyplot.scatter(gasparticles.x.value_in(units.kpc), gasparticles.y.value_in(units.kpc),
                   c=c[2], alpha=1, s=1, lw=0)
    if len(starsparticles)>0:
        pyplot.scatter(starsparticles.x.value_in(units.kpc), starsparticles.y.value_in(units.kpc),
                   c=c[3], alpha=1, s=1, lw=0)

    pyplot.savefig(filename)

def save_to_file(gas_particles, star_particles, index):
    filename = "data/gas_particle_data_i{0:04}.hdf5".format(index)
    filename2 = "data/star_particle_data_i{0:04}.hdf5".format(index)
    write_set_to_file(gas_particles, filename, "amuse")
    write_set_to_file(star_particles, filename2, "amuse")

def make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk):
    converter=nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    galaxy1 = new_galactics_model(n_halo,
                                  converter,
                                  #do_scale = True,
                                  bulge_number_of_particles=n_bulge,
                                  disk_number_of_particles=n_disk)
    galaxy2 = Particles(len(galaxy1))
    galaxy2.mass = galaxy1.mass
    galaxy2.position = galaxy1.position
    galaxy2.velocity = galaxy1.velocity


    
    galaxy1.rotate(0., numpy.pi/2, numpy.pi/4)
    galaxy1.position += [25.0, 25, 0] | units.kpc
    #galaxy1.velocity += [-3000.0, 0.0, -3000.0] | units.km/units.s
    galaxy1.velocity += [0.0, 0.0, 0.0] | units.km/units.s

    galaxy2.rotate(numpy.pi/4, numpy.pi/4, 0.0)
    galaxy2.position -= [25.0, 0, 0] | units.kpc
    galaxy2.velocity -= [0.0, 0.0, 0] | units.km/units.s

    return galaxy1, galaxy2, converter


def simulate_merger(galaxy1, galaxy2, converter, n_halo, t_endpoint,num_fig=10, 
                    N=100, Mcloud=100. | units.MSun, Rcloud=1. | units.parsec):


    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)
    gas=molecular_cloud(targetN=N,convert_nbody=conv,
                    base_grid=body_centered_grid_unit_cube).result
    
    gas2=Particles(len(gas))
    gas2.mass = gas.mass
    gas2.position = gas.position
    gas2.velocity = gas.velocity

    gas.position += [25.0, 25, 0] | units.kpc
    gas2.position -= [25.0, 0, 0] | units.kpc 

    gas.add_particles(gas2)
    gas.name = "gas"
    stars = Particles(0)
    stars.add_particles(galaxy1)
    stars.add_particles(galaxy2)



    # hydro = Hydro(Fi, gas, stars)
    #hydro.gas_particles.add_particles(gas)

    #converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    dynamics = BHTree(converter) #you can add number of workers here (1 or 2)

    
    dynamics.parameters.epsilon_squared = (100 | units.parsec)**2
    #print(dynamics.parameters)
  
    dynamics.particles.add_particles(stars)
    #set2 = dynamics.particles.add_particles(galaxy2)
    #dynamics.particles.move_to_center()
    #disk1 = set1[:n_halo]
    #disk2 = set2[:n_halo]

    hydro = Fi(conv, mode="openmp")
    #print(hydro.parameters)
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.radiation_flag = False
    hydro.parameters.gamma = 1
    hydro.parameters.isothermal_flag = True
    hydro.parameters.integrate_entropy_flag = False
    hydro.parameters.timestep = 1 | units.Myr 
    hydro.parameters.verbosity = 0
    hydro.parameters.eps_is_h_flag = False    # h_smooth is constant
    eps = 0.1 | units.au
    hydro.parameters.gas_epsilon = eps
    hydro.parameters.sph_h_const = eps

    hydro.particles.add_particles(gas)
    hydro.dm_particles.add_particles(stars.as_set())

    channel = {"from_stars": stars.new_channel_to(dynamics.particles),
            "to_stars": dynamics.particles.new_channel_to(stars)}

    channel.update({"from_gas": gas.new_channel_to(hydro.particles)})
    channel.update({"to_gas": hydro.particles.new_channel_to(gas)})
    channel.update({"from_stars_to_hydro": stars.new_channel_to(hydro.dm_particles)})
    channel.update({"from_hydro_to_stars": hydro.dm_particles.new_channel_to(stars)})


    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(hydro, (dynamics,))
    gravhydro.add_system(dynamics, (hydro,))
    gravhydro.timestep = 10 | units.Myr

    #make_plot(stars, gas,  filename="Galaxy_merger_t0Myr.png")
    save_to_file(gas, stars, 0)
    t_ends=numpy.linspace(t_endpoint.value_in(units.Myr)/num_fig,t_endpoint.value_in(units.Myr),num_fig) |units.Myr
    print(t_ends)
    for t_end in t_ends:
        print("is evolving to:{} Myr".format(t_end.value_in(units.Myr)) )
        gravhydro.evolve_model(t_end)
        channel["to_stars"].copy()
        channel["to_gas"].copy()
        #make_plot(stars,gas, filename="pictures/Galaxy_merger_t"+str(int(t_end.value_in(units.Myr)))+"Myr.png")
        save_to_file(gas, stars, t_end.value_in(units.Myr))

    dynamics.stop()
    hydro.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="M_galaxy", default = 1.0e9 | units.MSun,
                      help="Galaxy mass [%default]")
    result.add_option("-R", unit=units.kpc,
                      dest="R_galaxy", default = 3 | units.kpc,
                      help="Galaxy size [%default]")
    result.add_option("--n_bulge", dest="n_bulge", default = 10000,
                      help="number of stars in the bulge [%default]")
    result.add_option("--n_disk", dest="n_disk", default = 0,
                      help="number of stars in the disk [%default]")
    result.add_option("--n_halo", dest="n_halo", default = 10000,
                      help="number of stars in the halo [%default]")
    result.add_option("--t_end", unit=units.Myr,
                      dest="t_end", default = 3000|units.Myr,
                      help="End of the simulation [%default]")
    return result

if __name__ == '__main__':
    o, arguments  = new_option_parser().parse_args()
    galaxy1, galaxy2, converter = make_galaxies(o.M_galaxy, o.R_galaxy,
                                                o.n_halo, o.n_bulge, o.n_disk)
    
    simulate_merger(galaxy1, galaxy2, converter, o.n_halo, o.t_end,num_fig=30,N=10000, Mcloud=100000. | units.MSun, Rcloud=3000. | units.parsec)
