import numpy
import copy
from matplotlib import pyplot
from amuse.lab import *
from amuse.couple import bridge
from amuse.ext.composition_methods import *
from amuse.community.bhtree.interface import BHTree
from amuse.community.fi.interface import Fi
import initial_conditions
import os



def save_to_file(gas_particles, star_particles, index, data_path):
    filename = os.path.join(data_path, "gas_particle_data_i{0:04}.hdf5".format(index))
    filename2 = os.path.join(data_path, "star_particle_data_i{0:04}.hdf5".format(index))
    write_set_to_file(gas_particles, filename, "amuse")
    write_set_to_file(star_particles, filename2, "amuse")


def simulate_merger(stars, converter_stars, gas, converter_gas, t_endpoint, num_fig, data_path):
    dynamics = BHTree(converter_stars) 
    dynamics.parameters.epsilon_squared = (100 | units.parsec)**2 
    dynamics.particles.add_particles(stars)


    hydro = Fi(converter_gas, mode="openmp")

    hydro.parameters.use_hydro_flag = True
    hydro.parameters.radiation_flag = False
    hydro.parameters.gamma = 1
    hydro.parameters.isothermal_flag = True
    hydro.parameters.integrate_entropy_flag = False
    hydro.parameters.timestep = 1 | units.Myr 
    hydro.parameters.verbosity = 0
    hydro.parameters.eps_is_h_flag = False 
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

    save_to_file(gas, stars, 0, data_path)
    t_ends=numpy.linspace(t_endpoint.value_in(units.Myr)/num_fig,t_endpoint.value_in(units.Myr),num_fig) |units.Myr

    for t_end in t_ends:
        print("is evolving to:{} Myr".format(t_end.value_in(units.Myr)) )
        gravhydro.evolve_model(t_end)
        channel["to_stars"].copy()
        channel["to_gas"].copy()
        save_to_file(gas, stars, t_end.value_in(units.Myr), data_path)

    dynamics.stop()
    hydro.stop()



def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--M_dm", unit=units.MSun,
                      dest="M_stars", default = 1.0e9 | units.MSun,
                      help="Galaxy dark matter and stellar mass [%default]")
    result.add_option("--M_gas", unit=units.MSun,
                      dest="M_gas", default = 1.0e8 | units.MSun,
                      help="Galaxy gas mass [%default]")
    result.add_option("-R", unit=units.kpc,
                      dest="R_galaxy", default = 3 | units.kpc,
                      help="Galaxy size [%default]")
    result.add_option("--n_bulge", dest="n_bulge", default = 20000,
                      help="number of star particles in the bulge [%default]")
    result.add_option("--n_cloud", dest="n_cloud", default = 10000,
                      help="number of gas particles[%default]")
    result.add_option("--t_end", unit=units.Myr,
                      dest="t_end", default = 3000|units.Myr,
                      help="End of the simulation [%default]")
    result.add_option("--num_fig", dest="num_fig", default = 60,
                      help="Number of data snapshots made [%default]")
    result.add_option("--pos1", unit=units.kpc,
                      dest="galaxy1_position", default = [25.0, 25.0, 0] | units.kpc,
                      help="Position vector of first galaxy [%default]")
    result.add_option("--pos2", unit=units.kpc,
                      dest="galaxy2_position", default = [-25.0, 0, 0] | units.kpc,
                      help="Position vector of second galaxy [%default]")
    result.add_option("-v", unit=units.km/units.s,
                      dest="galaxy_velocity", default = [0.0, 0.0, 0.0] | units.km/units.s,
                      help="Velocity vector of galaxies [%default]")
    result.add_option("-p", dest="data_path", default =  "../../data/interim/",
                      help="Path where data should be saved")
    return result

if __name__ == '__main__':
    o, arguments  = new_option_parser().parse_args()
    stars, converter_stars = initial_conditions.make_galaxy_pair(o.M_galaxy, o.R_galaxy,
                                                o.n_bulge, o.galaxy1_position, o.galaxy2_position, o.galaxy_velocity)
    gas, converter_gas = initial_conditions.make_cloud_pair(o.M_gas, o.R_galaxy,
                                                o.n_cloud, o.galaxy1_position, o.galaxy2_position, o.galaxy_velocity)    
    simulate_merger(stars, converter_stars, gas, converter_gas, o.t_end, o.num_fig, o.data_path)
