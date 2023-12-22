from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model
import numpy
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

def make_galaxy_pair(M_galaxy, R_galaxy, n_bulge, galaxy1_position, galaxy2_position, galaxy_velocity):
    '''
    Creates two dwarf galaxies each made up of n_bulge particles with initial mass M_galaxy and radius R_galaxy, 
    initial positions are galaxy1_position and galaxy2_position, and they have opposing velocities of galaxy_velocity

    Returns dm_particles corresponding to those conditions as well as an appropriate converter
    '''
    converter=nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    galaxy1 = new_galactics_model(0,
                                  converter,
                                  bulge_number_of_particles=n_bulge,
                                  disk_number_of_particles=0)
    galaxy2 = Particles(len(galaxy1))
    galaxy2.mass = galaxy1.mass
    galaxy2.position = galaxy1.position
    galaxy2.velocity = galaxy1.velocity

    galaxy1.rotate(0., numpy.pi/2, numpy.pi/4)
    galaxy1.position += galaxy1_position 
    galaxy1.velocity += galaxy_velocity 

    galaxy2.rotate(numpy.pi/4, numpy.pi/4, 0.0)
    galaxy2.position += galaxy2_position 
    galaxy2.velocity -= galaxy_velocity
    stars = Particles(0)
    stars.add_particles(galaxy1)
    stars.add_particles(galaxy2)

    return stars, converter

def make_cloud_pair(Mcloud, Rcloud, n_cloud, galaxy1_position, galaxy2_position, galaxy_velocity):
    '''
    Creates the gas particles of two dwarf galaxies each made up of n_cloud particles with mass Mcloud and radius Rcloud, 
    initial positions are galaxy1_position and galaxy2_position, and they have opposing velocities of galaxy_velocity

    Returns gas_particles corresponding to those conditions as well as an appropriate converter
    '''
    converter = nbody_system.nbody_to_si(Mcloud,Rcloud)
    gas=molecular_cloud(targetN=n_cloud,convert_nbody=converter,
                    base_grid=body_centered_grid_unit_cube).result
    
    gas2=Particles(len(gas))
    gas2.mass = gas.mass
    gas2.position = gas.position
    gas2.velocity = gas.velocity

    gas.position += galaxy1_position 
    gas.velocity += galaxy_velocity 
    gas2.position += galaxy2_position 
    gas2.velocity -= galaxy_velocity

    gas.add_particles(gas2)
    gas.name = "gas"
    return gas, converter

# def make_galaxies_with_gas(M_galaxy, R_galaxy, n_bulge, Mcloud, Rcloud, n_cloud, galaxy1_position=[25.0, 25, 0] | units.kpc, galaxy2_position=[-25.0, 0, 0] | units.kpc, galaxy_velocity=[0.0, 0.0, 0.0] | units.km/units.s):
#     stars, converter_stars = make_galaxy_pair(M_galaxy, R_galaxy, n_bulge, galaxy1_position, galaxy2_position, galaxy_velocity)
#     gas, converter_gas = make_cloud_pair(Mcloud, Rcloud, n_cloud, galaxy1_position, galaxy2_position, galaxy_velocity)

    