from amuse.lab import read_set_from_file
from amuse.units import units
from making_channel_hydro_stars import make_plot
import numpy

def read_from_file(index):
    filename = "data/gas_particle_data_i{0:04}.hdf5".format(index.value_in(units.Myr))
    filename2 = "data/star_particle_data_i{0:04}.hdf5".format(index.value_in(units.Myr))
    gas_particles = read_set_from_file(filename, format="amuse")
    star_particles = read_set_from_file(filename2, format="amuse")
    return gas_particles, star_particles

num_fig = 30
t_endpoint = 3000 | units.Myr
t_ends = numpy.linspace(t_endpoint.value_in(units.Myr)/num_fig,t_endpoint.value_in(units.Myr),num_fig) |units.Myr
for t_end in t_ends:
    gas, stars = read_from_file(t_end)
    make_plot(stars, gas, filename="pictures/Galaxy_merger_t"+str(int(t_end.value_in(units.Myr)))+"Myr.png")

