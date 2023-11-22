import numpy
from matplotlib import pyplot 
from prepare_figure import single_frame
from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver
from hydrodynamics_class import Hydro

timeprint=[]
nprint=[]

def make_plot(gasparticles, starsparticles, filename):
    x_label = "X (pc)"
    y_label = "Y (pc)"
    figure = single_frame(x_label, y_label, logy=False, xsize=14, ysize=14)
    from distinct_colours import get_distinct
    c = get_distinct(3)
    if len(gasparticles)>0:
        pyplot.scatter(gasparticles.x.value_in(units.pc), gasparticles.y.value_in(units.pc),
                   c=c[0], alpha=1, s=1, lw=0)
    if len(starsparticles)>0:
        pyplot.scatter(starsparticles.x.value_in(units.pc), starsparticles.y.value_in(units.pc),
                   c=c[1], alpha=1, s=20, lw=0)
    pyplot.savefig(filename)

def run_molecular_cloud(N=100, Mcloud=100. | units.MSun,
                        Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)
    gas=molecular_cloud(targetN=N,convert_nbody=conv,
            base_grid=body_centered_grid_unit_cube).result
    gas.name = "gas"

    rho_cloud = Mcloud/Rcloud**3
    tff = 0.5427/numpy.sqrt(constants.G*rho_cloud)

    stars = Particles(0)
    hydro = Hydro(Fi, gas, stars)

    rho_cloud = 3.*Mcloud/(4.*numpy.pi*Rcloud**3)

    dt = 0.02 | units.Myr
    tend= 1.0 | units.Myr
    dt_diag = 0.1 | units.Myr
    t_diag = 0 | units.Myr

    i=0
    E0 = 0.0
    time = 0.0 | units.Myr

    while time < tend:
        time += dt
        print("Evolve to time=", time.in_(units.Myr))
        Mtot = 0|units.MSun
        if len(hydro.star_particles) > 0:
            Mtot = hydro.gas_particles.mass.sum() \
                    + hydro.star_particles.mass.sum()
        else:
            Mtot = hydro.gas_particles.mass.sum()

        if Mtot < Mcloud-(1.e-5|units.MSun):
            exit(-1)

        hydro.evolve_model(time)
        E = hydro.gas_particles.kinetic_energy() \
             + hydro.gas_particles.potential_energy() \
             + hydro.gas_particles.thermal_energy()
        E_th = hydro.gas_particles.thermal_energy()
        if i==0:
            E0 = E
        Eerr = (E-E0)/E0

        hydro.print_diagnostics()
        if time>t_diag:
            t_diag += dt_diag
            hydro.write_set_to_file(i)
            make_plot(hydro.gas_particles, hydro.star_particles,"gas_cloud_"+str(i))
            i=i+1
            
              

        timeprint.append(time.value_in(units.Myr))
        nprint.append(hydro.get_number_of_star_particles())
    hydro.stop()
    return gas

if __name__ in ("__main__","__plot__"):
    numpy.random.seed(3141)
    parts = run_molecular_cloud(1e5, Mcloud=1e6 | units.MSun,
                                Rcloud=30 | units.parsec)
    print(timeprint,nprint)
    pyplot.plot(timeprint,nprint)
    pyplot.savefig("stars_over_time.png")

