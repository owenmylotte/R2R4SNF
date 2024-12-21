import numpy as np
import openmc
import openmc.deplete
import time

# Core parameters
fueltemp = 900  # K
lbptemp = 600   # K
modtemp = 600   # K
refltemp = 600  # K
bartemp = 600   # K
cooltemp = 600  # K

# Add this line:
N_holes = 216  # Number of fuel holes per fuel block

R = 8.314  # m3Pa /K/mol
T = cooltemp
M = 4  # g/mol
P = 6.39e6  # Pa
rho = M * P / (R*T) / 1e6  # g/cm3

fuelpf = 0.35
new_fuelpf = 0.60  # Increase the packing fraction for FCM-HTGR
ratio_pf = new_fuelpf/fuelpf

# Material definitions
# Path to the cross section file
openmc.Materials.cross_sections = '/Users/owen/PycharmProjects/R2R4SNF/inputs/openmc/nuclear_data/endfb71_hdf5/cross_sections.xml'

kernel = openmc.Material(name='kernel')  # UN fuel
kernel.set_density('g/cm3', 14.3)
kernel.add_nuclide('U235', 0.1975)
kernel.add_nuclide('U238', 0.8025)
kernel.add_nuclide('N15', 0.99)
kernel.add_nuclide('N14', 0.01)
kernel.temperature = fueltemp

buff = openmc.Material(name='buffer')
buff.set_density('g/cm3', 1.0)
buff.add_nuclide('C0', 1.0)
# buff.add_s_alpha_beta('c_Graphite')
buff.temperature = fueltemp

iPyC = openmc.Material(name='iPyC')
iPyC.set_density('g/cm3', 1.9)
iPyC.add_nuclide('C0', 1.0)
# iPyC.add_s_alpha_beta('c_Graphite')
iPyC.temperature = fueltemp

SiC = openmc.Material(name='SiC')
SiC.set_density('g/cm3', 3.2)
SiC.add_nuclide('C0', 0.5)
SiC.add_element('Si', 0.5)
SiC.temperature = fueltemp

oPyC = openmc.Material(name='oPyC')
oPyC.set_density('g/cm3', 1.9)
oPyC.add_nuclide('C0', 1.0)
# oPyC.add_s_alpha_beta('c_Graphite')
oPyC.temperature = fueltemp

compact_matrix = openmc.Material(name='compact_matrix')  # change the compact matrix to SiC
compact_matrix.set_density('g/cm3', 3.2)
compact_matrix.add_nuclide('C0', 0.5)
compact_matrix.add_element('Si', 0.5)
compact_matrix.temperature = fueltemp

compact_graphite = openmc.Material(name='compact_graphite')
compact_graphite.set_density('g/cm3', 1.6479575)
compact_graphite.add_nuclide('C0', 0.999999540, percent_type='wo')
compact_graphite.add_nuclide('B10', 0.000000150, percent_type='wo')
compact_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
compact_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
compact_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
compact_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
compact_graphite.add_nuclide('K39', 0.000000040, percent_type='wo')
compact_graphite.add_nuclide('V51', 0.000000018, percent_type='wo')
# compact_graphite.add_s_alpha_beta('c_Graphite')
compact_graphite.temperature = fueltemp

block_graphite = openmc.Material(name='block_graphite')
block_graphite.set_density('g/cm3', 1.85)
block_graphite.add_nuclide('C0', 0.999999540, percent_type='wo')
block_graphite.add_nuclide('B10', 0.000000150, percent_type='wo')
block_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
block_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
block_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
block_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
block_graphite.add_nuclide('K39', 0.000000040, percent_type='wo')
block_graphite.add_nuclide('V51', 0.000000018, percent_type='wo')
# block_graphite.add_s_alpha_beta('c_Graphite')
block_graphite.temperature = modtemp

coolant = openmc.Material(name='coolant')
coolant.set_density('g/cm3', 0.02502)
coolant.add_nuclide('Na23', 1.0)
coolant.temperature = cooltemp

reflector_graphite = openmc.Material(name='reflector_graphite')
reflector_graphite.set_density('g/cm3', 1.85)
reflector_graphite.add_nuclide('C0', 0.999999540, percent_type='wo')
reflector_graphite.add_nuclide('B10', 0.000000150, percent_type='wo')
reflector_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
reflector_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
reflector_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
reflector_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
reflector_graphite.add_nuclide('K39', 0.000000040, percent_type='wo')
reflector_graphite.add_nuclide('V51', 0.000000018, percent_type='wo')
# reflector_graphite.add_s_alpha_beta('c_Graphite')
reflector_graphite.temperature = refltemp

permanent_graphite = openmc.Material(name='permanent_graphite')
permanent_graphite.set_density('g/cm3', 1.85)
permanent_graphite.add_nuclide('C0', 0.999999540, percent_type='wo')
permanent_graphite.add_nuclide('B10', 0.000000150, percent_type='wo')
permanent_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
permanent_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
permanent_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
permanent_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
permanent_graphite.add_nuclide('K39', 0.000000040, percent_type='wo')
permanent_graphite.add_nuclide('V51', 0.000000018, percent_type='wo')
# permanent_graphite.add_s_alpha_beta('c_Graphite')
permanent_graphite.temperature = refltemp

barrel = openmc.Material(name='barrel')
barrel.set_density('g/cm3', 7.94)
barrel.add_element('Ni', 32.5)
barrel.add_element('Cr', 21)
barrel.add_element('Fe', 39.5)
barrel.add_element('Al', 0.3)
barrel.add_element('Ti', 0.3)
barrel.add_nuclide('C0', 0.08)
barrel.temperature = bartemp

# Creating homogeneous TRISO particles
hete_triso_radius = [2125e-5, 3125e-5, 3475e-5, 3825e-5, 4225e-5]
hete_lbp_radius = [1000e-5, 1180e-5, 1410e-5]

r0 = 0.5  # fix ratio of each layer equally except kernel
r1 = 0.583  # ratio of compact_matrix
r2 = 0.583  # ratio of lbp_matrix

num_triso = 7091
num_lbp = 52742

# Calculate heterogeneous volumes
hete_kernel = np.pi*(4*hete_triso_radius[0]**3/3)*num_triso*15*10*ratio_pf
hete_buffer = np.pi*(4*(hete_triso_radius[1]**3-hete_triso_radius[0]**3)/3)*num_triso*15*10*ratio_pf
hete_iPyC = np.pi*(4*(hete_triso_radius[2]**3-hete_triso_radius[1]**3)/3)*num_triso*15*10*ratio_pf
hete_SiC = np.pi*(4*(hete_triso_radius[3]**3-hete_triso_radius[2]**3)/3)*num_triso*15*10*ratio_pf
hete_oPyC = np.pi*(4*(hete_triso_radius[4]**3-hete_triso_radius[3]**3)/3)*num_triso*15*10*ratio_pf
hete_PyC = hete_iPyC + hete_oPyC

hete_triso = hete_kernel + hete_buffer + hete_PyC + hete_SiC
hete_matrix = (np.pi*(0.635**2)*793) - hete_triso

# Calculate TRISO radii
triso_radius = []
R0 = ((hete_matrix*r1)/793/np.pi)**(0.5)
R1 = ((hete_matrix*r1 + hete_oPyC*r0)/793/np.pi)**0.5
R2 = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0)/793/np.pi)**0.5
R3 = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC*r0)/793/np.pi)**0.5
R4 = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC*r0 + hete_buffer*r0)/793/np.pi)**0.5
R5 = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC*r0 + hete_buffer*r0 + hete_kernel)/793/np.pi)**0.5
R6 = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC*r0 + hete_buffer + hete_kernel)/793/np.pi)**0.5
R7 = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC + hete_buffer + hete_kernel)/793/np.pi)**0.5
R8 = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC + hete_iPyC + hete_buffer + hete_kernel)/793/np.pi)**0.5
R9 = ((hete_matrix*r1 + hete_oPyC + hete_SiC + hete_iPyC + hete_buffer + hete_kernel)/793/np.pi)**0.5
R10 = ((hete_matrix + hete_oPyC + hete_SiC + hete_iPyC + hete_buffer + hete_kernel)/793/np.pi)**0.5

triso_radius.extend([R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10])
print("\ntriso_radius :")
print(triso_radius)

# Creating TRISO particles
triso_spheres = [openmc.model.Cylinder(r=r) for r in triso_radius]

# Create a rectangular outer boundary with reflective conditions
rad1 = 20.6756  # Based on the fuel block edge length
xmin_plane = openmc.XPlane(x0=-rad1, boundary_type='reflective')
xmax_plane = openmc.XPlane(x0=rad1, boundary_type='reflective')
ymin_plane = openmc.YPlane(y0=-rad1, boundary_type='reflective')
ymax_plane = openmc.YPlane(y0=rad1, boundary_type='reflective')
fuel_block_z0 = openmc.ZPlane(z0=-1/2, boundary_type='reflective')
fuel_block_z1 = openmc.ZPlane(z0=1/2, boundary_type='reflective')

# Define the outer box region
outer_box_region = (+xmin_plane & -xmax_plane &
                   +ymin_plane & -ymax_plane &
                   +fuel_block_z0 & -fuel_block_z1)

# Create TRISO cells
old_triso_cells = []
old_triso_univ = []
for i in range(1):
    old_triso_cells.append([
        openmc.Cell(fill=compact_matrix, region=-triso_spheres[0] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=oPyC, region=+triso_spheres[0] & -triso_spheres[1] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=SiC, region=+triso_spheres[1] & -triso_spheres[2] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=iPyC, region=+triso_spheres[2] & -triso_spheres[3] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=buff, region=+triso_spheres[3] & -triso_spheres[4] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=kernel, region=+triso_spheres[4] & -triso_spheres[5] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=buff, region=+triso_spheres[5] & -triso_spheres[6] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=iPyC, region=+triso_spheres[6] & -triso_spheres[7] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=SiC, region=+triso_spheres[7] & -triso_spheres[8] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=oPyC, region=+triso_spheres[8] & -triso_spheres[9] & +fuel_block_z0 & -fuel_block_z1),
        openmc.Cell(fill=compact_matrix, region=+triso_spheres[9] & -triso_spheres[10] & +fuel_block_z0 & -fuel_block_z1)
    ])
    old_triso_univ.append(openmc.Universe(cells=old_triso_cells[i]))

# Create TRISO cells and universes
triso_cells = []
triso_outer_cells = []
triso_univ = []
for i in range(1):
    triso_cells.append(openmc.Cell(fill=old_triso_univ[i], region=-triso_spheres[-1] & +fuel_block_z0 & -fuel_block_z1))
    triso_outer_cells.append(openmc.Cell(fill=block_graphite, region=+triso_spheres[-1] | -fuel_block_z0 | +fuel_block_z1))
    triso_univ.append(openmc.Universe(cells=[triso_cells[i], triso_outer_cells[i]]))

# Create the hexagonal fuel block
fuel_block = openmc.model.hexagonal_prism(edge_length=20.6756, orientation='x')
fuel_block_region = fuel_block & outer_box_region

# Create coolant channels
small_coolant = openmc.model.Cylinder(r=0.635)
small_coolant_cells = [
    openmc.Cell(fill=coolant, region=-small_coolant),
    openmc.Cell(fill=block_graphite, region=+small_coolant)
]
small_coolant_univ = openmc.Universe(cells=small_coolant_cells)

large_coolant = openmc.model.Cylinder(r=0.794)
large_coolant_cells = [
    openmc.Cell(fill=coolant, region=-large_coolant),
    openmc.Cell(fill=block_graphite, region=+large_coolant)
]
large_coolant_univ = openmc.Universe(cells=large_coolant_cells)

hole = openmc.model.Cylinder(r=0.635)
hole_cells = [
    openmc.Cell(fill=block_graphite, region=-hole),
    openmc.Cell(fill=block_graphite, region=+hole)
]
hole_univ = openmc.Universe(cells=hole_cells)

# Create lattice universes
a = triso_univ[0]
b = large_coolant_univ
c = small_coolant_univ
d = hole_univ

# Create fuel block lattice
fuel_block_lattice = openmc.HexLattice(name='fuel_block_lattice')
fuel_block_lattice.orientation = 'x'
fuel_block_lattice.center = (0, 0)
fuel_block_lattice.pitch = [1.8796]
fuel_block_lattice.universes = [
    [a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a],
    [b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a],
    [a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b],
    [a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a],
    [b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a],
    [a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b],
    [a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a],
    [b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a],
    [a]+[c]+[a]+[c]+[a]+[c]+[a]+[c]+[a]+[c]+[a]+[c],
    [d]*6,
    [d]
]

# Create an outer cell filled with reflector material
outer_cell = openmc.Cell(fill=block_graphite, region=outer_box_region & ~fuel_block)

# Create the fuel block cells and universe
fuel_block_cells = openmc.Cell(name='fuel_block_cells',
                              fill=fuel_block_lattice,
                              region=fuel_block_region)

# Create a universe for outside the lattice
fuel_block_outer_cell = openmc.Cell(fill=block_graphite)
fuel_block_outer_univ = openmc.Universe(cells=[fuel_block_outer_cell])

# Set the outer universe for the lattice
fuel_block_lattice.outer = fuel_block_outer_univ

# Modify the fuel block cells definition to use the outer region
fuel_block_cells = openmc.Cell(
    name='fuel_block_cells',
    fill=fuel_block_lattice,
    region=fuel_block_region
)

# Create the root universe
root_universe = openmc.Universe(cells=[fuel_block_cells, outer_cell])

# Create and export the geometry
geom = openmc.Geometry(root_universe)
geom.export_to_xml()

# Export materials
mats = list(geom.get_all_materials().values())
openmc.Materials(mats).export_to_xml()

# Settings
settings = openmc.Settings()
settings.batches = 50
settings.inactive = 10
settings.particles = 500

# Create uniform spatial source with proper 3D bounds
lower_left = [-rad1, -rad1, -0.5]  # Adding z-coordinate
upper_right = [rad1, rad1, 0.5]    # Adding z-coordinate
uniform_dist = openmc.stats.Box(lower_left, upper_right)
settings.source = openmc.Source(space=uniform_dist)

# Depletion settings
kernel.volume = N_holes * ((np.pi * (triso_spheres[5].r**2 - triso_spheres[4].r**2) * 1.0)) * (19 + 15) + ((np.pi * (triso_spheres[5].r**2 - triso_spheres[4].r**2) * 1.0)) * 12

# Create model
model = openmc.model.Model(geometry=geom, materials=mats, settings=settings)
model.export_to_xml()

# Set up depletion calculation
chain_file = '../nuclear_data/chain_endfb71_pwr.xml'
operator = openmc.deplete.Operator(model, chain_file, diff_burnable_mats=False)

# Calculate power and time steps
power_density = 5  # W/cm3
linear_power = power_density * (np.pi * (rad1**2))  # W/cm
BU_final = 80  # MWd/kgU
N_days = (BU_final * operator.heavy_metal * 1000) / linear_power
end_step = N_days - 3803

# Define time steps
time_operation = ([3/24]*4) + ([4/24]*3) + ([6/24]*4) + ([12/24]*2) + ([1]) + ([3]*1) + \
                ([7]*3) + ([14]*6) + ([28]*8) + ([29]) + ([31]*2) + ([31]) + ([30]) + \
                ([31]) + ([30]) + ([31]) + ([31]) + ([30]) + ([31]) + ([59]) + \
                ([61]*2) + ([62]) + ([61]*2) + ([60]) + ([61]*2) + ([62]) + ([61]*2) + \
                ([120]*20) + ([end_step])

time_cooling = [0.04, 1.0, 7.0, 30.0, 30.0, 300.0] + ([365]*4) + \
              [36.5, 146.0, 91.25, 91.25, 730, 730, 547.5, 1277.5, 1825, 1095, 730] + \
              [1825]*5 + [3650]*5 + [9125]*6 + [18250]*5 + [36500]*5 + [91250]*6 + \
              [182500] + [365000]*7 + [3650000] + [10950000, 18250000, 36500000, 109500000, 182500000]

time_steps = time_operation + time_cooling
power = [linear_power] * len(time_operation) + [0] * len(time_cooling)


# Plot settings
allcolors = {
    buff: 'blue', iPyC: 'yellow', SiC: 'green', oPyC: 'yellow',
    compact_matrix: 'cyan', compact_graphite: 'purple', block_graphite: 'yellow',
    coolant: 'lightblue', permanent_graphite: 'green',
    barrel: 'black', kernel: 'red', reflector_graphite: 'grey'
}

# Create plots
plot1 = openmc.Plot.from_geometry(geom)
plot1.filename = 'core000'
plot1.origin = [0, 0, 0]
plot1.width = [250, 250]
plot1.pixels = (10500, 10500)
plot1.color_by = 'material'
plot1.basis = 'xy'
plot1.colors = allcolors

plot2 = openmc.Plot.from_geometry(geom)
plot2.filename = 'core008'
plot2.origin = [0, 17.905594838485538*2, 0]
plot2.width = [17.905594838485538*2+1, 17.905594838485538*2+1]
plot2.pixels = (8000, 8000)
plot2.color_by = 'material'
plot2.basis = 'xy'
plot2.colors = allcolors

# Export plots
plots = openmc.Plots([plot1, plot2])
plots.export_to_xml()

# Create integrator and run depletion
integrator = openmc.deplete.CECMIntegrator(operator, time_steps, power, timestep_units='d')

print("Running depletion calculation...")
print(f"Total simulation days: {N_days}")
print(f"Operation days: {sum(time_operation)}")

start_time = time.time()
integrator.integrate()
print(f"Calculation completed in {time.time() - start_time:.2f} seconds")
print(f"Heavy metal mass: {operator.heavy_metal}")
