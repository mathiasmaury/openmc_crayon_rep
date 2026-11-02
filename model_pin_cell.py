import numpy as np
import openmc

def build_model():
    """
    Construit et retourne la géométrie, les matériaux, les settings et les tallies de la cellule REP.
    - Conditions 'reflective' en X et Y pour simuler le réseau infini d'assemblages.
    - Conditions 'vacuum' aux extrémités Z (-20 cm à +20 cm) pour créer le profil axial.
    """
    # =========================================================================
    # 1. MATÉRIAUX
    # =========================================================================
    
    # Combustible : UO2 enrichi à 3.5% en U235, densité 10.4 g/cm3
    fuel = openmc.Material(name='UO2')
    fuel.set_density('g/cm3', 10.4)
    fuel.add_element('U', 1, enrichment=3.5)
    fuel.add_element('O', 2)
    fuel.temperature = 900
    
    # Spécification du volume pour le calcul de déplétion (cm3)
    # Hauteur axiale h_val = 40.0 cm (-20 cm à +20 cm)
    r_fuel_val = 0.4096
    h_val = 40.0
    fuel.volume = np.pi * (r_fuel_val**2) * h_val

    # Gaine (cladding) : Zircaloy-4, densité 6.55 g/cm3
    clad = openmc.Material(name='Zircaloy-4')
    clad.set_density('g/cm3', 6.55)
    clad.add_element('Zr', 1.0)
    clad.temperature = 300

    # Modérateur / caloporteur : eau claire
    water = openmc.Material(name='water')
    water.set_density('g/cm3', 0.74)
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    water.temperature = 300

    materials = openmc.Materials([fuel, clad, water])

    # =========================================================================
    # 2. GÉOMÉTRIE (RÉSEAU INFINI RADIAL X/Y + FUITES AXIALES Z)
    # =========================================================================
    
    r_fuel = openmc.ZCylinder(r=r_fuel_val)
    r_clad = openmc.ZCylinder(r=0.4750)
    pitch = 1.26
    h_z = h_val / 2.0  # 20.0 cm

    # Conditions réfléchissantes en X et Y (réseau infini)
    min_x = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
    max_x = openmc.XPlane(x0=pitch/2, boundary_type='reflective')
    min_y = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
    max_y = openmc.YPlane(y0=pitch/2, boundary_type='reflective')

    # Conditions de vide en Z pour créer les fuites axiales et le profil en cosinus
    min_z = openmc.ZPlane(z0=-h_z, boundary_type='vacuum')
    max_z = openmc.ZPlane(z0=h_z, boundary_type='vacuum')

    # Régions géométriques
    fuel_region = -r_fuel & +min_z & -max_z
    clad_region = +r_fuel & -r_clad & +min_z & -max_z
    water_region = +r_clad & +min_x & -max_x & +min_y & -max_y & +min_z & -max_z

    # Cellules
    fuel_cell = openmc.Cell(name='Fuel', fill=fuel, region=fuel_region)
    clad_cell = openmc.Cell(name='Clad', fill=clad, region=clad_region)
    water_cell = openmc.Cell(name='Water', fill=water, region=water_region)

    root_universe = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])
    geometry = openmc.Geometry(root_universe)

    # =========================================================================
    # 3. SETTINGS
    # =========================================================================
    
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = 120
    settings.inactive = 20
    settings.particles = 50000

    # Source spatiale uniforme dans le volume du combustible
    bounds = [-r_fuel_val, -r_fuel_val, -h_z, r_fuel_val, r_fuel_val, h_z]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source = openmc.IndependentSource(space=uniform_dist)

    # =========================================================================
    # 4. TALLIES (SPECTRE & MAILLAGE 3D POUR PARAVIEW / PLOTS)
    # =========================================================================
    
    tallies = openmc.Tallies()

    # A. Tally Spectre énergétique par cellule (groupe CCFE-709)
    energy_bins = openmc.mgxs.GROUP_STRUCTURES['CCFE-709']
    energy_filter = openmc.EnergyFilter(energy_bins)
    cell_filter_spec = openmc.CellFilter([fuel_cell, clad_cell, water_cell])
    
    tally_spec = openmc.Tally(name='flux_spectrum')
    tally_spec.filters = [energy_filter, cell_filter_spec]
    tally_spec.scores = ['flux']
    tallies.append(tally_spec)

    # B. Maillage 3D raffiné (10 x 10 x 50 mailles axiales)
    mesh = openmc.RegularMesh()
    mesh.dimension = [10, 10, 50]
    mesh.lower_left = [-pitch/2, -pitch/2, -h_z]
    mesh.upper_right = [pitch/2, pitch/2, h_z]
    mesh_filter = openmc.MeshFilter(mesh)

    # B.1. Tally Flux Global (sur tout le maillage)
    tally_flux_mesh = openmc.Tally(name='flux_mesh_3d')
    tally_flux_mesh.filters = [mesh_filter]
    tally_flux_mesh.scores = ['flux']
    tallies.append(tally_flux_mesh)

    # B.2. Tally Puissance (kappa-fission restreint strictement au combustible)
    fuel_cell_filter = openmc.CellFilter(fuel_cell)
    tally_power_mesh = openmc.Tally(name='power_mesh_3d')
    tally_power_mesh.filters = [mesh_filter, fuel_cell_filter]
    tally_power_mesh.scores = ['kappa-fission']
    tallies.append(tally_power_mesh)

    # Instanciation de l'objet Model d'OpenMC
    model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)

    return {
        'model': model,
        'materials': materials,
        'geometry': geometry,
        'settings': settings,
        'tallies': tallies,
        'fuel_cell': fuel_cell,
        'clad_cell': clad_cell,
        'water_cell': water_cell,
        'fuel': fuel,
        'pitch': pitch
    }

if __name__ == '__main__':
    data = build_model()
    data['model'].export_to_xml()
    print("[Model] Fichiers XML générés avec succès via model_pin_cell.py !")