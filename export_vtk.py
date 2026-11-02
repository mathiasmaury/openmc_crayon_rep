import os
import numpy as np
import openmc

def export_geometry_vtk(vtk_output_path):
    """
    Génère un fichier VTK représentant le champ 3D des matériaux (Géométrie).
    - Material 1 : Combustible UO2
    - Material 2 : Gaine Zircaloy-4
    - Material 3 : Modérateur H2O
    """
    pitch = 1.26
    h_z = 20.0  # -20 cm à +20 cm
    r_fuel = 0.4096
    r_clad = 0.4750

    # Résolution de la grille 3D pour la géométrie (10 x 10 x 50)
    nx, ny, nz = 10, 10, 50
    x = np.linspace(-pitch/2, pitch/2, nx + 1)
    y = np.linspace(-pitch/2, pitch/2, ny + 1)
    z = np.linspace(-h_z, h_z, nz + 1)

    # Centres des mailles
    x_c = 0.5 * (x[:-1] + x[1:])
    y_c = 0.5 * (y[:-1] + y[1:])
    
    # Création du champ de matériaux
    mat_ids = np.zeros((nx, ny, nz), dtype=int)
    for ix in range(nx):
        for iy in range(ny):
            r = np.sqrt(x_c[ix]**2 + y_c[iy]**2)
            if r <= r_fuel:
                mat_id = 1  # Combustible
            elif r <= r_clad:
                mat_id = 2  # Gaine
            else:
                mat_id = 3  # Modérateur
            mat_ids[ix, iy, :] = mat_id

    os.makedirs(os.path.dirname(vtk_output_path), exist_ok=True)

    with open(vtk_output_path, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('Geometrie Materiaux Crayon REP\n')
        f.write('ASCII\n')
        f.write('DATASET RECTILINEAR_GRID\n')
        f.write(f'DIMENSIONS {nx+1} {ny+1} {nz+1}\n')
        
        f.write(f'X_COORDINATES {nx+1} float\n')
        f.write(' '.join(map(str, x)) + '\n')
        
        f.write(f'Y_COORDINATES {ny+1} float\n')
        f.write(' '.join(map(str, y)) + '\n')
        
        f.write(f'Z_COORDINATES {nz+1} float\n')
        f.write(' '.join(map(str, z)) + '\n')
        
        f.write(f'CELL_DATA {nx*ny*nz}\n')
        f.write('SCALARS material_id int 1\n')
        f.write('LOOKUP_TABLE default\n')
        f.write(' '.join(map(str, mat_ids.flatten(order='F'))) + '\n')


def export_flux_vtk(sp_filename, vtk_output_path):
    """
    Lit le StatePoint statique et écrit un fichier VTK avec le flux et la puissance.
    """
    sp = openmc.StatePoint(sp_filename)
    
    # Extraction du flux
    tally_flux = sp.get_tally(name='flux_mesh_3d')
    mesh_obj = tally_flux.find_filter(openmc.MeshFilter).mesh
    flux = tally_flux.get_values(scores=['flux']).reshape(mesh_obj.dimension[::-1]).transpose(2, 1, 0)
    
    # Extraction de la puissance (kappa-fission dans le combustible)
    tally_power = sp.get_tally(name='power_mesh_3d')
    power = tally_power.get_values(scores=['kappa-fission']).reshape(mesh_obj.dimension[::-1]).transpose(2, 1, 0)

    nx, ny, nz = mesh_obj.dimension
    lx, ly = mesh_obj.lower_left, mesh_obj.upper_right

    x = np.linspace(lx[0], ly[0], nx + 1)
    y = np.linspace(lx[1], ly[1], ny + 1)
    z = np.linspace(lx[2], ly[2], nz + 1)

    os.makedirs(os.path.dirname(vtk_output_path), exist_ok=True)

    with open(vtk_output_path, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('Flux et Puissance OpenMC 3D Statique\n')
        f.write('ASCII\n')
        f.write('DATASET RECTILINEAR_GRID\n')
        f.write(f'DIMENSIONS {nx+1} {ny+1} {nz+1}\n')

        f.write(f'X_COORDINATES {nx+1} float\n')
        f.write(' '.join(map(str, x)) + '\n')

        f.write(f'Y_COORDINATES {ny+1} float\n')
        f.write(' '.join(map(str, y)) + '\n')

        f.write(f'Z_COORDINATES {nz+1} float\n')
        f.write(' '.join(map(str, z)) + '\n')

        f.write(f'CELL_DATA {nx*ny*nz}\n')

        # Champ 1 : Flux
        f.write('SCALARS flux float 1\n')
        f.write('LOOKUP_TABLE default\n')
        f.write(' '.join(map(str, flux.flatten(order='F'))) + '\n')

        # Champ 2 : Puissance
        f.write('SCALARS power float 1\n')
        f.write('LOOKUP_TABLE default\n')
        f.write(' '.join(map(str, power.flatten(order='F'))) + '\n')


def export_all_vtk():
    print("[Export VTK] Génération des fichiers VTK pour le projet final...")
    
    # 1. Génération du VTK de Géométrie (Matériaux)
    geom_vtk_path = "vtk_exports/geometry_3d.vtk"
    export_geometry_vtk(geom_vtk_path)
    print(f"   -> Fichier géométrie généré : {geom_vtk_path}")

    # 2. Génération du VTK de Flux/Puissance Statique
    static_dir = "outputs/static"
    if os.path.exists(static_dir):
        sp_files = [f for f in os.listdir(static_dir) if f.startswith("statepoint")]
        if sp_files:
            latest_sp = os.path.join(static_dir, sorted(sp_files)[-1])
            flux_vtk_path = "vtk_exports/flux_3d_static.vtk"
            export_flux_vtk(latest_sp, flux_vtk_path)
            print(f"   -> Fichier flux/puissance statique généré : {flux_vtk_path}")
            
            # Copie vers le répertoire de partage WSL /mnt/c/ si disponible
            if os.path.exists("/mnt/c/openmc_data"):
                export_geometry_vtk("/mnt/c/openmc_data/geometry_3d.vtk")
                export_flux_vtk(latest_sp, "/mnt/c/openmc_data/flux_3d_static.vtk")
                print("   -> Fichiers également copiés vers '/mnt/c/openmc_data/'.")
        else:
            print("   -> Aucun fichier statepoint trouvé dans 'outputs/static/'. Lance d'abord run_static.py.")
    else:
        print("   -> Dossier 'outputs/static/' non trouvé. Lance d'abord run_static.py.")

if __name__ == '__main__':
    export_all_vtk()