import os
import matplotlib.pyplot as plt
import numpy as np
import openmc
import openmc.deplete

def plot_static_results():
    print("[Plot] Traitement des résultats statiques...")
    
    # 1. Recherche automatique du StatePoint statique
    static_dir = 'outputs/static'
    if not os.path.exists(static_dir):
        print(f"   -> Dossier '{static_dir}' non trouvé. Exécute d'abord run_static.py.")
        return

    sp_files = [f for f in os.listdir(static_dir) if f.startswith('statepoint')]
    if not sp_files:
        print(f"   -> Aucun fichier StatePoint trouvé dans '{static_dir}'.")
        return

    sp_path = os.path.join(static_dir, sorted(sp_files)[-1])
    sp = openmc.StatePoint(sp_path)

    # A. Affichage du keff
    keff = sp.keff
    print(f"   -> k_eff : {keff.nominal_value:.5f} ± {keff.std_dev:.5f}")

    # B. Spectre énergétique en Léthargie (E * Phi(E))
    try:
        tally_spec = sp.get_tally(name='flux_spectrum')
        energy_filter = tally_spec.find_filter(openmc.EnergyFilter)
        energy_bins = energy_filter.bins

        # Énergies centrales géométriques et largeurs dE
        energy_centers = np.sqrt(energy_bins[:, 0] * energy_bins[:, 1]) #centre d'intervalle d'énergie = moyenne géomtrique
        dE = energy_bins[:, 1] - energy_bins[:, 0] #1ere colonne de energy_bins = bornes inférieures d'énergies des groupes, 2e colonne = bornes sup

        # Extraction robuste par CellFilter
        cell_filter = tally_spec.find_filter(openmc.CellFilter) #le filtre sur le type de cellule (fuel, clad ou water)
        cell_ids = cell_filter.bins

        flux_fuel_raw = tally_spec.get_slice(scores=['flux'], filters=[openmc.CellFilter], filter_bins=[(cell_ids[0],)]).mean.ravel() # cell_ids[0] = id de la cellule "fuel"
        flux_clad_raw = tally_spec.get_slice(scores=['flux'], filters=[openmc.CellFilter], filter_bins=[(cell_ids[1],)]).mean.ravel()
        flux_water_raw = tally_spec.get_slice(scores=['flux'], filters=[openmc.CellFilter], filter_bins=[(cell_ids[2],)]).mean.ravel()

        # Conversion en flux par unité de léthargie : (flux / dE) * E
        flux_fuel = (flux_fuel_raw / dE) * energy_centers
        flux_clad = (flux_clad_raw / dE) * energy_centers
        flux_water = (flux_water_raw / dE) * energy_centers

        plt.figure(figsize=(9, 5))
        plt.loglog(energy_centers, flux_fuel, label='Combustible ($UO_2$)', color='firebrick', lw=2)
        plt.loglog(energy_centers, flux_clad, label='Gaine (Zircaloy-4)', color='forestgreen', ls='--', lw=1.5)
        plt.loglog(energy_centers, flux_water, label='Modérateur ($H_2O$)', color='navy', lw=2)
        plt.title('Spectre du flux neutronique en léthargie ($E \cdot \Phi(E)$)')
        plt.xlabel('Énergie [eV]')
        plt.ylabel('Flux par unité de léthargie [arb. units]')
        plt.grid(True, which='both', ls='--', alpha=0.5)
        plt.legend()
        plt.tight_layout()
        plt.savefig('figures/spectre_neutronique.png', dpi=300)
        plt.close()
    except Exception as e:
        print(f"   -> Erreur sur le tracé du spectre : {e}")

    # C. Profil 1D Axial (Z) : Flux uniquement
    try:
        tally_flux = sp.get_tally(name='flux_mesh_3d')
        mesh_flux = tally_flux.find_filter(openmc.MeshFilter).mesh
        flux_3d = tally_flux.get_values(scores=['flux']).reshape(mesh_flux.dimension[::-1]).transpose(2, 1, 0)

        # Moyenne axiale le long de Z (moyenné sur les axes X=0 et Y=1)
        flux_1d_z = flux_3d.mean(axis=(0, 1)) #moyenne sur X et Y pour obtenir le profil axial en Z
        z_coords = np.linspace(mesh_flux.lower_left[2], mesh_flux.upper_right[2], mesh_flux.dimension[2]) #50 points le long de l'axe Z

        plt.figure(figsize=(8, 4))
        plt.plot(z_coords, flux_1d_z, color='tab:blue', lw=2, label='Flux moyen')
        plt.xlabel('Position Z (Hauteur) [cm]')
        plt.ylabel('Flux moyen [n / cm² / s]')
        plt.title('Distribution spatiale 1D Axiale (Z) du flux neutronique')
        plt.grid(True, ls='--', alpha=0.5)
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig('figures/profil_flux_puissance_1d.png', dpi=300)
        plt.close()
    except Exception as e: #si erreur lors de l'execution du bloc try ci dessus
        print(f"   -> Erreur sur le tracé du profil axial (Z) : {e}")


def plot_depletion_results():
    print("[Plot] Traitement des résultats de déplétion...")
    depletion_file = 'outputs/depletion/depletion_results.h5'
    if not os.path.exists(depletion_file):
        print(f"   -> Fichier '{depletion_file}' non trouvé. Exécute d'abord run_depletion.py.")
        return

    results = openmc.deplete.Results(depletion_file)
    
    time_days = results.get_times() #temps en jour
    time_months = time_days / 30.4167 #temps en mois

    # A. Évolution du keff
    _, keff_matrix = results.get_keff()
    keff_values = keff_matrix[:, 0]
    keff_stdev = keff_matrix[:, 1]

    plt.figure(figsize=(8, 4))
    plt.plot(time_months, keff_values, marker='o', color='firebrick', lw=2, label='$k_{eff}$')
    plt.fill_between(time_months, keff_values - keff_stdev, keff_values + keff_stdev, color='firebrick', alpha=0.2)
    plt.title("Chute de réactivité au cours du temps (12 mois / 6 pas)")
    plt.xlabel("Temps [Mois]")
    plt.ylabel("Facteur de multiplication ($k_{eff}$)")
    plt.grid(True, ls='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("figures/keff_depletion.png", dpi=300)
    plt.close()

    # B. Évolution Isotopique
    try:
        first_step = results[0] #pour accéder aux données du premier pas de déplétion, pour inspecter la structure initiale du modèle (notamment la liste des matériaux consommables)
        burnable_mat_id = list(first_step.index_mat.keys())[0] #Récupère l'identifiant numérique (l'ID du matériau OpenMC) du combustible soumis à la déplétion.

        _, atoms_u235 = results.get_atoms(burnable_mat_id, "U235") #tableau des temps ignoré renvoyé en première position
        _, atoms_pu239 = results.get_atoms(burnable_mat_id, "Pu239")

        plt.figure(figsize=(8, 4))
        plt.plot(time_months, atoms_u235, marker='s', label='$^{235}$U (Consommation)', color='navy', lw=2)
        plt.plot(time_months, atoms_pu239, marker='^', label='$^{239}$Pu (Création)', color='forestgreen', lw=2)
        plt.title("Évolution isotopique du combustible au cours du cycle")
        plt.xlabel("Temps [Mois]")
        plt.ylabel("Nombre d'atomes")
        plt.yscale('log')
        plt.grid(True, which='both', ls='--', alpha=0.5)
        plt.legend()
        plt.tight_layout() #ajuste automatiquement les marges et les espacements d'une figure pour que tous les éléments s'y insèrent proprement sans se chevaucher.
        plt.savefig("figures/isotopes_depletion.png", dpi=300)
        plt.close()
    except Exception as e:
        print(f"   -> Erreur sur l'extraction des isotopes : {e}")


if __name__ == '__main__':
    os.makedirs('figures', exist_ok=True)
    plot_static_results()
    plot_depletion_results()
    print("[Plot] Tous les graphiques ont été générés avec succès dans le dossier 'figures/' !")