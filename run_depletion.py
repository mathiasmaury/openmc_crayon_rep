import os
import openmc
import openmc.deplete
from model_pin_cell import build_model

def run_depletion_simulation():
    print("[Depletion] Importation du modèle central...")
    model_data = build_model()
    
    # Ajustement des settings Monte Carlo pour la déplétion
    settings = model_data['settings']
    settings.batches = 50
    settings.inactive = 10
    settings.particles = 10000

    model = model_data['model']

    output_dir = 'outputs/depletion'
    os.makedirs(output_dir, exist_ok=True)

    # Chemin absolu vers le fichier de chaîne de décroissance/transmutation
    chain_file = '/root/openmc_projects/openmc_data/chain_casca.xml'

    print("[Depletion] Initialisation de l'opérateur de déplétion...")
    operator = openmc.deplete.CoupledOperator(model, chain_file)

    # Puissance thermique nominale du crayon : 17.3 kW (17 300 W)
    power = 17300.0  # Watts
    
    # 1 an = 365 jours. 6 pas de 2 mois (~60.8333 jours par pas)
    # 60.833333 jours * 86400 secondes/jour = 5 256 000 secondes
    time_steps = [365.0 / 6.0 * 86400.0] * 6  # 6 pas de 2 mois

    integrator = openmc.deplete.PredictorIntegrator(operator, time_steps, power, timestep_units='s')

    print("[Depletion] Lancement du calcul d'épuisement temporel (1 an / 6 pas de 2 mois)...")
    
    cwd = os.getcwd()
    os.chdir(output_dir)
    try:
        integrator.integrate()
    finally:
        os.chdir(cwd)

    print(f"[Depletion] Calcul de déplétion terminé avec succès ! Résultats stockés dans '{output_dir}'.")

if __name__ == '__main__':
    run_depletion_simulation()