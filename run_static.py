import os
import openmc
from model_pin_cell import build_model

def run_static_simulation():
    print("[Static] Importation du modèle central...")
    
    # 1. Récupération du modèle complet depuis model_pin_cell.py
    data = build_model()
    model = data['model']
    
    # 2. Création du dossier de sortie propre
    output_dir = 'outputs/static'
    os.makedirs(output_dir, exist_ok=True)
    
    # 3. Exportation des fichiers XML et exécution d'OpenMC dans outputs/static/
    cwd = os.getcwd()
    os.chdir(output_dir)
    try:
        model.export_to_xml()
        print(f"[Static] Fichiers XML générés dans '{output_dir}'.")
        print("[Static] Lancement du calcul statique OpenMC...")
        openmc.run()
    finally:
        os.chdir(cwd)

    print(f"[Static] Calcul statique terminé avec succès ! Résultats stockés dans '{output_dir}'.")

if __name__ == '__main__':
    run_static_simulation()