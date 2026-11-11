# OpenMC PWR Pin-Cell Simulation & Neutronics Analysis

Projet complet de modélisation neutronique d'un crayon de combustible REP avec OpenMC.

## Fonctionnalités & Structure du projet
- `model_pin_cell.py` : Définition de la géométrie 3D du crayon, des matériaux ($UO_2$, Zircaloy-4, $H_2O$) et des tallies.
- `run_static.py` & `run_depletion.py` : Exécution des simulations Monte Carlo en criticité statique et déplétion isotopique.
- `plot_results.py` : Extraction des spectres neutroniques, de l'évolution du $k_{\text{eff}}$ et de l'inventaire lourd ($^{235}\text{U}$, $^{239}\text{Pu}$).
- `export_vtk.py` : Génération des cartes 3D de matériaux au format VTK pour visualisation dans ParaView.

## Résultats
Les graphiques et exportations générés sont automatiquement sauvegardés dans `/figures` et `/vtk_exports`.

## Préréquis
- Python 3.10+
- OpenMC avec données ENDF/B-VIII.1
- NumPy, Matplotlib, h5py