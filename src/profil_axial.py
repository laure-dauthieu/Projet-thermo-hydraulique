"""Séance 2 : profil axial de puissance linéique et de température de l'eau
dans le sous-canal moyen (flux en cosinus, z = 0 au milieu du cœur)."""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

from parametres import H_ACTIVE, P_PRIMAIRE
from geometrie import m_sc, q_lin, h_entree

H = H_ACTIVE

# 1. Axe vertical : z va de -H/2 (bas du cœur) à +H/2 (haut du cœur) [m]
z = np.linspace(-H / 2, H / 2, 200)

# 2. Puissance linéique en cosinus [W/m]
# La moyenne d'un cosinus sur [-H/2, H/2] vaut 2/π de son maximum,
# donc pour retrouver la moyenne q_lin il faut q'_max = q_lin × π/2.
q_max = q_lin * np.pi / 2
q = q_max * np.cos(np.pi * z / H)

# 3. Enthalpie de l'eau h(z) par bilan d'énergie [J/kg]
# Entre le bas du cœur et la cote z, l'eau a reçu ∫ q'(z') dz'.
# Divisé par le débit m_sc, ça donne l'augmentation d'enthalpie :
#   h(z) = h_entree + q'_max H / (π m_sc) × [sin(π z / H) + 1]
h = h_entree + q_max * H / (np.pi * m_sc) * (np.sin(np.pi * z / H) + 1)

# 4. Température de l'eau T(z) : CoolProp convertit (h, P) -> T [K]
T = np.array([PropsSI('T', 'H', hi, 'P', P_PRIMAIRE, 'Water') for hi in h])

if __name__ == "__main__":
    print(f"q'_max          = {q_max/1e3:.1f} kW/m")
    print(f"T entrée (bas)  = {T[0] - 273.15:.1f} °C")
    print(f"T sortie (haut) = {T[-1] - 273.15:.1f} °C")

    # 5. Graphiques : z en ordonnée (comme dans le cœur, bas en bas)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6), sharey=True)

    ax1.plot(q / 1e3, z)
    ax1.set_xlabel("Puissance linéique q'(z) [kW/m]")
    ax1.set_ylabel("Cote z [m] (0 = milieu du cœur)")
    ax1.set_title("Profil de puissance (cosinus)")
    ax1.grid(True)

    ax2.plot(T - 273.15, z, color="tab:red")
    ax2.set_xlabel("Température de l'eau T(z) [°C]")
    ax2.set_title("Échauffement de l'eau")
    ax2.grid(True)

    fig.suptitle("Sous-canal moyen — REP 900 MWe")
    fig.tight_layout()

    # Dossier figures/ à la racine du projet (créé s'il n'existe pas)
    dossier_figures = Path(__file__).resolve().parent.parent / "figures"
    dossier_figures.mkdir(exist_ok=True)
    fig.savefig(dossier_figures / "profil_axial.png", dpi=150)

    plt.show()
