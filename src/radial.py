"""Séance 3 : températures radiales (paroi, gaine, pastille) le long du crayon moyen."""

import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

from parametres import P_PRIMAIRE, D_GAINE, E_GAINE, K_GAINE, D_PASTILLE_UO2, H_JEU
from geometrie import A, D_h, m_sc
from profil_axial import z, q, h, T

# --- Étape 1 : convection eau -> paroi de la gaine ---
# Propriétés de l'eau en chaque z (à partir de l'enthalpie h(z))
mu = np.array([PropsSI('V', 'H', hi, 'P', P_PRIMAIRE, 'Water') for hi in h])
k_eau = np.array([PropsSI('L', 'H', hi, 'P', P_PRIMAIRE, 'Water') for hi in h])
Pr = np.array([PropsSI('PRANDTL', 'H', hi, 'P', P_PRIMAIRE, 'Water') for hi in h])

# ρv = ṁ_sc / A est constant le long du canal : seule la viscosité varie
Re = m_sc * D_h / (A * mu)
Nu = 0.023 * Re**0.8 * Pr**0.4                # Dittus-Boelter (fluide chauffé)
h_conv = Nu * k_eau / D_h                     # [W/m²/K]

# Loi de Newton : q' traverse la surface extérieure de la gaine (π D_gaine par mètre)
dT_conv = q / (h_conv * np.pi * D_GAINE)      # [K]
T_paroi = T + dT_conv                         # [K]

# --- Étape 2 : conduction à travers la gaine (cylindre creux) ---
# En régime permanent, sans source de chaleur dans la gaine, la même puissance q'
# traverse chaque cylindre de rayon r : q' = -2π r k dT/dr.
# On intègre de r_int à r_ext  ->  ΔT = q' / (2π k) · ln(r_ext / r_int)
r_ext_gaine = D_GAINE / 2
r_int_gaine = r_ext_gaine - E_GAINE
dT_gaine = q / (2 * np.pi * K_GAINE) * np.log(r_ext_gaine / r_int_gaine)   # [K]
T_gaine_int = T_paroi + dT_gaine              # face intérieure de la gaine [K]

# --- Étape 3 : traversée du jeu pastille-gaine (rempli d'hélium) ---
# Le jeu est très mince (quelques dizaines de µm) : on ne le décrit pas couche par couche,
# on le résume par une conductance globale h_jeu [W/m²/K] (conduction dans l'hélium
# + rayonnement + contacts pastille-gaine), comme une « convection » à travers un film.
# La chaleur traverse la surface de la pastille (π D_pastille par mètre) :
#   ΔT_jeu = q' / (h_jeu π D_pastille)
r_pastille = D_PASTILLE_UO2 / 2
e_jeu = r_int_gaine - r_pastille              # largeur du jeu à froid [m]
dT_jeu = q / (H_JEU * np.pi * D_PASTILLE_UO2)  # [K]
T_pastille_surf = T_gaine_int + dT_jeu        # surface de la pastille d'UO2 [K]

if __name__ == "__main__":
    i_max = np.argmax(T_paroi)
    print(f"h_conv min/max = {h_conv.min()/1e3:.1f} / {h_conv.max()/1e3:.1f} kW/m²/K")
    print(f"ΔT_conv max    = {dT_conv.max():.1f} K")
    print(f"T_paroi max    = {T_paroi.max() - 273.15:.1f} °C à z = {z[i_max]:.2f} m")

    j_max = np.argmax(T_gaine_int)
    print(f"ΔT_gaine max   = {dT_gaine.max():.1f} K")
    print(f"T_gaine_int max = {T_gaine_int.max() - 273.15:.1f} °C à z = {z[j_max]:.2f} m")

    k_max = np.argmax(T_pastille_surf)
    print(f"e_jeu (à froid) = {e_jeu*1e6:.0f} µm")
    print(f"ΔT_jeu max     = {dT_jeu.max():.1f} K")
    print(f"T_pastille_surf max = {T_pastille_surf.max() - 273.15:.1f} °C à z = {z[k_max]:.2f} m")
