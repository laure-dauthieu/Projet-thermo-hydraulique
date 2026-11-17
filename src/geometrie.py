import math
from CoolProp.CoolProp import PropsSI
from parametres import *

# 1. Section de passage du sous-canal [m²], 
# (surface où l’eau peut circuler entre 4 crayons : le carré du réseau moins le disque du crayon)
A = PAS**2 - math.pi * D_GAINE**2 / 4

# 2. Diamètre hydraulique [m]  ->  D_h = 4A / (π D_gaine)
D_h = 4 * A / (math.pi * D_GAINE)

# 3. Puissance linéique moyenne [W/m], 
# (On divise la puissance totale par la longueur totale de tous les crayons mis bout à bout)
N_CRAYONS = N_ASSEMBLAGES * N_CRAYONS_PAR_ASSEMBLAGE
q_lin = P_TH / (N_CRAYONS * H_ACTIVE)

# 4. Débit du cœur par bilan d'énergie [kg/s], 
# (puissance = débit × variation d’enthalpie, la chaleur du coeur sert à faire monter l'eau de 286°C à 323°C)
h_entree = PropsSI('H', 'T', T_ENTREE, 'P', P_PRIMAIRE, 'Water')   # [J/kg]
h_sortie = PropsSI('H', 'T', T_SORTIE, 'P', P_PRIMAIRE, 'Water')    # [J/kg]
m_coeur = P_TH / (h_sortie - h_entree)

# 5. Débit, vitesse et Reynolds dans un sous-canal
m_sc = m_coeur / N_CRAYONS # m_coeur / N_CRAYONS
T_moy = (T_ENTREE + T_SORTIE) / 2
rho = PropsSI('D', 'T', T_moy, 'P', P_PRIMAIRE, 'Water')   # [kg/m³]
mu = PropsSI('V', 'T', T_moy, 'P', P_PRIMAIRE, 'Water')   # [Pa.s]
v = m_sc / (rho * A)
Re = rho * v * D_h / mu

# Affichage (seulement si on lance ce fichier directement, pas lors d'un import)
if __name__ == "__main__":
    print(f"A      = {A:.3e} m²")
    print(f"D_h    = {D_h*1e3:.2f} mm")
    print(f"q_lin = {q_lin:.2f} W/m")
    print(f"m_coeur = {m_coeur:.2f} kg/s")
    print(f"m_sc   = {m_sc:.3f} kg/s")
    print(f"v      = {v:.2f} m/s")
    print(f"Re     = {Re:.2e}")