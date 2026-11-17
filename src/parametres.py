"""Paramètres de référence d'un REP 900 MWe (ordres de grandeur, unités SI)."""

# --- Cœur ---
P_TH = 2785e6                 # puissance thermique du cœur [W]
P_PRIMAIRE = 155e5            # pression primaire [Pa]
T_ENTREE = 286 + 273.15       # température d'entrée cœur [K]
T_SORTIE = 323 + 273.15       # température de sortie cœur [K]

N_ASSEMBLAGES = 157           # nombre d'assemblages [-]
N_CRAYONS_PAR_ASSEMBLAGE = 264  # crayons combustible par assemblage 17x17 [-]
H_ACTIVE = 3.66               # hauteur active du combustible [m]

# --- Géométrie du crayon et du réseau ---
D_GAINE = 9.5e-3              # diamètre extérieur de la gaine [m]
E_GAINE = 0.57e-3             # épaisseur de la gaine [m]
D_PASTILLE_UO2 = 8.19e-3      # diamètre de la pastille UO2 [m]
PAS = 12.6e-3                 # pas du réseau carré [m]

# --- Propriétés des matériaux ---
K_GAINE = 16.0                # conductivité thermique du Zircaloy [W/m/K] (ordre de grandeur vers 350 °C)
