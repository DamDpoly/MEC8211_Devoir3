"""
Labo3-equation_de_la_chaleur.py : Vérification de code par la Méthode des Solutions Manufacturées (MMS)
Créé en février 2025
Auteur: David Vidal
     
Modifié par: Damien Defrance
Pour le devoir 2 de MEC8211

"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Définition des variables symboliques
r, t, R, l, C_e, D_eff, k = sp.symbols('r t R l C_e D_eff k')

# solution MMS
C_MMS = (((r**2)/(R**2)) * (C_e + sp.exp(-l * t)) - sp.exp(-l * t))

# Appliquer l'opérateur sur la solution MMS
C_t = sp.diff(C_MMS, t)
C_r = sp.diff(C_MMS, r)
C_rr = (1/r) * sp.diff(C_r * r, r)

# Terme source
source = C_t - C_rr * D_eff + k * C_MMS

# Conditions aux limites et initiales
C_initial = C_MMS.subs(t, 0)
C_boundary_R = C_MMS.subs(r, R)
dCdr_boundary_r0 = C_r.subs(r, 0)

# Afficher les résultats
print("Dérivée en temps :")
print(C_t)
print("\nDérivée première :")
print(C_r)
print("\nLaplacien :")
print(C_rr)
print("\nTerme source :")
print(source)
print("\nCondition initiale C(r, 0) :")
print(C_initial)
print("\nCondition frontière C(R, t) :")
print(C_boundary_R)
print("\nCondition frontière Neumann dC/dx(t,0) :")
print(dCdr_boundary_r0)

# Créer une fonction appelable pour l'expression symbolique
f_C_MMS_symbolic = sp.lambdify([r, t, R, l, C_e, D_eff, k], C_MMS, "numpy")
f_source_symbolic = sp.lambdify([r, t, R, l, C_e, D_eff, k], source, "numpy")

# taille du domaine
tmin, tmax = 0, 1
rmin, rmax = 0, 1
nt, nr = 50, 50

# Établir une grille régulière de points d'interpolation
tdom = np.linspace(tmin, tmax, nt)
rdom = np.linspace(rmin, rmax, nr)
ti, ri = np.meshgrid(tdom, rdom, indexing='ij')

# Évaluer la fonction MMS et le terme source sur le maillage
z_MMS = f_C_MMS_symbolic(ri, ti, 1, 1, 1, 1, 1)
z_source = f_source_symbolic(ri, ti, 1, 1, 1, 1, 1)

# Tracer les résultats
plt.figure()
plt.contourf(ri, ti, z_MMS, levels=50)
plt.colorbar()
plt.title('Solution Manufacturée')
plt.xlabel('r')
plt.ylabel('t')
plt.show()

plt.figure()
plt.contourf(ri, ti, z_source, levels=50)
plt.colorbar()
plt.title('Terme Source')
plt.xlabel('r')
plt.ylabel('t')
plt.show()

