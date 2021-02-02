# -*- coding: utf-8 -*-

"""
Exercici 2 d'avaluació continuada.
 - Mòdul: Fresnel.py
 - Revisió: 15/09/2020

 NOTES:
    - En tot el codi es treballa en unitats del SI (inclosos angles en radians)
"""

import numpy as np
from numpy import sin, arcsin, cos, tan, arctan
import matplotlib.pyplot as plt

"""
Variables globals
    ld: float                # longitud d'ona calculada a partir del DNI
    n: float                 # índex de refracció de ld en un vidre N-BK7 (font: http://refractiveindex.info)
"""
ld = (20 * (39418004 % 17) + 400) * 1e-9        # ld = 440 nmç
n1 = 1.
n2 = 1.5263


def isBrewster(phi1, n2, n1=1.):
    """Comprova si s'incideix en l'angle de Brewster."""
    if phi1 == arctan(n2/n1):
        return True
    else:
        return False


def isOverCritical(phi1, n1, n2):
    """Comprova si s'incideix per sobre de l'angle crític."""
    if n2 < n1 and phi1 > arcsin(n2/n1):
        return True
    else:
        return False


def isNormal(phi1):
    """Comprova si hi ha incidència normal."""
    if phi1 == 0.:
        return True
    else:
        return False


def refractionAngle(phi1, n1, n2):
    """Retorna l'angle que forma amb la normal el raig transmès en una refracció."""
    if isNormal(phi1):                  # incidència normal
        return 0.
    elif isOverCritical(phi1, n1, n2):          # reflexió total
        pass
    else:                                       # transmissió estàndard
        return arcsin(n1/n2*sin(phi1))


def tParallel(phi1, n1, n2):
    """Retorna coeficient de transmissió en la component paral·lela."""
    phi2 = refractionAngle(phi1, n1, n2)

    if phi2 == 0:
        return 2*n1 / (n1+n2)
    elif isOverCritical(phi1, n1, n2):
        return 0.
    else:
        return 2*sin(phi2)*cos(phi1) / (sin(phi1 + phi2)*cos(phi2 - phi1))


def tPerpendicular(phi1, n1, n2):
    """Retorna coeficient de transmissió en la component perpendicular."""
    phi2 = refractionAngle(phi1, n1, n2)

    if phi2 == 0:
        return 2*n1 / (n1+n2)
    elif isOverCritical(phi1, n1, n2):
        return 0.
    else:
        return 2*sin(phi2)*cos(phi1) / sin(phi1 + phi2)


def rParallel(phi1, n1, n2):
    """Retorna coeficient de reflexió en la component paral·lela."""
    phi2 = refractionAngle(phi1, n1, n2)

    if isBrewster(phi1, n1, n2):                   # condició d'incidència en l'angle de Brewster
        return 0.
    elif isOverCritical(phi1, n1, n2):
        return -1.
    elif phi2 == 0:
        return (n1-n2) / (n1+n2)
    else:
        return tan(phi2 - phi1) / tan(phi2 + phi1)


def rPerpendicular(phi1, n1, n2):
    """Retorna coeficient de reflexió en la component paral·lela."""
    phi2 = refractionAngle(phi1, n1, n2)

    if isOverCritical(phi1, n1, n2):
        return 1.
    elif phi2 == 0:
        return (n1-n2) / (n1+n2)
    else:
        return sin(phi2 - phi1) / sin(phi2 + phi1)


phi = np.linspace(0, np.pi/2., num=500)
fresnel_coeffs = np.zeros((4, 500))
for n, phi1 in enumerate(phi):
    phi2 = refractionAngle(phi1, 1.000, n)
    fresnel_coeffs[0, n] = tParallel(phi1, n1, n2)
    fresnel_coeffs[1, n] = tPerpendicular(phi1, n1, n2)
    fresnel_coeffs[2, n] = rParallel(phi1, n1, n2)
    fresnel_coeffs[3, n] = rPerpendicular(phi1, n1, n2)


plt.plot(phi, fresnel_coeffs[0, :], color="red")
plt.plot(phi, fresnel_coeffs[1, :], color="green")
plt.plot(phi, fresnel_coeffs[2, :], color="blue")
plt.plot(phi, fresnel_coeffs[3, :], color="orange")

plt.plot(phi, fresnel_coeffs[0, :]-fresnel_coeffs[1, :], color="cyan")
plt.plot(phi, fresnel_coeffs[2, :]-fresnel_coeffs[3, :], color="yellow")

plt.show()
