# -*- coding: utf-8 -*-

"""
Exercici 2 d'avaluació continuada.
 - Mòdul: Fresnel.py
 - Revisió: 02/02/2021

 NOTES:
    - En tot el codi es treballa en unitats del SI (inclosos angles en radians)
"""

import numpy as np
from numpy import sin, arcsin, cos, tan, arctan
import matplotlib.pyplot as plt

"""
Variables globals
    n1: float                # índex de refracció del medi 1 (aire 1.000)
    n2: float                # índex de refracció del medi 2 (vidre N-BK7 per lambda = 440 nm, font: http://refractiveindex.info)
"""
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


def refractionAngle(phi1, n1, n2):
    """Retorna l'angle que forma amb la normal el raig transmès en una refracció."""
    if phi1 == 0.:                              # incidència normal
        return 0.
    elif isOverCritical(phi1, n1, n2):          # reflexió total
        pass
    else:                                       # transmissió estàndard
        return arcsin(n1/n2*sin(phi1))


def tParallel(phi1, n1, n2):
    """Retorna coeficient de transmissió en la component paral·lela."""
    phi2 = refractionAngle(phi1, n1, n2)

    if phi2 == 0:                               # incidència normal
        return 2*n1 / (n1+n2)
    elif isOverCritical(phi1, n1, n2):          # reflexió total
        return 0.
    else:                                       # transmissió estàndard
        return 2*sin(phi2)*cos(phi1) / (sin(phi1 + phi2)*cos(phi2 - phi1))


def tPerpendicular(phi1, n1, n2):
    """Retorna coeficient de transmissió en la component perpendicular."""
    phi2 = refractionAngle(phi1, n1, n2)

    if phi2 == 0:                               # incidència normal
        return 2*n1 / (n1+n2)
    elif isOverCritical(phi1, n1, n2):          # reflexió total
        return 0.
    else:                                       # transmissió estàndard
        return 2*sin(phi2)*cos(phi1) / sin(phi1 + phi2)


def rParallel(phi1, n1, n2):
    """Retorna coeficient de reflexió en la component paral·lela."""
    phi2 = refractionAngle(phi1, n1, n2)

    if isBrewster(phi1, n1, n2):                 # incidència en l'angle de Brewster
        return 0.
    elif isOverCritical(phi1, n1, n2):           # reflexió total
        return -1.
    elif phi2 == 0:                              # incidència normal
        return (n1-n2) / (n1+n2)
    else:                                        # reflexió estàndard
        return tan(phi2 - phi1) / tan(phi2 + phi1)


def rPerpendicular(phi1, n1, n2):
    """Retorna coeficient de reflexió en la component paral·lela."""
    phi2 = refractionAngle(phi1, n1, n2)

    if isOverCritical(phi1, n1, n2):             # reflexió total
        return 1.
    elif phi2 == 0:                              # incidència normal
        return (n1-n2) / (n1+n2)
    else:                                        # reflexió estàndard
        return sin(phi2 - phi1) / sin(phi2 + phi1)


N = 500                                          # nombre d'iteracions (resoució en phi1)
phi = np.linspace(0, np.pi/2., num=N)            # valors de phi1
t_parallel = np.zeros(N)                         # coef. de transmissió paral·lela
t_perpendicular = np.zeros(N)                    # coef. de transmissió perpendicular
r_parallel = np.zeros(N)                         # coef. de reflexió paral·lela
r_perpendicular = np.zeros(N)                    # coef. de reflexió perpendicular

"""Calcula els coeficients de Fresnel per als angles d'incidència de l'array phi"""
for i, phi1 in enumerate(phi):
    t_parallel[i] = tParallel(phi1, n1, n2)
    t_perpendicular[i] = tPerpendicular(phi1, n1, n2)
    r_parallel[i] = rParallel(phi1, n1, n2)
    r_perpendicular[i] = rPerpendicular(phi1, n1, n2)

"""
Calcula diferència entre els coeficients de transmissió paral·lels i perpendiculars
íd. per als de reflexió
"""
t_diff = t_parallel - t_perpendicular
r_diff = r_parallel - r_perpendicular


"""
Aproxima els coeficients de transmissió i reflexió per la mitjana en les components
i estudia en quin rang l'error és inferior al 5%
"""
t_approx = .5 * (t_parallel+t_perpendicular)
r_approx = .5 * (r_parallel+r_perpendicular)




plt.plot(phi, t_parallel, color="red")
plt.plot(phi, t_perpendicular, color="green")
plt.plot(phi, t_approx, color="blue")
"""
plt.plot(phi, r_parallel, color="blue")
plt.plot(phi, r_perpendicular, color="orange")
"""
"""
plt.plot(phi, t_diff, color="cyan")
plt.plot(phi, r_diff, color="yellow")
"""
plt.show()
