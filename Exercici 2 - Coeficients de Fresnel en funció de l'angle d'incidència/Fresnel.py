# -*- coding: utf-8 -*-

"""
Exercici 2 d'avaluació continuada.
 - Mòdul: Fresnel.py
 - Revisió: 15/09/2020

 NOTES:
    - En tot el codi es treballa en unitats del SI (inclosos angles en radians)
"""

import numpy as np
from numpy import sin, asin, cos, tan, atan

"""
Variables globals
    ld: float                # longitud d'ona calculada a partir del DNI
    n: float                 # índex de refracció de ld en un vidre N-BK7 (font: http://refractiveindex.info)
"""
ld = (20 * (39418004 % 17) + 400) * 1e-9        # ld = 440 nm
n = 1.5263


def refractionAngle(phi1, n2, n1=1.):
    """
    Retorna l'angle que forma amb la normal el raig transmès en una refracció.

    Paràmetres
        phi1: float                             # angle d'incidència
        n1: float                               # índex de refracció del medi 1 (default aire: 1.)
        n2: float                               # índex de refracció del medi 2
    Retorna
        refractionAngle: float                  # angle de refracció
    """
    return asin(n1/n2*sin(phi1))


def tParallel(phi1, phi2):
    """
    Retorna coeficient de transmissió en la component paral·lela.

    Paràmetres
        phi1: float                             # angle d'incidència
        phi2: float                             # angle de refracció
    Retorna
        tParallel: float                        # coef. de transmissió paral·lela
    """
    return 2*sin(phi2)*cos(phi1) / (sin(phi1 + phi2)*cos(phi2 - phi1))


def tPerpendicular(phi1, phi2):
    """
    Retorna coeficient de transmissió en la component perpendicular.

    Paràmetres
        phi1: float                             # angle d'incidència
        phi2: float                             # angle de refracció
    Retorna
        tPerpendicular: float                   # coef. de transmissió perpendicular
    """
    return 2*sin(phi2)*cos(phi1) / sin(phi1 + phi2)


def rParallel(phi1, phi2):
    """
    Retorna coeficient de reflexió en la component paral·lela.

    Paràmetres
        phi1: float                             # angle d'incidència
        phi2: float                             # angle de refracció
    Retorna
        rParallel: float                        # coef. de reflexió paral·lela
    """
    if phi1+phi2 == np.pi/2.:                   # condició d'incidència en l'angle de Brewster
        return 0.
    else:
        return tan(phi2 - phi1) / tan(phi2 + phi1)
