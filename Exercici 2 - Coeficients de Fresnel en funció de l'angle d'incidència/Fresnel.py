"""
Exercici 2 d'avaluació continuada.

 - Mòdul: Fresnel.py
 - Revisió: 03/02/2021
"""


import numpy as np
from numpy import sin, arcsin, cos, tan, arctan
import matplotlib.pyplot as plt


"""
Funcions
"""


def csvwrite(nom_out, *cols):
    """Escriu el contingut d'un conjunt d'arrays en un fitxer .csv per columnes."""
    lens = set([len(col) for col in cols])
    if len(lens) != 1:
        raise IndexError("Els arrays no tenen la mateixa mida.")

    else:
        with open(nom_out, "w") as fout:
            for i in range(min(lens)):
                for j in range(len(cols) - 1):
                    fout.write("%e," % (cols[j][i]))
                fout.write("%e\n" % (cols[j+1][i]))


def isBrewster(phi1, n1, n2):
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


"""
Execució del programa
"""

n1 = 1.                                          # índex de refracció del medi 1 (aire 1.000)
n2 = 1.5263                                      # índex de refracció del medi 2 (vidre N-BK7 per lambda = 440 nm, font: http://refractiveindex.info)
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
t_approx = .5 * (t_parallel+t_perpendicular)                                    # valor aproximat de t
r_approx = .5 * (r_parallel+r_perpendicular)                                    # valor aproximat de r

t_parallel_errs = np.abs((t_approx-t_parallel) / t_parallel)                    # error relatiu respecte t paral·lel
t_perpendicular_errs = np.abs((t_approx-t_perpendicular) / t_perpendicular)     # error relatiu respecte t perpendicular
r_parallel_errs = np.abs((r_approx-r_parallel) / r_parallel)                    # error relatiu respecte r paral·lel
r_perpendicular_errs = np.abs((r_approx-r_perpendicular) / r_perpendicular)     # error relatiu respecte r perpendicular

t_phi1_intervals = []                                                           # intervals d'error acotat per t
r_phi1_intervals = []                                                           # intervals d'error acotat per r

err_bound = .05                                                                 # error relatiu màxim
current_phi_range = []
for (phi1, t_parallel_err, t_perpendicular_err) in zip(phi, t_parallel_errs, t_perpendicular_errs):
    if t_parallel_err < err_bound and t_perpendicular_err < err_bound:
        current_phi_range.append(phi1)
    else:
        if len(current_phi_range) > 0:
            t_phi1_intervals.append([np.min(current_phi_range), np.max(current_phi_range)])
            current_phi_range.clear()
        else:
            continue

for (phi1, r_parallel_err, r_perpendicular_err) in zip(phi, r_parallel_errs, r_perpendicular_errs):
    if r_parallel_err < err_bound and r_perpendicular_err < err_bound:
        current_phi_range.append(phi1)
    else:
        if len(current_phi_range) > 0:
            r_phi1_intervals.append([np.min(current_phi_range), np.max(current_phi_range)])
            current_phi_range.clear()
        else:
            continue

print("Rang(s) angular(s) en què t_approx té un error menor al 5%:", end=" ")
for interval in t_phi1_intervals:
    print("[%.2f, %.2f]" % (interval[0], interval[1]), end=" ")
print()

print("Rang(s) angular(s) en què r_approx té un error menor al 5%:", end=" ")
for interval in r_phi1_intervals:
    print("[%.2f, %.2f]" % (interval[0], interval[1]), end=" ")
print()


"""
Escriptura de les dades en un fitxer .csv
"""
csvwrite("./Fresnel.csv", phi, t_parallel, t_perpendicular, r_parallel, r_perpendicular, t_diff, r_diff)


"""
Representació gràfica dels resultats
"""
plt.plot(np.degrees(phi), t_parallel, color="red")
plt.plot(np.degrees(phi), t_perpendicular, color="green")
plt.plot(np.degrees(phi), r_parallel, color="blue")
plt.plot(np.degrees(phi), r_perpendicular, color="orange")

plt.show()
