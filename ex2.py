# Ex 2: Équation d’état de l’eau à partir de la dynamique moléculaire

# Afin de modéliser les planètes de type Jupiter, Saturne, ou même des exo-planètes très massives (dites « super-Jupiters »), la connaissance de l’équation d’état des composants est nécessaire. Ces équations d’état doivent être valables jusqu’à plusieurs centaines de méga-bar ; autrement dit, celles-ci ne sont en aucun cas accessibles expérimentalement. On peut cependant obtenir une équation d’état numériquement à partir d’une dynamique moléculaire.

# Le principe est le suivant : on place dans une boite un certain nombre de particules régies par les équations microscopiques
# (Newton par exemple, ou même par des équations prenant en considération la mécanique quantique) puis on laisse celles-ci évoluer dans la boite ;
#  on calcule à chaque pas de temps l’énergie interne à partir des intéractions électrostatiques et la pression à partir du tenseur des contraintes.
# On obtient en sortie l’évolution du système pour une densité fixée (par le choix de taille de la boite) et une température fixée (par un algorithme de thermostat que nous ne détaillerons pas ici).

# On se propose d’analyser quelques fichiers de sortie de tels calculs pour l’équation d’état de l’eau à très haute pression.
# Les fichiers de sortie sont disponibles ici; leur nom indique les conditions thermodynamiques correspondant au fichier, p.ex. 6000K_30gcc.out pour T=6000
# K et ρ=30 gcc.
# Le but est, pour chaque condition température-densité, d’extraire l’évolution de l’énergie et de la pression au cours du temps, puis d’en extraire la valeur moyenne ainsi que les fluctuations. Il arrive souvent que l’état initial choisi pour le système ne corresponde pas à son état d’équilibre, et qu’il faille donc « jeter » les quelques pas de temps en début de simulation qui correspondent à cette relaxation du système. Pour savoir combien de temps prend cette relaxation, il sera utile de tracer l’évolution au cours du temps de la pression et l’énergie pour quelques simulations. 
# Une fois l’équation d’état P(ρ,T) et E(ρ,T) extraite, on pourra tracer le réseau d’isothermes.

import numpy as N
import matplotlib.pyplot as P

tolerance = 1e-8 

class Simulation:

    def __init__(self, temp, dens, path):

        self.temp = float(temp)
        self.dens = float(dens)
        tmp = N.loadtxt(path, skiprows=1).T
        self.pot = tmp[0]
        self.kin = tmp[1]
        self.tot = self.pot + self.kin
        self.press = tmp[2]

    def __str__(self):

        return "Simulation at {:.0f} g/cc and {:.0f} K ; {:d} timesteps". \
            format(self.dens, self.temp, len(self.pot))

    def thermo(self, skipSteps=0):
        
        return {'T': self.temp,
                'rho': self.dens,
                'E': self.tot[skipSteps:].mean(),
                'P': self.press[skipSteps:].mean(),
                'dE': self.tot[skipSteps:].std(),
                'dP': self.press[skipSteps:].std()}

    def plot(self, skipSteps=0):
       
        fig, (axen, axpress) = P.subplots(2, sharex=True)
        axen.plot(list(range(skipSteps, len(self.tot))), self.tot[skipSteps:],
                  'rd--')
        axen.set_title("Internal energy (Ha)")
        axpress.plot(list(range(skipSteps, len(self.press))), self.press[skipSteps:],
                     'rd--')
        axpress.set_title("Pressure (GPa)")
        axpress.set_xlabel("Timesteps")

        P.show()


def mimic_simulation(filename):
    with open(filename, 'w') as f:
        f.write("""Potential energy (Ha)	Kinetic Energy (Ha)	Pressure (GPa)
-668.2463567264        	0.7755612311   		9287.7370229824
-668.2118514558        	0.7755612311		9286.1395903265
-668.3119088218        	0.7755612311		9247.6604398856
-668.4762735176        	0.7755612311		9191.8574820856
-668.4762735176        	0.7755612311		9191.8574820856
""")


def test_Simulation_init():
    mimic_simulation("equationEtat_simuTest.out")
    s = Simulation(10, 10, "equationEtat_simuTest.out")
    assert len(s.kin) == 5
    assert abs(s.kin[2] - 0.7755612311) < tolerance
    assert abs(s.pot[1] + 668.2118514558) < tolerance


def test_Simulation_str():
    mimic_simulation("equationEtat_simuTest.out")
    s = Simulation(10, 20, "equationEtat_simuTest.out")
    assert str(s) == "Simulation at 20 g/cc and 10 K ; 5 timesteps"


def test_Simulation_thermo():
    mimic_simulation("equationEtat_simuTest.out")
    s = Simulation(10, 20, "equationEtat_simuTest.out")
    assert abs(s.thermo()['T'] - 10) < tolerance
    assert abs(s.thermo()['rho'] - 20) < tolerance
    assert abs(s.thermo()['E'] + 667.56897157674) < tolerance
    assert abs(s.thermo()['P'] - 9241.0504034731) < tolerance
    assert abs(s.thermo(3)['E'] + 667.7007122865) < tolerance
    assert abs(s.thermo(3)['P'] - 9191.8574820856) < tolerance



if __name__ == '__main__':
    a0 = 0.52918      
    amu = 1.6605      
    k_B = 3.16681e-6
    nk_GPa = a0 ** 3 * k_B * 2.942e4 / 6 / amu
    nsteps = 200
    temps = [6000, 20000, 50000]
    colors = {6000: 'r', 20000: 'b', 50000: 'k'}
    denss = [7, 15, 25, 30]
    keys = ['T', 'rho', 'E', 'dE', 'P', 'dP']
    eos = dict.fromkeys(keys, N.zeros(0))

    for t, rho in [(t, rho) for t in temps for rho in denss]:
        filenm = "outputs/{}K_{:0>2d}gcc.out".format(t, rho)
        s = Simulation(t, rho, filenm)
        for key in keys:
            eos[key] = N.append(eos[key], s.thermo(nsteps)[key])

    fig, (axen, axpress) = P.subplots(2, sharex=True)
    fig.suptitle("High-pressure equation of state for water", size='x-large')
    axen.set_title("Energy")
    axen.set_ylabel("U / NkT")
    axpress.set_title("Pressure")
    axpress.set_ylabel("P / nkT")
    axpress.set_xlabel("rho (g/cc)")
    for t in temps:
        sel = eos['T'] == t
        axen.errorbar(x=eos['rho'][sel], y=eos['E'][sel] / k_B / t,
                      yerr=eos['dE'][sel] / k_B / t, fmt=colors[t] + '-')
        axpress.errorbar(x=eos['rho'][sel],
                         y=eos['P'][sel] / eos['rho'][sel] / nk_GPa / t,
                         yerr=eos['dP'][sel] / eos['rho'][sel] / nk_GPa / t,
                         fmt=colors[t] + '-',
                         label="{} K".format(t))
    axpress.legend(loc='best')
    P.show()