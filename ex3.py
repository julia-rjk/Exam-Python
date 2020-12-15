import numpy as numpy
import math


class Ville:
    # initialisation d’une ville sans destination
    def __init__(self):
        self.destinations = numpy.array([]).reshape(-1, 2)
        #  on se place dans le cas où les coordonnées des destinations sont entières, comprises entre 0 (inclus) et TAILLE = 50 (exclus)
        self.TAILLE = 50

    # création de n destinations aléatoires
    def aleatoire(self, n):
        self.destinations = numpy.random.randint(self.TAILLE, size=(n, 2))

    # retourne le nombre total (entier) de trajets :(𝑛−1)!/2 (utiliser math.factorial() ).
    def nb_trajet(self):
        length = len(self.destinations)
        print(length)
        if length > 2:
            return int(math.factorial(length - 1) / 2)
        if length > 0:
            return 1
        return 0

    # retourne la distance (Manhattan) entre les deux destinations de numéro i et j
    # 𝑑(𝐴,𝐵) =|𝑥𝐵−𝑥𝐴|+|𝑦𝐵−𝑦𝐴|
    def distance(self, i, j):
        return numpy.abs(self.destinations[i] - self.destinations[j]).sum()

    # Plus proche voisin:
    # retourne la destination la plus proche de la destination (au sens de Ville.distance()), hors les destinations de la liste exclus
    def plus_proche(self, i, exclus=[]):
        length = len(self.destinations)
        villeAutour = []

        for j in range(length):
            dest = self.destinations[j]
            if(numpy.any(dest != i) and dest not in exclus):
                villeAutour.append(dest)

        distances = []
        for k in range(len(villeAutour)):
            dist = villeAutour[k]
            distances.append(self.distance(i, dist))

        return villeAutour[numpy.argmin(distances)]

    def optimisation_trajet(self, trajet):
        """
        Retourne le trajet le plus court de tous les trajets « voisins » à
        `trajet` (i.e. résultant d'une simple interversion de 2 étapes).
        """

        size = len(self.destinations)
        trajets = [trajet.interversion(i, j)
                   for i in range(size) for j in range(i+1, size)]
        longueurs = [t.longueur() for t in trajets]
        opt = trajets[numpy.argmin(longueurs)]
        if opt.longueur() > trajet.longueur():
            opt = trajet

        return opt

    def getDestinations(self):
        return self.destinations


class Trajet:

    # initialisation sur une ville. Si la liste etapes n’est pas spécifiée, le trajet par défaut est celui suivant les destinations de ville.
    def __init__(self, ville, etapes=None):
        self.ville = ville

        # self.etapes = ville.destinations
        if etapes is None:
            self.etapes = numpy.arange(len(self.ville.destinations))
        else:
            self.etapes = numpy.array(etapes)

    # retourne la longueur totale du trajet bouclé(i.e. revenant à son point de départ).
    def longueur(self):
        nbEtapes = len(self.etapes) - 1
        for i in range(nbEtapes):
            longueur = sum(self.ville.distance(
                self.etapes[i], self.etapes[i+1]))

        longueur += self.ville.distance(self.etapes[-1], self.etapes[0])

        return longueur

    def getVille(self):
        return self.ville


def main():
    ville = Ville()
    ville.aleatoire(50)
    ville.nb_trajet()
    dest = ville.getDestinations()
    print(dest[0])
    print("Ville la plus proche : ", ville.plus_proche(dest[0]))
    # ahi = ville.trajet_voisins()
    # print(" ahi : " , ahi.getVille().getDestinations())
    return 0


main()
