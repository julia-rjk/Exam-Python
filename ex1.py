# Ecrire une fonction integ_romberg(f, a, b, epsilon=1e-6)
# permettant de calculer l’intégrale numérique de la fonction f entre les bornes a et b avec une précision epsilon selon la méthode de Romberg 
# (https://fr.wikipedia.org/wiki/M%C3%A9thode_de_Romberg).

# Il s’agit d’une méthode qui permet d’améliorer les méthodes usuelles de calcul numérique des intégrales, comme la méthode des Trapèzes.
# On montre qu’en combinant judicieusement les valeurs obtenues par la méthode des Trapèzes pour différentes 
#subdivisons de l’intervalle d’intégration, on augmente l’ordre de convergence (sa vitesse de convergence). 

import numpy as numpy

def integ_romberg(f, a, b, epsilon):
    h = b - a
    r=[]
    r.append([(h/2.0) * (f(a)+f(b))])
    for i in range(1,epsilon+1):
        h = h/2
        sum = 0
        for k in range(1,pow(2,i) ,2):
            sum = sum + f(a+ k * h)
        
        rowi = [0.5 * r[i-1][0] + sum * h]
        for j in range(1, i+1):
            rij = rowi[j-1] + (rowi[j-1]- r[i-1][j-1])/(4 ** j-1)
            rowi.append(rij)
        r .append(rowi)
    return r




def main():
    f = lambda x: x**2
    print(integ_romberg(f, 10, 20, 1))

main()