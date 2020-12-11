# Groupe 11.58 : Simulation informatique

import math
import matplotlib.pyplot as plt
import numpy as np


# Paramètres
L_barge = 0.43       # Longueur de la barge [m]  /!\ GEOMETRIE DE NOTRE BARQUE DIFFICILE A MODELISER,
lar_barge = 0.43     # Largeur de la barge [m]       CAR ELLE EST COMPOSEE DE BOUTEILLES (FLOTTEURS)
hc = 0.05            # Hauteur immergée [m]
dens_eau = 1000      # Densité de l'eau [kg/m³]
g = 9.81             # Constante de pesanteur [m/s²]
hcg = 0.146          # Hauteur du centre de gravité [m]
d = 0.53             # Distance à laquelle objet est déplacé [m]
m_charge = 0.275     # Masse de l'objet déplacé [kg]
D = 0.5              # Coefficient d'amortissement [kg/s]
I = 3.06             # Moment d'inertie [kg.m²]
h_barge = 0.15       # Hauteur de la barge [m]
m_tot = 3.8          # Masse totale du système [kg]
h_cable = 0          # Distance à laquelle grappin est descendu [m]


# Fonction qui calcule l'angle d'inclinaison selon la formule trouvée dans notre modélisation physique
l = 0.35
def calcul_angle_instant():
    return -math.degrees(math.atan((6*hc*d*m_charge)/(l**2*m_tot)))


###########################
### PREMIERE SIMULATION ###
###########################

# Initialisation et création tableaux
dt = 0.001     # pas (dt) [s]
end = 60.0     # durée [s]
x_0 = 0.0      # position initiale [m]
v_0 = 0.0      # vitesse initiale [m/s]

t = np.arange(0, end, dt)
C = np.empty_like(t)          
v_ang = np.empty_like(t)
a = np.empty_like(t)
theta = np.empty_like(t)


# Acceleration, vitesse, position [angulaires] et couples
def simulation():
    """
    pre: /
    post: Exécute une simulation jusqu'à t=end par pas de dt=step.
          Remplit les listes x, v, a des positions, vitesses et accélérations.
    """
    # Conditions initiales
    v_ang[0] = 0
    theta[0] = 0
    C[0] = -d*m_charge*g      # couple engendré par la charge (le reste est à l'équilibre en t=0)
    a[0] = (-d*m_charge*g) /I
    
    for i in range(len(t)-1):
        # Position angulaire
        theta[i+1] = theta[i] + v_ang[i]*dt
        # Calcul de la somme des couples
        C[i+1] = (L_barge*lar_barge*hc*dens_eau*g)*abs((math.sin(math.radians(theta[i]))*hcg)- \
                 (L_barge**2)*math.sin(math.radians(theta[i]))/(12*hc)) - (d*m_charge*g)
        # Accélération
        a[i+1] = (C[i]+ (-D*v_ang[i]))/I
        # Vitesse
        v_ang[i+1] = v_ang[i] + a[i]*dt

       
# Angle maximal
def angle_maximal():
    theta_submersion = math.atan((h_barge-hc)/(L_barge/2))
    theta_soulevement = math.atan(hc/(L_barge/2))
    if theta_submersion > theta_soulevement:
        return math.degrees(theta_soulevement)
    else:
        return math.degrees(theta_submersion)

a_max = angle_maximal()
a_max2 = -(angle_maximal())
x = [0, end]
y = [a_max, a_max]
x2 = [0,end]
y2 = [a_max2, a_max2]
    
#Calcul Hc   
def calcul_hc():
    hc = m_tot / (dens_eau*L_barge*lar_barge)
    return hc



### PREMIERE FENETRE DE GRAPHIQUES : accélération, vitesse, position [angulaires] et couple 
def graphiques():
    
    plt.figure("Graphiques 1")
    plt.subplot(4,1,1)
    plt.title("Accélération, vitesse, position [angulaires] et couples")
    plt.plot(t,C, color = 'darkturquoise', label="Couples", linewidth=2)
    plt.ylabel('[N*m]')
    plt.legend()
    plt.subplot(4,1,2)
    plt.plot(t,v_ang, color = 'steelblue', label="Vitesse angulaire", linewidth=2)
    plt.ylabel('[rad/s]')
    plt.legend()
    plt.subplot(4,1,3)
    plt.plot(t,a, color ='rebeccapurple', label="Accéleration angulaire", linewidth=2)
    plt.ylabel('[rad/s²]')
    plt.legend()
    plt.subplot(4,1,4)
    plt.plot(t, theta, color = 'crimson', label="Angle d'inclinaison", linewidth=2)
    plt.ylabel('[degrés]')
    plt.plot(x,y, '--', color='mediumslateblue', label="Angle maximal")
    plt.plot(x2,y2, '--', color='mediumslateblue')
    plt.xlabel('temps (sec)')
    plt.legend()
    plt.show()

### DEUXIEME FENETRE DE GRAPHIQUES : diagramme de phase  
def diagramme_phase():
    simulation()
    plt.figure("Graphique 2")
    plt.subplot(1,1,1)
    plt.plot(-theta,(-v_ang)*180/math.pi, color = "firebrick")
    plt.title(" Diagramme de phase ")
    plt.xlabel('angle inclinaison [°]')
    plt.ylabel('vitesse angulaire [°/s]')
    plt.show()
    


### EXCECUTION PREMIERE SIMULATION
if __name__ == "__main__":
    angle_maximal()
    simulation()    # exécute première simulation
    graphiques()    # dessine les graphiques correspondants
    diagramme_phase()
    print(theta[59000]) # affiche l'angle final




###########################
### DEUXIEME SIMULATION ###
###########################

# Calcul de l'angle d'inclinaison en fonction de la distance à laquelle l'objet est déplacé et de la masse de
# cet objet, fonctionne par itération: 
def calcul(d_moving, m_moving):
    """ pre : d_moving (float) est la distance à laquelle se trouve l'objet déplacé
              m_moving (float) est la masse de l'objet déplacé
        post : retourne theta, angle d'inclinaison.(en fonction également d'autres paramètres définis plus haut """
    
    # Fonctionne par ajout de petite itération jusqu'à trouver angle correct à 0.001 degrés près
    theta = 0.0001
    cgx = ((d_moving+math.sin(theta)*h_cable)*m_moving)/m_tot
    cpx = (((L_barge/2)**2)*math.tan(theta))/(3*hc)
    # Objectif : centre de gravité et centre de poussée alignés
    while abs(cgx-cpx) > 0.001:                    
        theta += 0.0001
        cgx = ((d_moving+math.sin(theta)*h_cable)*m_moving)/m_tot
        cpx = (((L_barge/2)**2)*math.tan(theta))/(3*hc)
    return theta
    

# Initialisation et création du tableau
dd = 0.01 # pas de distance [m]
end_ = 2  # [m]
d_axe = np.arange(0, end_ , dd)
moving_theta = np.empty_like(d_axe) 
d_ = [0, end_]
y_ = [a_max, a_max]


# Calcul de l'angle d'inclinaison avec déplacement et masse déplacée
def moving_charge():
    d_moving = 0
    for i in range(len(d_axe)-1):
        moving_theta[i+1] = math.degrees( calcul(d_moving, m_moving))
        d_moving += dd
    
# Masses à déterminer :
masse_moving1 = 0.100
masse_moving2 = 0.150
masse_moving3 = 0.200
masse_moving4 = 0.250
masse_moving5 = 0.500


### TROISIEME FENETRE DE GRAPHIQUES : Angle d'inclinaison avec différentes masses 
def graph_moving():
    
    plt.figure("Graphique 3")
    plt.subplot(1,1,1)
    global m_moving
    
    # Simulation avec différentes masses
    
    m_moving = masse_moving1
    moving_charge()
    plt.plot(d_axe,moving_theta, label= ("m = "+ str(masse_moving1) +" kg"),color = "lightseagreen", linewidth=3)
    plt.legend()
    
    m_moving = masse_moving2
    moving_charge()
    plt.plot(d_axe,moving_theta, label= ("m = "+ str(masse_moving2) + " kg"), color = 'steelblue', linewidth=3)
    plt.legend()

    m_moving = masse_moving3
    moving_charge()
    plt.plot(d_axe,moving_theta, label= ("m = "+ str(masse_moving3) +" kg" ), color = 'darkslateblue', linewidth=3)
    plt.legend()
    
    m_moving = masse_moving4
    moving_charge()
    plt.plot(d_axe,moving_theta, label= ("m = "+ str(masse_moving4) +" kg" ), color = 'mediumpurple', linewidth=3)
    plt.legend()
    
    m_moving = masse_moving5
    moving_charge()
    plt.plot(d_axe,moving_theta, label= ("m = "+ str(masse_moving5) +" kg" ), color = 'plum', linewidth=3)
    plt.legend()
    
    # Angle maximal
    plt.plot(d_,y_, '--', color='grey', label="angle maximal", linewidth=2)
    plt.legend()
    
    plt.title("Inclinaison par rapport à certaines masses")
    plt.xlabel("distance [m]")
    plt.ylabel("inclinaison [degrés]")
    plt.grid()
    plt.show()

### EXCECUTION DEUXIEME SIMULATION
if __name__ == "__main__":
    graph_moving()   # excécution simulations + graphiques
    



    
    
    
    
    
    
    
    
    
    
    
    
    
