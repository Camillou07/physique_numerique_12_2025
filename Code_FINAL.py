##########Projet de PHYS NUM : 24 Simulations cosmologiques avec l'approximation de Zel'Dovich####################



import numpy as np 
import matplotlib.pyplot as plt 
from scipy import integrate
from  scipy.integrate import quad


#_________________CONSTANTES______________________

H0=  2.27e-18 # constante de Hubble en sec-1 
omega_r_0=9e-5 # densité de radiation 
omega_m_0=0.3 # densité de matière totale 
omega_d_0=0.7 # densité d'énergie noire 


t_0= 4.35e17 # temps en sec =  aujourd'hui, noté avec 0 en cosmologie dans toutes les variables
t_b= 1e13   # en sec correspond à 300 000 ans = date de diffusion du fond diffus cosmologique

h= 4.4e13   # en sec, pas de temps 
a_b= 1e-4    # doit être petit le facteur d'échelle tend vers 0   


#friedman à t=0
G=6.67e-11 # constante de gravitation 
rho_m_0=3e-27 # densité moyenne de matière totale 
rho_r_0=8.3e-31 # densité  moyenne de radiation
rho_0=2.45e-27# densité moyenne d'énergie noire 



####### __________ETAPE 1 : TRACER FACTEUR ECHELLE AU COURS DU TEMPS _______________########





# Fonction  f_________________
def function_for_euler(a):
    return H0*a*np.sqrt(omega_r_0*(a**-4) + omega_m_0*(a**-3) + omega_d_0)



#Fonction euler optimisée____________
def euler_opt(t_b,t_0,a_b,h,function_for_euler):
    N=int((t_0-t_b)/h) 
    liste_x = np.linspace(t_b, t_0, N+1) # permet de créer déjà le tableau
    liste_a=np.zeros(N+1)
    liste_a[0]=a_b
    for i in range(N):
        liste_a[i+1]=liste_a[i]+h*function_for_euler(liste_a[i])
    return liste_x, liste_a
    

#on associe l1 aux tableau des temps et l2 à celui des valeurs de a______________
l1,l2=euler_opt(t_b,t_0,a_b,h,function_for_euler)


#Tracer le graphique___________________________________________
plt.figure(figsize=(8,5))
plt.plot(l1/(3.15e16),l2)
plt.xlabel("Time t ( Gyr) ")
plt.ylabel("Scale Factor a(t)")
plt.title("Plot of the scale factor evolution over time a(t)")
plt.show()



# Liste des temps pour lesquels on  veut connaître a(t)___________
Lt = [t_b, 1e14, 5e14, 1e15, 5e15, 1e16,5e16, 1e17, t_0]

# Affichage  des valeurs a(t)
for t in Lt:
    a_val = np.interp(t, l1, l2)
    print(f"t = {t:.3g} → a(t) = {a_val:.3g}")

#liste d'exemple  des valeurs de a qu'on veut étudier :
list_ex_a= [0.0001, 0.0116,0.0159,0.0206,0.0478,0.0733,0.208,0.331,1.03]






#######__________ ETAPE 2 : TRACER PARAMETRES DE DENSITÉ AU COURS DU TEMPS__________ #########




# Densités en fonction du facteur d'échelle l2 (a(t))
rho_m_t = rho_m_0 * l2**-3
rho_r_t = rho_r_0 * l2**-4
rho_d_t = rho_0 * np.ones_like(l2)  # énergie noire constante

# Hubble au temps t
H_t = H0 * np.sqrt(omega_r_0*l2**-4 + omega_m_0*l2**-3 + omega_d_0)

# Paramètres de densité Omega_i(t)
fonction_omega_m = rho_m_t / (3*H_t**2/(8*np.pi*G))
fonction_omega_r = rho_r_t / (3*H_t**2/(8*np.pi*G))
fonction_omega_d = rho_d_t / (3*H_t**2/(8*np.pi*G))



# Tracer la densité Oméga en fonction du temps 


plt.plot(l1/(3.15e16),fonction_omega_m, label="Omega_m", color="black")
plt.plot(l1/(3.15e16),fonction_omega_r, label="Omega_r", color="red")
plt.plot(l1/(3.15e16),fonction_omega_d, label="Omega_d", color="blue")
plt.xscale("log")   
plt.yscale("log")
plt.xlabel("Time (Gyr) ( in log) ")
plt.ylabel("Density Omega (in log) ")
plt.title("Evolution of omega as a function of the time t ")
plt.legend()  # affiche la légende
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()



#Tracer oméga en fonction du facteur d'échelle a 

plt.plot(l2,fonction_omega_m, label="Omega_m", color="black")
plt.plot(l2,fonction_omega_r, label="Omega_r", color="red")
plt.plot(l2,fonction_omega_d, label="Omega_d", color="blue")
plt.xscale("log")   
plt.yscale("log")
plt.xlabel(" Scale factor a (in log) ")
plt.ylabel(" Density Omega (in log)")
plt.title("Evolution of omega as a function of the scale factor a ")
plt.legend()  # affiche la légende
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()


# on peut tracer oméga en fonction du temps mais avec le temps pas en log

plt.plot(l1/(3.15e16),fonction_omega_m, label="Omega_m", color="black")
plt.plot(l1/(3.15e16),fonction_omega_r, label="Omega_r", color="red")
plt.plot(l1/(3.15e16),fonction_omega_d, label="Omega_d", color="blue")
#plt.xscale("log")   
plt.yscale("log")
plt.xlabel("Time(Gyr) ")
plt.ylabel(" Density Omega (in log)")
plt.title(" Evolution of omega as a function of the time t (t not in log)")
plt.legend()  # affiche la légende
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()



# on construit aussi rho en fonction du temps 

plt.plot(l1/(3.15e16),rho_m_t, label="rho_m", color="black")
plt.plot(l1/(3.15e16),rho_r_t, label="rho_r", color="red")
plt.plot(l1/(3.15e16),rho_d_t , label="rho_d", color="blue")
plt.xscale("log")   
plt.yscale("log")
plt.xlabel("Time (Gyr) (in log) ")
plt.ylabel(" Mean density rho (in log)")
plt.title("Evolution of rho as a function of the time t ")
plt.legend()  # affiche la légende
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()


# on construit aussi rho en fonction du facteur d'échelle a 

plt.plot(l2,rho_m_t, label="rho_m", color="black")
plt.plot(l2,rho_r_t, label="rho_r", color="red")
plt.plot(l2,rho_d_t , label="rho_d", color="blue")
plt.xscale("log")   
plt.yscale("log")
plt.xlabel("Scale Factor a (in log) ")
plt.ylabel("Mean density rho (in log) ")
plt.title(" Evolution of rho as a function of the scale factor ")
plt.legend()  # affiche la légende
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()






####### ETAPE 3 : "VIF DU PROJET" ####_____ Résolution des équations données dans le projet et établir la toile cosmique_________########





###________________________________________________ 1) On détermine D(t) par résolution d'une équa diff=> solution donnéee___________________


#on définit une fonction qui sera celle dans l'integrale_________
def fonction_for_rectangular_integration(y):
    y_safe = max(y, 1e-4)   # on remplace les valeurs trop petites par un petit seuil
    return (omega_m_0*(y_safe**-1) + omega_d_0*(y_safe**2) + (1-omega_m_0-omega_d_0))**-1.5 


#______verif______
for a_test in [0.1, 0.5, 1.0]:
    print(a_test, fonction_for_rectangular_integration(a_test))
#_________________


#Intégration rectangles sur un intervalle a,b______
def rectangular_integration(a,b,n):
    interval=(b-a)/n
    air_before=0
    for i in range(n):# boucle for sur tous les segments jusq'au n-ième
        air_before=air_before+interval*fonction_for_rectangular_integration(a+i*interval)
        if a+i*interval>=b:
            break    
    return air_before 

'''
#exemple d'application## de la fonction sur différentes valeurs de a aux extrémités de l'intervalle 
air_before_list = [] 
for a_val in l2:
        air_before = rectangular_integration(10**-4, a_val, 1000)
        air_before_list.append(air_before)
print(air_before_list)
'''

#fonction H_____________
def H(a):
    H_val = H0 * np.sqrt(omega_r_0*(a**-4) + omega_m_0*(a**-3) + omega_d_0)
    if H_val <= 0:
        print("Attention H(a) <= 0 pour a =", a)
    return H_val



#différentes valeurs de a 

#liste des temps interessants 
Lt = [t_b, 1e14, 5e14, 1e15, 5e15, 1e16,5e16, 1e17, t_0]
#liste des a interessants 
list_ex_a= [0.0001, 0.0116,0.0159,0.0206,0.0478,0.0733,0.208,0.331,1.03]


#fonction qui détermine liste des aires d'avant interessantes 
list_before=[]
for i in list_ex_a:
    alpha=rectangular_integration(0.1,i, 1000)
    list_before.append(alpha)
print(list_before)


#fonction qui donne une liste des valeurs de D associées
list_D_values=[]
for j,k in zip(list_ex_a,list_before) :
    beta=(H(j)/H0) * k
    list_D_values.append(beta)
print(list_D_values)

'''
#vérification du facteur de croissance


#_______verif______

D_list = [(H(a)/H0)*aire for a, aire in zip(l2, air_before_list)]
plt.plot(l2, D_list)
plt.xlabel("a")
plt.ylabel("D(a)")
plt.title("Facteur de croissance linéaire")
plt.show()
#__________________

'''


#____ définition des constantes___
L=2000#  MPC/h 'taille de Univers'
N=1000 # nombre de cellules par côtés
dx=L/N
dy=L/N



###_______________________________ 2) déterminer q (immobile dans espace initial) = tableau de valeurs définies selon les deux axes x et y ___________________

x=np.linspace(-L/2,L/2,N) # coordonnée espace dans une direction 
y=np.linspace(-L/2,L/2,N)
X,Y=np.meshgrid(x,y)# tableaux qui regroupe toutes les coordonnées sur x et sur y 
coords = np.stack((X, Y), axis=-1)

print(x)


## ______________________________3) On détermine le gradient du potentiel ________________________________


# fonctions random pour les distributions gausienne 

psi=np.random.normal(0,1,(N,N))  #tableau 2D 

#---- création des fonction phi à partir de psi----

phi=psi*-( (3/2)*H0**2* omega_m_0)**-1


# ----Transfo de fourier 2D----

F_k= np.fft.fft2(phi) # tableau 2D 


#----Grille des valeurs de k-----

kx=2*np.pi *np.fft.fftfreq(N,dx)
ky=2*np.pi *np.fft.fftfreq(N,dy)
kxv, kyv = np.meshgrid(kx, ky, indexing='ij')
k=np.sqrt(kxv**2+kyv**2) # normalisation car pas homogène tableau 
k[0,0]=np.inf


#----dériver dans espace de fourier donc multiplier par ik----

grad_phi_k_x = 1j * (kxv/k) * F_k # tableau 2D 
grad_phi_k_y = 1j * (kyv/k) * F_k


# -----on revient dans l'espace réel----

grad_phi_x = np.fft.ifft2(grad_phi_k_x).real # tablaeau 2D 
grad_phi_y = np.fft.ifft2(grad_phi_k_y).real 


#_____on empile juste les deux tableaux 2D____

gradient_potential= np.stack((grad_phi_x, grad_phi_y), axis=-1)

# ----multiplie les valeurs obtenues par D puis on ajoute le tableau de départ pour chaque coordonée____ 


#liste des temps interessants 
Lt = [t_b, 1e14, 5e14, 1e15, 5e15, 1e16,5e16, 1e17, t_0]
#liste des a interessants 
list_ex_a= [0.0001, 0.0116,0.0159,0.0206,0.0478,0.0733,0.208,0.331,1.03]


# plot final pour afficher les variations de densité 



#t=t_0
Xnew3 = coords + list_D_values[8]* gradient_potential   # tableau 2D avec deux coordonnées dans chaque case 
 
#On trace le résultat 

Xf3 = Xnew3[:,:,0].flatten()
Yf3 = Xnew3[:,:,1].flatten()

plt.figure(figsize=(6,5))
plt.hist2d(Xf3, Yf3, bins=500, cmap='plasma')
plt.colorbar(label='Density of the points')
plt.xlabel("Final Position X [Mpc/h]")
plt.ylabel(" Finale Position Y [Mpc/h]")
plt.title("Densité de matière simulée avec l'approximation de Zel'dovich ")
plt.show() 



#t=1e16
Xnew = coords + list_D_values[5]* gradient_potential   # tableau 2D avec deux coordonnées dans chaque case 
 
#On trace le résultat 

Xf = Xnew[:,:,0].flatten()
Yf = Xnew[:,:,1].flatten()

plt.figure(figsize=(6,5))
plt.hist2d(Xf, Yf, bins=500, cmap='plasma')
plt.colorbar(label='Density of the points')
plt.xlabel("Final Position X [Mpc/h]")
plt.ylabel(" Finale Position Y [Mpc/h]")
plt.title("Densité de matière simulée avec l'approximation de Zel'dovich ")
plt.show() 


