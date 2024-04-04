# Film
Le but de la classe Film est de générer un film magnétique représenté par une matrice et de calculer son énergie totale dans le but de la minimiser par la suite. Chaque élément d'une matrice Film représente un grain magnétique du film et un grain peut être repéré par ses coordonnées $(x,y)$ dans la matrice ou par son indice $i$. 
## Attributs
- nx : nombre de colonnes de grains dans la matrice Film.
- ny : nombre de lignes de grains dans la matrice Film.
- M_norm : module de l'aimantation des grains (la même partout car il s'agit d'une aimantation intrinsèque aux grains qui sont tous identiques).
- theta_def : direction des aimantations par défaut.
- theta_data : tableau de nx*ny cases ayant pour valeur la direction de l'aimantation. Ainsi, chaque case représente un grain magnétique et prend pour valeur la direction d'aimantation de ce grain.
- name : nom que l'on souhaite donner au Film. C'est le seul attribut publique car on veut pouvoir l'utiliser lorsque l'on crée une copie d'un Film.
## Constructeurs et destructeurs
- Film(name, nx, ny, M_norm, theta_def) : matrice de nx colonnes et ny lignes, toutes ayant pour valeur theta_def.
- Film(name, nx, ny, M_norm) : constructeur aléatoire. Matrice de nx colonnes et ny lignes avec des valeurs de theta prises aléatoirement entre $0$ et $2\pi$. Si on ne rentre pas une valeur de M_norm, il prendra 1 pour valeur par défaut.
- ~Film() : destructeur.
## Opérateurs
Soit un film F.
- F(i) : accès à la valeur de l'élément $i$ de F.
- F(x,y) : accès à la valeur de l'élément de la colonne x et la ligne y de F.
- F(i)=new_val : changer la valeur de l'élément $i$ de F en new_val.
- F(x,y)=new_val : changer la valeur de l'élément de la colonne x et la ligne y de F en new_val.  

Remarque : on prend des conditions aux limites périodiques.
## Méthodes
- find_i(x,y) : retourne l'indice i correspondant à la colonne x et la ligne y.
- find_xy(i) : retourne les coordonnées $(x,y)$ correspondant à l'indice $i$.
- show() : affiche le Film sous forme de matrice.
- write(fich_name) : écrit dans un fichier les valeurs de theta. fich_name est le nom que l'on veut donner au fichier de valeurs, il n'est pas obligatoire de le rentrer (il prendra le nom "theta_data" par défaut).
- E_ech(A) : calcule l'énergie d'échange totale du Film. A est la valeur de la constante d'échange.  
Pour l'instant on est à deux dimensons : $m_x=\sin\theta$, $m_y=0$ et $m_z=\cos\theta$.  
Comme on considère un film plan $(x,y)$, la composante du gradient selon $z$ est nulle.
- E_d(Nxx, Nyy, Nzz) : calcule l'énergie du champ démagnétisant totale du Film. Nxx, Nyy et Nzz sont les éléments du tenseur $[N_d]$. On doit avoir Nxx+Nyy+Nzz=1.  
- E_ani(K) : calcule l'énergie d'anisotropie totale du Film. K est la constante d'anisotropie.
- E_z(phi, H_a_norm) : calcule l'énergie totale due à la présence d'un champ magnétique extérieur de norme H_a_norm et faisant un angle phi avec l'axe $O_z$.
On est aussi à deux dimensions pour $\vec{H_a}$ : $H_{a_x}=H_a\sin(\varphi)$, $H_{a_y}=0$ et $H_{a_z}=H_a\cos(\varphi)$.
- E_tot(A, Nxx, Nyy, Nzz, K, phi, H_a_norm) : calcule l'énergie totale du Film.
- E_tot(A, Nxx, Nyy, Nzz, K) : calcule l'énergie totale du Film hors champ magnétique extérieur.

Remarque : on calcule les énergies totales en sommant sur tous les éléments de la matrice Film.

- write(fich_name, A, Nxx, Nyy, Nzz, K, phi, H_a_norm) : même fonction que la première fonction write mais note en plus les constantes du système au début du fichier.

# Magnetization_Plot
Ce fichier Python contient deux fonctions pour faire des plots des directions d'aimantation à partir des fichiers de données des angles générés à l'aide de la fonction write de la classe Film.
- f(theta, title=' ') : pour faire des plots des directions d'aimantation. theta est le fichier de données à afficher et title et le titre qu'on veut éventuellement donner à la figure.
- gif(files) : pour créer des GIF montrant l'évolution des aimantations d'un film magnétique. files est une liste contenant tous les fichiers de données des angles qu'on veut mettre dans le GIF.
