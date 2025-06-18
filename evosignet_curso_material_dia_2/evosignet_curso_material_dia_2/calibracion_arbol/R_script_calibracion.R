####################### CALIBRACIÓN DE ÁRBOLES ##############################
# Vamos a usar los paquetes ape y ggtree para calibrar árboles filogenéticos
# usando evidencia fósil

# install.packages("ape")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")
library(ape)
library(ggtree)

# Antes, hemos guardado el árbol SpeciesTree_rooted.txt de Orthofinder como
# tree_evosignet.nwk
tree <- read.tree("tree_evosignet.nwk")
is.binary(tree)
is.ultrametric(tree)

# Al ser un árbol de prueba, las longitudes de las ramas no son las que deberían,
# vamos a asignar un valor de 1 a cada una
length(tree$edge.length)
tree$edge.length <- rep(1, 20)
write.tree(tree, "new_revised_species_tree_complete_lenghts.txt")

# Para ver la notación de los nodos internos, vamos a hacer una representación
# rápida con ggtree
new_tree <- read.tree("new_revised_species_tree_complete_lenghts.txt")
p <- ggtree(new_tree) + geom_text(aes(label=node), hjust=-.05) + geom_tiplab()
plot(p)

# Con esta información y la información sobre el registro fósil del artículo
# https://www.nature.com/articles/s41467-021-22044-z ya podemos hacer la calibración
# Primero, hay que crear tres vectores, que, en el mismo orden, muestren los
# tiempos máximos y mínimos de aparición para los distintos nodos internos
# que conozcamos en nuestro árbol. Podemos fijar el age.max del nodo raíz como
# la máxima profundidad que vayamos a permitir en el árbol
node <- c(12, 21, 16, 19, 20, 14, 13)
age.min <- c(1600, 700, 470, 130, 124, 745, 729)
age.max <- c(1900, 1105, 540, 385,  130, 1350, 1105)

# Este último vector no se usa, pero hay que incluirlo para evitar errores de dimensiones
soft.bounds <- rep(TRUE, length(node))

# Creamos una tabla con los cuatro vectores
mycalibration <- data.frame(node, age.min, age.max, soft.bounds) 

# Y ya aplicamos el modelo relaxed
set.seed(355)
mytimetree <- chronos(new_tree, lambda = 4, model = "relaxed", 
                      calibration = mycalibration, control = chronos.control())

set.seed(355)
mytimetree2 <- chronos(new_tree, lambda = 9, model = "relaxed", 
                      calibration = mycalibration, control = chronos.control())

# Para algunos valores de lambda, el modelo puede no converger
set.seed(355)
mytimetree3 <- chronos(new_tree, lambda = 1, model = "relaxed", 
                       calibration = mycalibration, control = chronos.control())

# Para comprobar las diferencias, podemos hacer un gráfico de las longitudes
# del árbol usando diferentes gammas
# Para 4 y 9:
plot(mytimetree$edge.length,mytimetree2$edge.length, 
     xlim = c(0,2000), ylim = c(0,2000),
     pch=21,bg="grey",cex=1.2,bty="n",
     xlab=expression(paste(lambda,"= 1")),
     ylab=expression(paste(lambda,"= 0.1")))
lines(c(0,2000),c(0,2000))
legend("topleft","1:1 line",lty="solid",bty="n")
grid()
title(main=expression(paste(
  "Comparison of edge lengths with two different ",
  lambda," values")))

# Vemos que las longitudes son prácticamente idénticas

# Por último, representamos el árbol datado
plot_tree <- ggtree(mytimetree) + geom_tiplab()
plot_tree2 <- revts(plot_tree)
plot_tree3 <- plot_tree2 + theme_tree2(bgcolor = "transparent") + 
  xlim_expand(c(-2000, 700), panel = "Tree")
plot_tree3

# Guardamos el árbol
write.tree(mytimetree, "calibrated_tree.nwk")
