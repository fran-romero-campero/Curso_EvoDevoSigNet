################### COMPARACIÓN DE ÁRBOLES POR TOPOLOGÍA #######################
# Vamos a usar el paquete TreeDist
install.packages("TreeDist")
library(TreeDist)
library(ape)

#  De forma rápida, se pueden generar dos árboles distintos y 
# calcular su distancia por Robinson-Foulds
set.seed(725)

treeA <- rtree(8)
plot(treeA)
treeB <- rtree(8)
plot(treeB)

RobinsonFoulds(treeA, treeB)


# Controlando más los árboles, se va a usar el mismo ejemplo que en TreeDist
# https://cran.r-project.org/web/packages/TreeDist/vignettes/Robinson-Foulds.html
tree1 <- ape::read.tree(text='(A, (B, (C, (D, (E, (F, (G, H)))))));')
tree2 <- ape::read.tree(text='(A, (B, (C, (D, (E, (G, (F, H)))))));')
tree3 <- ape::read.tree(text='(A, (B, (C, (E, (D, (F, (G, H)))))));')


######### Robinson-Foulds
# Distancia entre el árbol 1 y 2
par(mfrow=c(1,2))
plot(tree1)
plot(tree2)

RobinsonFoulds(tree1, tree2) # Distancia de 2 con RB

# Distancia entre el árbol 1 y 2
par(mfrow=c(1,2))
plot(tree1)
plot(tree3)

RobinsonFoulds(tree1, tree2) # También distancia de 2 con RB

# Sin embargo, la división en conflicto en el segundo par porta más información, 
# por lo que la distancia debería ser mayor

VisualizeMatching(InfoRobinsonFoulds, tree1, tree2, 
                  Plot = TreeDistPlot, prune = 12)

VisualizeMatching(InfoRobinsonFoulds, tree1, tree3, 
                  Plot = TreeDistPlot, prune = 8)


######### Information-corrected Robinson-Foulds

# Por tanto, una distancia Information-corrected Robinson-Foulds puede aportar
# más información en casos como este

InfoRobinsonFoulds(tree1, tree2) # Distancia de 6.92
InfoRobinsonFoulds(tree1, tree3) # Distancia de 11.06

# De esta manera, se penaliza más diferencias en divisiones más centrales,
# diviisiones pares que separen, por ejemplo, 4 hojas de otras 4 hojas


# No obstante, también presenta limitaciones, por ejemplo, que satura 
# rápido, pudiendo alcanzar la mayor distancia moviendo una única hoja
# (H)

tree1 <- ape::read.tree(text='(1, (2, (3, (4, (5, (6, (7, 8)))))));')
tree2 <- ape::read.tree(text='(8, (1, (2, (3, (4, (5, (6, 7)))))));')

VisualizeMatching(RobinsonFouldsMatching, tree1, tree2, Plot = TreeDistPlot)


# Para esto se han desarrollado distancias Robinson-Foulds generalizadas, 
# que trazan un mapeo entre los splits de ambos árboles y calculan su similitud,
# aunque no sean idénticas

tree1 <- ape::read.tree(text='((A, B), ((C, (D, E)), (F, (G, (H, I)))));')
tree2 <- ape::read.tree(text='((A, B), ((C, D, (E, I)), (F, (G, H))));')

VisualizeMatching(RobinsonFouldsMatching, tree1, tree2, Plot = TreeDistPlot)
RobinsonFoulds(tree1, tree2) # distancia de 9, cerca de la máxima que es 11 
# para esta topología

# Muy penalizado, ya que simplemente eliminando la hoja I salen dos árboles 
# idénticos excepto por la resolución de un nodo interno
VisualizeMatching(RobinsonFouldsMatching,
                  drop.tip(tree1, 'I'),
                  drop.tip(tree2, 'I'))


####### Shared Phylogenetic Information

# Si usamos una distancia generalizada (Shared phylogenetic information) que
# no asigne probabilidades de 0 a las que no son idénticas, observamos la 
# similitud oculta en el caso anterior
VisualizeMatching(SharedPhylogeneticInfo, tree1, tree2)
SharedPhylogeneticInfo(tree1, tree2) # 12.32 bits de distancia

# Podemos calcular el contenido de información total de un árbol:
SplitwiseInfo(tree1)

# A pesar de esto, la Shared phylogenetic information asigna valores de 0 a 
# pares de splits incompatibles entre sí, que puede ser indeseable en ciertas 
# condiciones

tree1 <- ape::read.tree(text='((F,(G,(H,(I,J)))),(E,(D,(C,(B,A)))));')
tree2 <- ape::read.tree(text='((F,(G,(H,(I,A)))),(E,(D,(C,(B,J)))));')
SharedPhylogeneticInfo(tree1, tree2)
VisualizeMatching(SharedPhylogeneticInfo, tree1, tree2,
                  Plot = TreeDistPlot, matchZeros = FALSE, prune = c(5, 18))


#######   Mutual Clustering Information

# Como alternativa, Mutual Clustering Information reconoce la similitud en 
# la estructura del árbol incluso cuando cada posible emparejamiento de splits
# entra en conflicto
VisualizeMatching(MutualClusteringInfo, tree1, tree2,
                  Plot = TreeDistPlot, matchZeros = FALSE, prune = c(9, 18))

MutualClusteringInfo(tree1, tree2) # 1.34

# Como antes, la Mutual Clustering Information total de un árbol se puede calcular:
ClusteringEntropy(tree1) # 6.15



