# branchcutting
R package with implementation of the branchcutting algorithm for delimitation of species on a single phylogenetic tree.

### **Package installation**
It can be installed using devtools:

```
library(devtools)
devtools::install_github("onmikula/branchcutting", dependencies = TRUE)
```

### **Example usage**
When installed you can test its basic functionality:

```
library(branchcutting)
bcspecies <- bcut(dendromus)
bcspecies
plot(bcspecies, "species")
delim <- export(bcspecies)
head(delim)
```

Or, you can explore its functions as usual in R:

```
help(package="branchcutting")
?branchcutting
?plot.bcut
```
