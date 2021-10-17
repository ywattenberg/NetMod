library(sna)

obsMat = as.matrix(read.csv("matrix.csv"))
cguRec1 <- cug.test(obsMat, grecip, cmod="edges", reps=3000)
cguRec2 <- cug.test(obsMat, grecip, cmod="dyad.census", reps=3000)
cguTrans <- cug.test(obsMat, gtrans, cmod="dyad.census", reps=3000)

print(cguRec1)
print(cguRec2)
print(cguTrans)
