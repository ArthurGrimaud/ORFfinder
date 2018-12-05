from myBio import *


dicoSequence = readFasta("humaMitoGenome(Short)")

orf = getAllOrfCoor(dicoSequence,">NC_012920.1 Homo sapiens mitochondrion, complete genome")
print(orf)

orfFiltered = orfFilter(orf,dicoSequence[">NC_012920.1 Homo sapiens mitochondrion, complete genome"],50,200)
print("FILTER \n",orfFiltered)
