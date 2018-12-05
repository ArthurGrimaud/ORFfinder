from myBio import *


dicoSequence = readFasta("humaMitoGenome(Short)")

# #display(posStrand1orfCoor,posStrand2orfCoor,posStrand3orfCoor,seqTest,negStrand1orfCoor,negStrand2orfCoor,negStrand3orfCoor)

getAllOrfCoor(dicoSequence,">NC_012920.1 Homo sapiens mitochondrion, complete genome")
