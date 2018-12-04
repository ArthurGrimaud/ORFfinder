from myBio import *

seqTest = "cgttgtatgcgtgtgacgtactgtgtgatgcggcgttgacggtgtgtcgtaggctatgagctagtggtcgaatgggcatgcgctcgcgacgctgcactcgcagcgtcgacgctgcagctgcatgcgttggggacagaatgcatgccgctgcagctgcagtcgc"
seqTestNeg = Anti_sens(seqTest)



posStrand1StartCodon = startStopFinder(seqTest)
posStrand2StartCodon = startStopFinder(seqTest,2)
posStrand3StartCodon = startStopFinder(seqTest,3)

posStrand1StopCodon = startStopFinder(seqTest,codon="stop")
posStrand2StopCodon = startStopFinder(seqTest,2,codon="stop")
posStrand3StopCodon = startStopFinder(seqTest,3,codon="stop")

posStrand1orfCoor = coordOrfFinder(posStrand1StartCodon,posStrand1StopCodon)
posStrand2orfCoor = coordOrfFinder(posStrand2StartCodon,posStrand2StopCodon)
posStrand3orfCoor = coordOrfFinder(posStrand3StartCodon,posStrand3StopCodon)


print(posStrand1orfCoor)
print(posStrand2orfCoor)
print(posStrand3orfCoor)
