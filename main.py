from myBio import *

seqTest = "cgttgtatgcgtgtgacgtactgtgtgatgcggcgttgacggtgtgtcgtaggctatgagctagtggtcgaatgggcatgcgctcgcgacgctgcactcgcagcgtcgacgctgcagctgcatgcgttggggacagaatgcatgccgctgcagctgcagtcgc"



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

print(orfFilter(posStrand1orfCoor,seqTest))
print(orfFilter(posStrand2orfCoor,seqTest))
print(orfFilter(posStrand3orfCoor,seqTest))

posStrand3orfCoorFiltered = orfFilter(posStrand3orfCoor,seqTest)

# print(coorToSequence(posStrand3orfCoorFiltered[0],seqTest))

display(posStrand1orfCoor,posStrand2orfCoor,posStrand3orfCoor,seqTest)
