from myBio import *

seqTest = "cgttgtatgcgtgtgacgtactgtgtgatgcggcgttgacggtgtgtcgtaggctatgagctagtggtcgaatgggcatgcgctcgcgacgctgcactcgcagcgtcgacgctgcagctgcatgcgttggggacagaatgcatgccgctgcagctgcagtcgc"



posStrand1StartCodon = []
posStrand2StartCodon = []
posStrand3StartCodon = []

negStrand1StartCodon = []
negStrand2StartCodon = []
negStrand3StartCodon = []



posStrand1StartCodon = startStopFinder(seqTest)
posStrand2StartCodon = startStopFinder(seqTest,2)
posStrand3StartCodon = startStopFinder(seqTest,3)

posStrand1StopCodon = startStopFinder(seqTest,codon="stop")
posStrand2StopCodon = startStopFinder(seqTest,2,codon="stop")
posStrand3StopCodon = startStopFinder(seqTest,3,codon="stop")

print(coordOrfFinder(posStrand1StartCodon,posStrand1StopCodon))
