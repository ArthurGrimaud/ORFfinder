#Projet OBI 12/18
#jeremy bulle / arthur grimaud
#Librairie pour creation du programme ORF finder





def readFasta(file):
    fastaSeq = {}

    file=open(file,"r")
    for line in file:
        if ">" in line:
            nameLine = line.replace("\n","")
            fastaSeq[nameLine]= ""
        if ">" not in line:
            fastaSeq[nameLine] = fastaSeq[nameLine] + line.replace("\n","")
    return fastaSeq

def writeFasta(dico,fileN):
    file = open(fileN, "w+")
    for item in dico.items():
        for i in item :
            print(i)
            file.write(str(i)+'\n')
        file.write('\n')



def convertGeneticTableFile(file):
    geneticCodeTable =[]
    file=open(file,"r")
    for line in file:
        if "=" in line:
            cleanLine =""
            cleanLine =  str(line.split("= ")[1]).replace("\n","")
            geneticCodeTable.append(cleanLine)
    return(geneticCodeTable)



def translate(seq,geneticTable):
    aaSeq =""
    interAA1 =[]
    interAA2 =[]
    indexAA3 = 0
    for codon in range(0,len(seq)-3,3):
        interAA1 =[]
        for ia in range(len(geneticTable[2])): #Premier nucleotide du codon
            if geneticTable[2][ia] == seq[codon]:
                interAA1.append(ia)
                interAA2 =[]
                for ib in range(interAA1[0],interAA1[len(interAA1)-1]): #Second nucleotide du codon
                    if geneticTable[3][ib] == seq[codon+1]:
                        interAA2.append(ib)
                        for ic in interAA2: #Troisieme nucleotide du codon
                            if geneticTable[4][ic] == seq[codon+2]:
                                indexAA3 = ic

        aaSeq = aaSeq + str(geneticTable[0][indexAA3])

    return aaSeq
