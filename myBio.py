# coding: utf-8
#Projet OBI 12/18
#jeremy bulle / arthur grimaud
#Librairie pour creation du programme ORF finder


def getAllOrfCoor(sequenceDic,sequenceName): #retourne une liste de listes de coordonnées
    allCoor = []
    sequence = sequenceDic[sequenceName]             #au format [S+1,S+2,S+3,S-1,S-2,S-3] (S = strand)
    sequenceRev = Anti_sens(sequenceDic,sequenceName)
    for i in range(1,4): #pour les 3 ORF des brins positifs
        print(coordOrfFinder(startStopFinder(sequence,i),startStopFinder(sequence,i,codon="stop")))
        allCoor.append(coordOrfFinder(startStopFinder(sequence,i),startStopFinder(sequence,i,codon="stop")))
    for i in range(1,4):
        allCoor.append(coordOrfFinder(startStopFinder(sequenceRev,i),startStopFinder(sequenceRev,i,codon="stop")))
        print(coordOrfFinder(startStopFinder(sequenceRev,i),startStopFinder(sequenceRev,i,codon="stop")))

    return allCoor





def display(S1coor,S2coor,S3coor,sequence,S4coor,S5coor,S6coor):
    allPosStrandCoor = [S1coor,S2coor,S3coor]
    allNegStrandCoor = [S4coor,S5coor,S6coor]
    seqLine = ""
    for strandCoor in allPosStrandCoor:
        for j in range(len(sequence)):
            done = False
            for coor in strandCoor:
                if coor[0] <= j <= coor[1]+2 and done == False:
                    seqLine = seqLine + sequence[j]
                    done = True
                elif done == False:
                    seqLine = seqLine + "-"
                    done = True
        print(seqLine)
        seqLine = ""
    print(sequence)
    # for strandCoor in allNegStrandCoor:
    #     for j in range(len(sequence)):
    #         done = False
    #         for coor in strandCoor:
    #             if coor[0] <= j <= coor[1]+2 and done == False:
    #                 seqLine = seqLine + sequence[j]
    #                 done = True
    #             elif done == False:
    #                 seqLine = seqLine + "-"
    #                 done = True
    #     print(seqLine)
    #     seqLine = ""





def coorToSequence(coor,sequence):
    seqGene = ""
    for n in range(coor[0],coor[1]+3):
        seqGene = seqGene + sequence[n]
    return seqGene




def orfFilter(orfCoorList,sequence,minLength = 10,maxLength = 1500):
    orfCoorFiltered = []
    for coor in range(len(orfCoorList)):
        if minLength < orfCoorList[coor][1]-orfCoorList[coor][0] < maxLength:
            orfCoorFiltered.append(orfCoorList[coor])

    return orfCoorFiltered


def coordOrfFinder(startPos,stopPos):
    orfCoor = []
    for start in startPos:
        oneCoor = ()
        found = False
        for stop in stopPos:
            if start < stop and found == False:
                oneCoor = (start,stop)
                orfCoor.append(oneCoor)
                found = True
    return orfCoor



def startStopFinder(seq, readingFrame=1, codon = "start" ):
    strandCodon=[]
    if codon == "start":
        for pos in range(readingFrame,len(seq),3):
            if Is_codon_start(seq,pos):
                strandCodon.append(pos)
        return strandCodon
    elif codon =="stop":
        for pos in range(readingFrame,len(seq),3):
            if Is_codon_stop(seq,pos):
                strandCodon.append(pos)
        return strandCodon
    else:
        print("invalid codon argument (please enter: start or stop)")

def readFasta(file):
    """ lit un fichier fasta et crée un dict avec le nom des genes en tant que clé
    les nom des genes doit etre précédé par un >,
    tout ce qui suit après le premier retour à la ligne est considéré comme seq"""
    fastaSeq = {}

    file=open(file,"r")
    for line in file:
        if ">" in line:
            nameLine = line.replace("\n","")
            fastaSeq[nameLine]= ""
        if ">" not in line:
            fastaSeq[nameLine] = fastaSeq[nameLine] + line.replace("\n","")
    return fastaSeq

def Write_fasta(Dict_seq,seq,filename,x=2):
    """ ajoute à la fin du fichier, le contenu du dict de gène. la clé doit etre
    un nom de gène et le contenu une seq au format fasta
    le 4e arg(x), valeur par défault = 2, brin sens et anti-sens de copier
    si x !=2, uniquement brin sens"""
    chemin = "./"+filename
    fic = open(chemin)
    contenu_fic = fic.read()
    if x == 2 :
        brin_anti = Anti_sens(Dict_seq,seq)
        new_contenu = contenu_fic + "\n" + seq + '\n'+ Dict_seq[seq] + "\n" + seq + " " + "anti-sens" +'\n'+ brin_anti
    if x is not 2 :
        new_contenu = contenu_fic + "\n" + seq + '\n'+ Dict_seq[seq]
    fic.write(new_contenu)
    fic.close()



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

def enter_seq():
    adn = input("Entrez␣la␣chaine␣:␣")
    for i in range (len(adn) ):
        print("i:␣", i, "->" , adn[i])

def Compte_aa(seq_prot):
    """compte le nombre d'apparition d'un AA donné dans la séquence d'une prot"""
    print("quel AA?")
    aa = input()
    nb_aa = 0
    aa = aa.upper()
    for i in seq_prot:
        if i == aa :
            nb_aa = nb_aa + 1
    return nb_aa
    print ("il␣y␣a␣" , Compte_aa(protein) , "␣Cysteine(s)")

def Compte_all_aa(seq_prot):
    "affiche le nombre d'apparition de chaque AA dans une seq prot"
    List_aa =["E","D","A","R","N","C","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    for i in List_aa :
        Nb_aa = 0
        for j in seq_prot:
            if i == j:
                Nb_aa = Nb_aa +1
        print(i ,"->", Nb_aa)


def Is_dna(adn):
    """retourne tru si la sequence entree est une sequence d adn """
    seq = ""
    for i in adn:
        i = i.upper()
        seq = seq + i
    flag = 0
    for i in seq :
        if i=="A" or i=="T" or i=="G" or i=="C":
            flag = flag
        else :
            flag = flag +1
    if flag > 0 :
        return False
    else :
        return True

def One_word(seq, start,wlen):
    """seq est une string entrée, start est l'index de départ du mot
    dans la string
    wlen est un entier qui détermine la longueure du mot recherché
    One_word retourne le mot trouvé"""
    s=""
    j = 0
    num = start
    while j < wlen :
        try :
            s = s + seq[num]
            j = j +1
            num = num +1
        except IndexError :
            return s
    return s

def Count_word(seq, word):
    """ seq est une string entrée, word est une string qui est recherchéen dans
    la seq
    Count_word retourne le nombre d'occurence de word dans la séquence et dans
    tous les cadres de lecture"""
    word = word.upper()
    wlen = len(word)
    flag = 0
    for i in range(len(seq)):
        Codon = One_word(seq,i,wlen)
        if Codon == word :
            flag = flag +1
    print(word, "est présent", flag, "fois dans la séquence")

def Is_codon_start(seq,pos):
    seq = seq.upper()
    Codon = One_word(seq,pos,3)
    if Codon == "ATG" :
        return True
    else :
        return False

def Is_codon_stop(seq,pos):
        seq = seq.upper()
        Codon = One_word(seq,pos,3)
        if Codon=="TAA" or Codon=="TAG" or Codon=="TGA":
            return True
        else :
            return False

def Is_gene(seq):
    for i in range(len(seq)):
        if Is_codon_start(seq,i)==True:
            new_seq = seq[i:len(seq)]
            for j in range(0,len(new_seq),3):
                if Is_codon_stop(seq,j)==True:
                    print ("il y a un gene, il commence en position", i,"et fini en", j+3)
                    seq_gene = new_seq[0:j]
                    print(seq_gene)

def Anti_sens(dict_seq,seq_name):
    """retourne la seq d origine et la converti en brin anti sens, les sens de lecture es 5-3"""
    seq_anti = ""
    temp_anti=""
    temp = dict_seq[seq_name]
    temp_anti = temp[::-1]
    for n in temp_anti :
        if n == "A" :
            n = "T"
        elif n == "T":
            n = "A"
        elif n == "G":
            n = "C"
        elif n == "C":
            n = "G"
        seq_anti = seq_anti + n
    return seq_anti


if __name__ == '__main__':

    test = "AATAACTT"
    Count_word(test,"AA")
