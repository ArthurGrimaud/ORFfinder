#Projet OBI 12/18
#jeremy bulle / arthur grimaud
#Librairie pour creation du programme ORF finder







def atgFinder(seq, readingFrame=1):
    strandStartCodon=[]
    for pos in range(readingFrame,len(seq),3):
        if Is_codon_start(seq,pos):
            strandStartCodon.append(pos)
    return strandStartCodon

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

def Write_fasta(Dict_seq,seq,filename,x=2):
    """ 4e arg defautl = 2, brin sens et anti-sens de copier, si x !=2, uniquement brin sens"""
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

def Compte_aa(prot):
    print("quel AA?")
    aa = input()
    nb_aa = 0
    aa = aa.upper()
    for i in protein:
        if i == aa :
            nb_aa = nb_aa + 1
    return nb_aa
    print ("il␣y␣a␣" , Compte_aa(protein) , "␣Cysteine(s)")

def Compte_all_aa(seq_prot):
    List_aa =["E","D","A","R","N","C","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    for i in List_aa :
        Nb_aa = 0
        for j in seq_prot:
            if i == j:
                Nb_aa = Nb_aa +1
        print(i ,"->", Nb_aa)


def Is_dna(adn):
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

    Dict = Read_fasta("seqalone.txt")
    print (Dict)
    bla = Dict[">gi|28302128|ref|NM_000518.4| Homo sapiens hemoglobin subunit beta (HBB), mRNA"]
    print(bla)
    Write_fasta(Dict,">gi|28302128|ref|NM_000518.4| Homo sapiens hemoglobin subunit beta (HBB), mRNA","seqalonetest.txt")
