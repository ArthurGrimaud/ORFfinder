# coding: utf-8
#Projet OBI 12/18
#Jeremy bulle / Arthur grimaud
#Librairie pour creation du programme ORF finder

import tkinter.filedialog as tkFileDialog
from tkinter import *
import os

allCoor = []
allOrfDict = []

def getAllOrfCoor(sequenceDic,sequenceName):
    """retourne une liste de listes de coordonnées pour chaques orfFilter
    arg:
    SequenceDic : dictionnaire de sequenceDic
    sequenceName : clef de la sequence
    """
    global allCoor
    sequence = sequenceDic[sequenceName] #au format [S+1,S+2,S+3,S-1,S-2,S-3] (S = strand)
    sequenceRev = Anti_sens(sequenceDic,sequenceName)
    for i in range(4): #pour les 3 ORF des brins positifs
        allCoor.append(coordOrfFinder(startStopFinder(sequence,i),startStopFinder(sequence,i,codon="stop")))
    for i in range(4):
        allCoor.append(coordOrfFinder(startStopFinder(sequenceRev,i),startStopFinder(sequenceRev,i,codon="stop")))
    return allCoor

def readNCBIFeatures(window):
    """
    Prend en argument un fichier ncbi "FEATURES" et renvoie une liste dictionaires
    au format [{id:,start:,stop;,name:}]

    arg: window : pour ouvrir la fenetre de selection de fichier
    """
    listDic = []
    dico = {}
    currdir = os.getcwd()
    tempdir = tkFileDialog.askopenfilename(parent=window, initialdir=currdir, title='Please select a directory')
    if len(tempdir) > 0:
        file = tempdir

    fastaSeq = []

    file=open(file,"r")
    for line in file:
        if ">" in line:
            dico = {}
            splited = line.split(":")
            dico["id"] = splited[0].split("|")[1]
            dico["start"] = int(splited[1])
            dico["stop"] = int(splited[2].split(" ")[0])
            listDic.append(dico)
    print(listDic)
    return listDic

def compare(orf1,orf2,window,sequence):
    """compare deux jeux de coordonées et identifie les coordonées identiques
    l'identifiant correspondant au coordonnées identiques est affichée dans une listBox tkInter

    arg: orf1 et orf2 : données au format [{id:AAA,start:BBB,stop;CCC,name:DDD},{...},...]
         window : fenetre Tkinter parent
         sequence : sequence dont sont issue les coordonnées des ORF
    """
    simi=Tk()

    listsimi = Listbox(simi)
    listsimi.pack()
    labelmax = Label(simi,text="The following NCBI ORF has been found",background="yellow")
    labelmax.pack()

    for o1 in orf1:
        for o2 in orf2:

            invposstart = len(sequence) - int(o1["start"]) # les positions sur NCBI sont notée sur le brin complementaire
            invposstart += 1                               #alors que dans le programme elles sont sur le brin inverse complementaire.
            invposstop = len(sequence) - int(o1["stop"])   #les coordonnées sont donc "inversées" pour les rendres comparable à celles de NCBI
            invposstop += 1

            print (len(sequence))
            len(sequence) - int(o1["stop"])
            if int(o1["start"]) == int(o2["start"]) and int(o1["stop"]) == int(o2["stop"]):
                ID = str(o1["id"])
                listsimi.insert(END, ID)
            if invposstart == int(o2["start"]) and invposstop == int(o2["stop"]):
                ID = str(o1["id"])
                listsimi.insert(END, ID)

    simi.mainloop()


def writeInCsv(allCoorF,fileName,createCSV=True):
    """
    arg: allCoorF : les coordonnées des orfFilter
        fileName : le nom du fichier crée

    return: la liste des dictionaires avec les données des ORF

    """
    allOrfDict = []
    number = 0
    for frame in allCoorF:
        for coor in frame:
            number += 1
            oneOrfDic = {}
            oneOrfDic["id"]=number
            oneOrfDic["start"]=coor[0]
            oneOrfDic["stop"]=coor[1]
            oneOrfDic["name"]="unknow"
            allOrfDict.append(oneOrfDic)

    if createCSV:
        file = str(fileName.get())
        fichier = open(file, "a")

        for i in allOrfDict:
            for key in i.keys():
                fichier.write(key+";")
            fichier.write("\n")
            for key in i.keys():
                fichier.write(str(i[key])+";")
            fichier.write("\n")

        fichier.close()
    return allOrfDict

def getLengths(allCoor):
    """
    Prend en argument les coordonnées des ORF et retourne
    une liste de la taille de ses orf avec leur positions associés
    """
    orfLength =[]
    orfPos=[]
    for frame in allCoor:
        for coor in frame:
            orfLength.append(coor[1]-coor[0])
            orfPos.append(coor)
    print(orfLength)
    print(orfPos)
    return (orfLength,orfPos)

def getReadingFrame(allCoor,seqCoor):
    """
    Retourne le cadre de lecture d'un ORF en fonction de ses coordonnées
    et de la liste de coordonnées dont il est issue
    """
    for frame in range(len(allCoor)):
        for pos in allCoor[frame]:
            if seqCoor == pos:
                if 0<=frame<=2:
                    return frame+1
                else:
                    return int((frame-2)/-1)

def getLongestORF(allCoor,sequence):
    """
    Affiche dans une fenetre tkInter
    le plus long orf de la liste de longueur d'orf
    """
    seqCoor = getLengths(allCoor)[1][getLengths(allCoor)[0].index(max(getLengths(allCoor)[0]))]
    sequenceDisplay = Tk()

    seqStart= Label(sequenceDisplay,text="Start Position :  "+str(seqCoor[0]))
    seqStart.pack()
    seqStop= Label(sequenceDisplay,text="Stop Position :  "+str(seqCoor[1]))
    seqStop.pack()
    seqStop= Label(sequenceDisplay,text="Frame :  "+str(getReadingFrame(allCoor,seqCoor)))
    seqStop.pack()
    seqLength= Label(sequenceDisplay,text="length :  "+str(seqCoor[1]-seqCoor[0])+" bp")
    seqLength.pack()

    if getReadingFrame(allCoor,seqCoor) >0:
        seqDNA = Label(sequenceDisplay,text=coorToSequence(seqCoor,sequence))
        seqDNA.pack()
        seqAA= Label(sequenceDisplay,text=translate(coorToSequence(seqCoor,sequence),convertGeneticTableFile("geneticCode")))
        seqAA.pack()
    else:
        seqDNA = Label(sequenceDisplay,text=coorToSequence(seqCoor,Anti_sens_withSequence(sequence)))
        seqDNA.pack()
        seqAA= Label(sequenceDisplay,text=translate(coorToSequence(seqCoor,Anti_sens_withSequence(sequence)),convertGeneticTableFile("geneticCode")))
        seqAA.pack()

    sequenceDisplay.mainloop()

def getTopORF(allCoor,sequence,top):
    if top > 1 :
        top = top/100
    for frame in allCoor :
        for coord in allCoor[frame]:
            long = allCoor[frame][1]-allCoor[frame][0] + 1
            liste_long.append((allCoor[0],allCoor[1],long))
    liste_long_decr = sorted(liste_long, lambda colonnes: colonnes[2], reverse = True)
    orf_tot = len(liste_long_decr)
    nb_target = int(round(top*orf_tot))
    for nb in (0,nb_target,1):
        liste_target.append(liste_long_decr[nb])

    return liste_target




def coorToSequence(coor,sequence):
    """Permet d'obtenir la sequence de nucleotide correspondante aux coordonnées d'un ORF
    """
    seqGene = ""
    for n in range(coor[0]-1,coor[1]+3):
        seqGene = seqGene + sequence[n]
    return seqGene

def orfFilter(orfCoorList,sequence,minLength = 10,maxLength = 1500):
    """Crée une nouvelle liste d'ORF filtée suivant une taille max et min à partir d'une liste de coordonnées d'ORF
    """
    orfCoorFiltered = []
    allOrfCoorFiltered =[]
    for list in orfCoorList:
        orfCoorFiltered = []
        for coor in list:
            if minLength < coor[1]-coor[0] < maxLength:
                orfCoorFiltered.append(coor)
        allOrfCoorFiltered.append(orfCoorFiltered)

    return allOrfCoorFiltered

def coordOrfFinder(startPos,stopPos):
    """Retourne une liste de coordonnées des orf en fonction des position des
    codons start et stop
    arg : startPos : liste des positions des codons start
          stopPos : liste des positions des condons stop
    return : liste de tuple des coordonnées
    """
    orfCoor = []
    for start in startPos:
        oneCoor = ()
        found = False
        for stop in stopPos:
            if start < stop and found == False:
                oneCoor = (start+1,stop+3)
                orfCoor.append(oneCoor)
                found = True
    return orfCoor

def startStopFinder(seq, readingFrame=1, codon = "start" ):
    """Pour un cadre de lecture donnée renvoie la liste des coordonnées des codons Start ou Stop
        arg : seq: sequence ADN
              readingframe : cadre de lecture pris en compte (1 par default)
              codon: codons recherchés (start ou stop)
    """
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
    """
    converti une table de code genetique au format ncbi au format .txt en listes python
    """
    geneticCodeTable =[]
    file=open(file,"r")
    for line in file:
        if "=" in line:
            cleanLine =""
            cleanLine =  str(line.split("= ")[1]).replace("\n","")
            geneticCodeTable.append(cleanLine)
    return(geneticCodeTable)



def translate(seq,geneticTable):
    """
    Pour chaques codons d'une sequence, determine l'acide aminée correspondant et retourne
    la sequence peptidique

    arg: seq : sequence à traduire
         geneticTable : code genetique utilisé pour la traduction
    """
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


def Is_codon_start(seq,pos):
    """retourne TRUE si le codon correspond a un codon start
    """
    seq = seq.upper()
    Codon = One_word(seq,pos,3)
    if Codon == "ATG" :
        return True
    else :
        return False

def Is_codon_stop(seq,pos):
    """retourne TRUE si le codon correspond a un codon stop"""
    seq = seq.upper()
    Codon = One_word(seq,pos,3)
    if Codon=="TAA" or Codon=="TAG" or Codon=="AGA" or Codon == "AGG":
        return True
    else :
        return False

def Is_gene(seq):
    """lecture donnée renvoie la liste des coordonnées des codons Start ou Stop"""
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

def Anti_sens_withSequence(seq):
    """retourne la seq d origine et la converti en brin anti sens, les sens de lecture es 5-3"""
    seq_anti = ""
    temp_anti = ""
    temp_anti = seq[::-1]
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
