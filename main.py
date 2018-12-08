import tkinter.filedialog as tkFileDialog
from tkinter import *
import os
from myBio import*

allCoor = []
allCoorF=[]
sequence = ""
dicoSeq ={}
###################################################################################


def selectFile(window,listbox):
    global dicoSeq
    currdir = os.getcwd()
    tempdir = tkFileDialog.askopenfilename(parent=window, initialdir=currdir, title='Please select a directory')
    if len(tempdir) > 0:
        file = tempdir
    dicoSeq = readFasta(file)
    for key in dicoSeq.keys():
        listbox.insert(END, key)


def printCoor(listbox,dicoSeq,window,listCoor): #Affiche les coordonnées dans une listeBox tkinter
    global sequence
    global allCoor
    sequenceName = listbox.get(listbox.curselection())
    sequence = dicoSeq[sequenceName]
    allCoor = getAllOrfCoor(dicoSeq,sequenceName)
    for i in allCoor:
        for j in i:
            listCoor.insert(END, j)

def printFilteredCoor(allCoor,min,max,listCoorF):
    global allCoorF
    minVal = int(min.get())
    maxVal = int(max.get())
    allCoorF = orfFilter(allCoor,sequence,minVal,maxVal)
    listCoorF.delete(0, last="end")
    for i in allCoorF:
        for j in i:
            listCoorF.insert(END, j)



def getReadingFrame(allCoor,seqCoor):
    for frame in range(len(allCoor)):
        for pos in allCoor[frame]:
            if seqCoor == pos:
                if 0<=frame<=2:
                    return frame+1
                else:
                    return int((frame-2)/-1)

def displaySelectedOrf(listCoorF,sequence,allCoorF):
    seqCoor = listCoorF.get(listCoorF.curselection())
    sequenceDisplay = Tk()

    seqStart= Label(sequenceDisplay,text="Start Position :  "+str(seqCoor[0]))
    seqStart.pack()
    seqStop= Label(sequenceDisplay,text="Stop Position :  "+str(seqCoor[1]))
    seqStop.pack()
    seqStop= Label(sequenceDisplay,text="Frame :  "+str(getReadingFrame(allCoorF,seqCoor)))
    seqStop.pack()
    seqLength= Label(sequenceDisplay,text="length :  "+str(seqCoor[1]-seqCoor[0])+" bp")
    seqLength.pack()

    if getReadingFrame(allCoorF,seqCoor) >0:
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



####################################################################################""

window = Tk()
window.configure(background='grey')


select = Button(window, text="Select fasta file (1)", command=lambda : selectFile(window,listbox),bg="firebrick3")
select.grid(row=1,column=0)

viewCoor = Button(window, text="Find ORF(2)", command=lambda : printCoor(listbox,dicoSeq,window,listCoor),bg="firebrick3")
viewCoor.grid(row=0,column=1)

selectSeq = Button(window, text="Filter ORF(3)", command=lambda : printFilteredCoor(allCoor,min,max,listCoorF),bg="firebrick3")
selectSeq.grid(row=0,column=2)

selectSeq = Button(window, text="Display(4)", command=lambda : displaySelectedOrf(listCoorF,sequence,allCoorF),bg="firebrick3")
selectSeq.grid(row=0,column=3)

selectSeq = Button(window, text="Store ORF in .cvs (4)", command=lambda : writeInCsv(allCoorF,fileName) ,bg="firebrick3")
selectSeq.grid(row=0,column=4)

labelmax = Label( window,text="Max Length (Nucleotides):",background="yellow")
labelmax.grid(row=3,column=1)
labelmin = Label( window,text="Min Length (Nucleotides):",background="orange")
labelmin.grid(row=2,column=1)

max = Entry( window,background="yellow")
max.insert(END, '100000000')
max.grid(row=3,column=2)

min = Entry( window,background="orange")
min.insert(END, '0')
min.grid(row=2,column=2)

fileName = Entry( window,background="green")
fileName.insert(END, "ORF")
fileName.grid(row=1,column=4)

listbox = Listbox(window)
listbox.grid(row=1,column=1)

listCoor = Listbox(window)
listCoor.grid(row=1,column=2)

listCoorF = Listbox(window)
listCoorF.grid(row=1,column=3)


window.mainloop()
