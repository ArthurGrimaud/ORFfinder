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

def displaySelectedOrf(listCoorF,sequence):
    seqCoor = listCoorF.get(listCoorF.curselection())
    sequenceDisplay = Tk()
    seqDNA = Label(sequenceDisplay,text=coorToSequence(seqCoor,sequence))
    seqDNA.pack()
    seqAA= Label(sequenceDisplay,text=translate(coorToSequence(seqCoor,sequence),convertGeneticTableFile("geneticCode")))
    seqAA.pack()
    sequenceDisplay.mainloop()





####################################################################################""

window = Tk()


select = Button(window, text="Selectionner un fichier fasta (1)", command=lambda : selectFile(window,listbox))
select.grid(row=0,column=0)

viewCoor = Button(window, text="Find ORF(2)", command=lambda : printCoor(listbox,dicoSeq,window,listCoor))
viewCoor.grid(row=0,column=1)

selectSeq = Button(window, text="Filter ORF(3)", command=lambda : printFilteredCoor(allCoor,min,max,listCoorF))
selectSeq.grid(row=0,column=2)

selectSeq = Button(window, text="Display", command=lambda : displaySelectedOrf(listCoorF,sequence))
selectSeq.grid(row=0,column=3)

labelmax = Label( window,text="Max")
labelmax.grid(row=3,column=1)
labelmin = Label( window,text="Min")
labelmin.grid(row=2,column=1)

max = Entry( window)
max.insert(END, '100000000')
max.grid(row=3,column=2)
min = Entry( window)
min.insert(END, '0')
min.grid(row=2,column=2)

listbox = Listbox(window)
listbox.grid(row=1,column=1)

listCoor = Listbox(window)
listCoor.grid(row=1,column=2)

listCoorF = Listbox(window)
listCoorF.grid(row=1,column=3)


window.mainloop()
