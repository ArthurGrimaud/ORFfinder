import tkinter.filedialog as tkFileDialog
from tkinter import *
import os
from myBio import*




def selectFile():
    currdir = os.getcwd()
    tempdir = tkFileDialog.askopenfilename(parent=window, initialdir=currdir, title='Please select a directory')
    if len(tempdir) > 0:
        print ("You chose", tempdir)
        file = tempdir
    dicoSeq = readFasta(file)
    listeSeq(dicoSeq)

def listeSeq(dicoSeq):
    listbox = Listbox(window)
    listbox.pack(side=BOTTOM)


    for key in dicoSeq.keys():
        listbox.insert(END, key)

    selectSeq = Button(window, text="Find ORF", command=lambda : getsequence(listbox,dicoSeq))
    selectSeq.pack(side=BOTTOM)

def getsequence(listbox,dicoSeq):
    sequenceName = listbox.get(listbox.curselection())
    sequence = dicoSeq[sequenceName]
    print(sequence)
    allCoor = getAllOrfCoor(dicoSeq,sequenceName)
    importseq = Label(window,text ="Done")
    importseq.pack()

############################""#main###############################""

window = Tk()
select = Button(window, text="Selectionner un fichier fasta", command=lambda : selectFile())
select.pack()

window.mainloop()
