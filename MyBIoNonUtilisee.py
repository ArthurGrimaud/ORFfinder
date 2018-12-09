
# def enter_seq():
#     adn = input("Entrez␣la␣chaine␣:␣")
#     for i in range (len(adn) ):
#         print("i:␣", i, "->" , adn[i])

# def Compte_aa(seq_prot):
#     """compte le nombre d'apparition d'un AA donné dans la séquence d'une prot"""
#     print("quel AA?")
#     aa = input()
#     nb_aa = 0
#     aa = aa.upper()
#     for i in seq_prot:
#         if i == aa :
#             nb_aa = nb_aa + 1
#     return nb_aa
#     print ("il␣y␣a␣" , Compte_aa(protein) , "␣Cysteine(s)")

# def Compte_all_aa(seq_prot):
#     "affiche le nombre d'apparition de chaque AA dans une seq prot"
#     List_aa =["E","D","A","R","N","C","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
#     for i in List_aa :
#         Nb_aa = 0
#         for j in seq_prot:
#             if i == j:
#                 Nb_aa = Nb_aa +1
#         print(i ,"->", Nb_aa)


# def Is_dna(adn):
#     """retourne tru si la sequence entree est une sequence d adn """
#     seq = ""
#     for i in adn:
#         i = i.upper()
#         seq = seq + i
#     flag = 0
#     for i in seq :
#         if i=="A" or i=="T" or i=="G" or i=="C":
#             flag = flag
#         else :
#             flag = flag +1
#     if flag > 0 :
#         return False
#     else :
#         return True

# def Count_word(seq, word):
#     """ seq est une string entrée, word est une string qui est recherchéen dans
#     la seq
#     Count_word retourne le nombre d'occurence de word dans la séquence et dans
#     tous les cadres de lecture"""
#     word = word.upper()
#     wlen = len(word)
#     flag = 0
#     for i in range(len(seq)):
#         Codon = One_word(seq,i,wlen)
#         if Codon == word :
#             flag = flag +1
#     print(word, "est présent", flag, "fois dans la séquence")
