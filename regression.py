import math
import random
import matplotlib.pyplot as plt

def Error(time,value,FORM,test,LEN) : #calcul l'erreur à minimisé
    J=0
    i=0
    while i<len(time) : #itère sur les valeur expérimentales
        if FORM =="exp" : #pour un modèle exponentiel
            J=pow((value[i]-value[0]*math.exp(test[0]*time[i])),2)+J #calcul de l'erreur
        if FORM =="linear" : #modèle linéaire
            J=pow((value[i]-(value[0]+test[0]*time[i])),2)+J
        if FORM =="poly" :
            F=value[0]
            f=0
            while(f<LEN) :
                F=F+test[f]*pow(time[i],f+1)
                f=f+1
            J=pow(value[i]-F,2)+J
        i=i+1
    return J

def Regression(time,value,FORM,step,LEN=1) : #la regression linéaire elle même
    totest=[] # list des combinaison de valeur à tester
    test=[] #la première combinaison de valeur que l'on initialise à 1
    tomin=[] #la liste contenant l'erreur de chaque valeur
    i=0
    while(i<LEN) : #initialisation de la première valeur
        test.append(1)
        i=i+1
    tomin.append(Error(time,value,FORM,test,LEN)) #on ajoute la valeur de l'erreur de test
    lencomb=1
    for x in test :
        lencomb=lencomb*3
    comb=[]
    for x in test : #on met en place une liste de longueur LEN qui permette de représenté toute les combinaison possible, on l'initialise à 
        comb.append(0) # 0 correspon à une valeur non modifier
        i=i+1
    A=True
    while(A) :
        totest.append(test)
        i=1
        while i<lencomb+1 : #on test toute les combinaison possible
            o=0
            while o<LEN :   #la combinaison en cours
                if (comb[o]+1)<3 :
                    comb[o]=comb[o]+1
                    break
                else :
                    comb[o]=0; 
                    o=o+1
            totest.append([])
            o=0
            while o<LEN : # on créer les valeur correspondant au combinaison
                if comb[o]==0 :
                    totest[i].append(test[o])
                if comb[o]==1 :
                    totest[i].append(test[o]+step)
                if comb[o]==2 :
                    totest[i].append(test[o]-step)                                        
                o=o+1
            tomin.append(Error(time,value,FORM,totest[i],LEN)) #on ajoute l'erreur à la liste des erreurs
            i=i+1
        mini=float('inf')
        for x in tomin :
            if x<mini :
                mini=x
        index=tomin.index(mini)
        test=totest[index]
        tomin=[]
        totest=[]
        tomin.append(mini)
        if index==0 :
            return test,tomin[0]

def plot(time,value,test,FORM,LEN=1) :
    a=[]
    for t in time :
        if FORM =="exp" :
            a.append(value[0]*math.exp(test[0]*t))
        if FORM =="linear" :
            a.append(value[0]+test[0]*t)
        if FORM =="poly" :
            F=value[0]
            f=0
            while(f<LEN) :
                F=F+test[f]*pow(t,f+1)
                f=f+1
            a.append(F)            
    plt.plot(time,value)
    plt.plot(time,a)
    plt.show()

#pharmacocinétique
#time = [0,0.5,1,1.5,2,3,5,7.5,10,15,20,24]
#value = [23.8,22.1,20.5,19,17.6,15.2,11.3,7.82,5.39,2.57,1.22,0.68]

# population italienne covid
time=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43]
value=[2313,2653,2548,3499,3593,3235,3527,4208,5324,5988,6554,5560,4783,5240,5198,6202,5907,5973,5215,4047,4053,4783,4669,4585,4807,4318,3599,3037,3834,4204,3950,4697,4094,3153,2973,2666,3787,3494,3492,3047,2251,2727,3370,2644]

test,epsi = Regression(time,value,"poly",0.001,4)
print("vecteur de a :")
print(test)
print("valeur de J :")
print(epsi)
norm=0
for a in test : #calcul de l'erreur relative
    norm=pow(a,2)+norm  #somme du carré des coordonné
epsi=math.sqrt(epsi)/math.sqrt(norm) #on divise la racine de j par la norme, qui est la racine de la somme des carré
print("valeur d'epsilon :")
print(epsi)
plot(time,value,test,"poly",4)
