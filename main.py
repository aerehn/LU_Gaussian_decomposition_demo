import time
# demo yhtälöryhmän ratkaisuajoista LU- ja Gaussin hajotelmia käyttämällä
#määritellään laskettavat matriisit
A = [[2, 3, 1],
     [4, 9, 6],
     [10, 24, 18]]
A2=[[9,3,7,5,8],
    [4,7,5,8,4],
    [43,6,3,67,8],
    [47,28,4,9,0],
    [67,38,2,4,87]]
A3=[[9,3,7,5,8,9,3,7,5,8],#10
    [4,7,5,8,4,23,4,5,67,4],#10
    [43,6,3,67,8,89,4,23,4,21],#10
    [47,28,4,9,0,45,32,45,32,54],#10,
    [67,38,2,4,87,34,32,5,333,323,],#10
    [34,4,5,3,4,5,3,4,232,23],#10
    [23,21,4,33,23,32,34,23,21,23],#10
    [3,2,5,3,4,6,3,4,3,67],#10
    [23,34,4,45,23,54,32,12,3,12],#10
    [84,238,456,234,12,653,2,343,234,32]]
B = [12, 39, 108]
B2=[34,2,4,5,6]
B3=[23,43,32,43,21,2,3,45,32,21]

B_vektorien_määrä=[50000, 100000, 150000, 200000]
 # ratkaisujen määrä

def print_2D(matrix):#voidaan tulostaa kaksiulotteinen matriisi
    for i in matrix:
        print(i)

def Gaussin_hajotelma(matrixA,matrixB):
    # Tehdään gaussin hajotelma
    n=len(matrixB)

    for k in range(n):
        for i in range(k + 1, n):
            r = matrixA[i][k] / matrixA[k][k]# lasketaan kerroin r

            for l in range(k, n):#tehdään tarvittavat vähennykset alemmista riveistä

                matrixA[i][l] = matrixA[i][l] - r * matrixA[k][l]#päivitetään A-matriisin arvot
            matrixB[i] = matrixB[i] - r * matrixB[k]#lasketaan uudet B-matriisin arvot

    # tehdään ratkaisulle oma vektorinsa
    x=[0 for i in range(n)]


    #yhdistetään A ja B takaisin sijoittamisen helpottamiseksi
    for i in range(n):
        matrixA[i].append(matrixB[i])


    #itse takaisin sijoitus
    for m in range(n-1,-1,-1):
        x[m] = matrixA[m][n]/matrixA[m][m]#lasketaan muutujan arvo
        for k in range(m-1,-1,-1):#vähennetään ylempien rivien viimeisistä sarakkeista niitä vastaava määrä viimeksi laskettua muuttujaa
            matrixA[k][n] -= matrixA[k][m]*x[m]
    #Esim. ensimmäisellä kierroksella lasketaan suoraan z:n arvo viimeiseltä riviltä 3/1=3. Sitten kahden ylemmän rivin viimeisestä arvosta eli muuttujien summista(B: arvot)
    #vähennetään sitä vastaavien rivien z:n kertoimien verran z:taa eli toisella rivillä 15-4*3=3 ja ekalla rivillä 12-3=9.
    return(x)


# LU-hajotelman funktio. Tuottaa vain L-matriisin
def L_hajotelma(matrixA):
    n=len(matrixA)
    #Alustetaan nxn kokoinen identiteettimatriisi L-matriisia varten
    L=[]
    for i in range(n):
        x=[0 for i in range(n)]# tehdään vaakarivi, jolla kaikki nollaa
        x[i]=1#lisätään x:n indeksiin i
        L.append(x)
    #Tehdään LU hajotelma U on tässä tapauksessa matrixA
    for k in range(n):
        for i in range(k + 1, n):
            r = matrixA[i][k] / matrixA[k][k]
            #Jaettava Amatrix:in arvo on samassa kohdassa kuin kertoimen tallennussijainti L:ssä
            #otetaan r:n arvot talteen
            L[i][k]=r
            for l in range(k, n):
                matrixA[i][l] = matrixA[i][l] - r * matrixA[k][l]
    return L

#tehdään sama temppu kuin L_hajotelma-metodissa, mutta nyt palautetaan muokattu matrixA
def U_hajotelma(matrixA):
    n = len(matrixA)
    for k in range(n):
        for i in range(k + 1, n):
            r = matrixA[i][k] / matrixA[k][k]
            for l in range(k, n):
                matrixA[i][l] = matrixA[i][l] - r * matrixA[k][l]
    return matrixA
# LU-hajotelma menetelmän takaisin sijoitus funktio
def LU_sijoitus(matrixL, matrixU, vectorB):
    n = len(vectorB)
    z = [0 for i in range(n)]
    x = [0 for i in range(n)]
    #eteenpäin sijoitus
    #lisätään B-vektorin arvot L-matriisin loppun
    for i in range(n):
        matrixL[i].append(vectorB[i])
    #print_2D(matrixL)
    for m in range(0,n):
        z[m] = matrixL[m][n]/matrixL[m][m]#lasketaan muutujan arvo
        for k in range(m,n):#vähennetään ylempien rivien viimeisistä sarakkeista niitä vastaava määrä viimeksi laskettua muuttujaa
            matrixL[k][n] -= matrixL[k][m]*z[m]
            #print()
            #print_2D(matrixL)

    # lisätään z-vektorin arvot U-matriisin loppun
    for i in range(n):
        matrixU[i].append(z[i])
    #takaisin sijoitus
    for m in range(n - 1, -1, -1):
        x[m] = matrixU[m][n] / matrixU[m][m]  # lasketaan muutujan arvo
        for k in range(m - 1, -1,-1):  # vähennetään ylempien rivien viimeisistä sarakkeista niitä vastaava määrä viimeksi laskettua muuttujaa
            matrixU[k][n] -= matrixU[k][m] * x[m]
    return x
    #takaisin sijoitus
#molemmat aikaajo-metodit toimivat käytännössä samalla tavalla. Ainoana erona on, että LU-hajotelmalla jokaisesta matriisia kohden hajotelmat lasketaan vain kerran.
#koska ei ole järkevää tehdä 200 000 erillistä b-vektoria niin lasketaan b-vektorille ratkaisu 200 000 kertaa.
def aja_LU_aika(MatrixA,vectorB,vektorien_määrä):
    ajat=[ 0 for i in range(len(vektorien_määrä)) ]
    #lastetaan ajat eri määrille b-vektoreita
    for i in range(len(vektorien_määrä)):
        if len(MatrixA) != len(vectorB):
            return 0
        start = time.time()#aloitetaan ajan mittaaminen
        U=U_hajotelma(MatrixA)#Lasketaan U-hajotelma
        L=L_hajotelma(MatrixA)#Lasketaan L-hajotelma
        ratkaisut = [0 for h in range(vektorien_määrä[i])]
        #Lasketaan ratkaisu b-vektorille vektorien_määrä indeksin i monta kertaa.
        for j in range(vektorien_määrä[i]):
            ratkaisut[j] = LU_sijoitus(L, U, vectorB)#ratkaisu haluttaisiin varmaankin myös tallentaa johonkin, jos sillä tehtäisiin myöhemmin jotain
        end=time.time()#katkaistaan aika
        ajat[i] = end-start#tallennetaan aika vektoriin ajat
    return ajat


def aja_Gaussin_aika(MatrixA,vectorB,vektorien_määrä):
    ajat = [0 for i in range(len(vektorien_määrä))]
    if len(MatrixA) != len(vectorB):
        return 0
    for i in range(len(vektorien_määrä)):
        start = time.time()#aloitetaan ajan mittaaminen
        ratkaisut = [0 for h in range(vektorien_määrä[i])]
        for j in range(vektorien_määrä[i]):
            ratkaisut[j] = Gaussin_hajotelma(MatrixA, vectorB)
        end=time.time()
        ajat[i] = end-start
    return ajat
#print("L-hajotelma:")
#L=L_hajotelma(A)
#print_2D(L)
#print("U-hajotelma:")
#U=U_hajotelma(A)
#print_2D(U)
#print("Vektori x:")
#print(LU_sijoitus(L, U, B))
#print()
#print("Gaussin hajotelmalla:")
#print(Gaussin_hajotelma(A,B))
start = time.time()
print("LU-hajotelman suoritusajat. 3x3,5x5,10x10 b-vektorien määrillä 50000, 100000, 150000 ja 200000")
print("kukin vaakarivi vastaa yhtä A-matriisin kokoa")
print(aja_LU_aika(A,B,B_vektorien_määrä))#matriisin koko:3x3
print(aja_LU_aika(A2,B2,B_vektorien_määrä))#matriisin koko:5x5
print(aja_LU_aika(A3,B3,B_vektorien_määrä))#matriisin koko:10x10
print()
print("Gaussin hajotelman ajat:")
print(aja_Gaussin_aika(A,B,B_vektorien_määrä))#matriisin koko:3x3
print(aja_Gaussin_aika(A2,B2,B_vektorien_määrä))#matriisin koko:5x5
print(aja_Gaussin_aika(A3,B3,B_vektorien_määrä))#matriisin koko:10x10
end = time.time()
print("ohjelmanajamiseen kätetty aika: "+str(end-start))