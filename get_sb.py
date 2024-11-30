import numpy as np
afile=open('POSCAR-conv','r')
jline=0
for line in afile:
    jline=jline+1
    if jline > 1 :
        if jline == 2 :
            ascale=float(line.split()[0])
            a1=np.zeros(3)
            a2=np.zeros(3)
            a3=np.zeros(3)
        if jline == 3 :
            a1[0]=float(line.split()[0])
            a1[1]=float(line.split()[1])
            a1[2]=float(line.split()[2])
        if jline == 4 :
            a2[0]=float(line.split()[0])
            a2[1]=float(line.split()[1])
            a2[2]=float(line.split()[2])
        if jline == 5 :
            a3[0]=float(line.split()[0])
            a3[1]=float(line.split()[1])
            a3[2]=float(line.split()[2])
        if jline == 6 :
            nspecies=len(line.split())
        if jline == 7 :
            nspecies=len(line.split())
#           print(int(line.split()[0]),int(line.split()[1]))
            natot=0
            for k in range(len(line.split())):
                natot=natot+int(line.split()[k])
#           print(natot)
            drt=np.zeros((natot,3))
            kount=0
        if jline == 8 :
#           print(line.split()[0])
            k1= len(line.split()[0])
        if jline >  8 :
            drt[kount,0]=float(line.split()[0])
            drt[kount,1]=float(line.split()[1])
            drt[kount,2]=float(line.split()[2])
            kount=kount+1
afile.close()

slist=["In", "As", "Sb"]
nspecies=len(slist)
itmp=int(7./93.*108)
#print(itmp)
nlist=[108, 108-itmp, itmp]
#print(slist)
#print(nlist)
for ii in range(1000):
    alist=[k for k in range(108)]
    blist=[k for k in range(108, natot)]
    brr = np.array(blist)
    np.random.shuffle(brr)
    blist=[brr[k] for k in range(108)]
    clist= alist+blist
#   print(clist)
    ypath = 'POSCAR_'+str(ii).zfill(4)
    file1 = open(ypath, 'w')
    s = str(ii) + " " + "\n" + "1." + "\n"
    s = s + str(a1[0]) + " " + str(a1[1]) + " " + str(a1[2]) + "\n"  \
          + str(a2[0]) + " " + str(a2[1]) + " " + str(a2[2]) + "\n"   \
          + str(a3[0]) + " " + str(a3[1]) + " " + str(a3[2]) + "\n"
    s = s + " ".join(slist) + "\n" + \
        " ".join(map(str, nlist)) + "\n" + "direct" + "\n"
    file1.write(s)
    nspecies = len(slist)
    k1 = 0
    for i in range(nspecies):
        for j in range(nlist[i]):
            k = clist[k1]
            file1.write(str(drt[k, 0]) + " " +
                        str(drt[k, 1]) + " " + str(drt[k, 2]) + "\n")
            k1 = k1 + 1
    file1.close()
