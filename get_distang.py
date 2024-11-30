import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


alist0 = []
afile=open('POSCAR_out_InAs0','r')
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
vol0 = np.dot(a3, np.cross(a1,a2))
if False:
    print( vol0, vol0/108.)

slist=["In", "As" ]
nlist=[108, 108 ]
nspecies=len(slist)
sym = []
for i in range(nspecies):
    for j in range(nlist[i]):
        sym.append(slist[i])
if False:
    print(len(sym))

for i in range(natot):
   for j in range(natot):
       if sym[i] == 'As' and sym[j] == 'In' :
           d1=drt[i,0]-drt[j,0]
           d2=drt[i,1]-drt[j,1]
           d3=drt[i,2]-drt[j,2]
           d1=d1-np.rint(d1)
           d2=d2-np.rint(d2)
           d3=d3-np.rint(d3)
           dx=d1*a1[0]+d2*a2[0]+d3*a3[0]
           dy=d1*a1[1]+d2*a2[1]+d3*a3[1]
           dz=d1*a1[2]+d2*a2[2]+d3*a3[2]
           r=np.sqrt(dx*dx+dy*dy+dz*dz)
#          print(r)
           if r < 3.:
               alist0.append(r)
#print(alist0)
print(alist0[0])
if True:
    clist0 = []
    for i in range(natot):
       if sym[i] == 'As' :
           for j in range(natot):
               if sym[j] == 'In' :
                   d1=drt[i,0]-drt[j,0]
                   d2=drt[i,1]-drt[j,1]
                   d3=drt[i,2]-drt[j,2]
                   d1=d1-np.rint(d1)
                   d2=d2-np.rint(d2)
                   d3=d3-np.rint(d3)
                   dx1=d1*a1[0]+d2*a2[0]+d3*a3[0]
                   dy1=d1*a1[1]+d2*a2[1]+d3*a3[1]
                   dz1=d1*a1[2]+d2*a2[2]+d3*a3[2]
                   rij=np.sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
                   if rij < 3.0 :
                       for k in range(natot):
                           if sym[k] == 'In' and k > j:
                               d1=drt[i,0]-drt[k,0]
                               d2=drt[i,1]-drt[k,1]
                               d3=drt[i,2]-drt[k,2]
                               d1=d1-np.rint(d1)
                               d2=d2-np.rint(d2)
                               d3=d3-np.rint(d3)
                               dx2=d1*a1[0]+d2*a2[0]+d3*a3[0]
                               dy2=d1*a1[1]+d2*a2[1]+d3*a3[1]
                               dz2=d1*a1[2]+d2*a2[2]+d3*a3[2]
                               rik=np.sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
                               if rik < 3.0 :
                                   angle=np.arccos((dx1*dx2+dy1*dy2+dz1*dz2)/(rij*rik))
                                   angle=angle*180./np.pi
                                   clist0.append(angle)
print(clist0[0])
alist = []
blist = []
dlist = []
elist = []
for ii in range(1000):
    fname='POSCAR_out_'+str(ii).zfill(4)
    afile=open(fname,'r')
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
#               print(int(line.split()[0]),int(line.split()[1]))
                natot=0
                for k in range(len(line.split())):
                    natot=natot+int(line.split()[k])
#               print(natot)
                drt=np.zeros((natot,3))
                kount=0
            if jline == 8 :
#               print(line.split()[0])
                k1= len(line.split()[0])
            if jline >  8 :
                drt[kount,0]=float(line.split()[0])
                drt[kount,1]=float(line.split()[1])
                drt[kount,2]=float(line.split()[2])
                kount=kount+1
    afile.close()
    vol = np.dot(a3, np.cross(a1,a2))
#   print( vol, vol/108.)
    if vol > vol0  and vol/vol0 > 1.05:
        continue
    else:
        print(vol,vol0, vol-vol0, vol/vol0)

    slist=["In", "As", "Sb"]
    itmp=int(7./93.*108)
    nlist=[108, 108-itmp, itmp]
    nspecies=len(slist)
    sym = []
    for i in range(nspecies):
        for j in range(nlist[i]):
            sym.append(slist[i])
#   print(len(sym))

    for i in range(natot):
       for j in range(natot):
           if sym[i] == 'As' and sym[j] == 'In' :
               d1=drt[i,0]-drt[j,0]
               d2=drt[i,1]-drt[j,1]
               d3=drt[i,2]-drt[j,2]
               d1=d1-np.rint(d1)
               d2=d2-np.rint(d2)
               d3=d3-np.rint(d3)
               dx=d1*a1[0]+d2*a2[0]+d3*a3[0]
               dy=d1*a1[1]+d2*a2[1]+d3*a3[1]
               dz=d1*a1[2]+d2*a2[2]+d3*a3[2]
               r=np.sqrt(dx*dx+dy*dy+dz*dz)
#              print(r)
               if r < 3.:
                   alist.append(r)
    for i in range(natot):
       for j in range(natot):
           if sym[i] == 'Sb' and sym[j] == 'In' :
               d1=drt[i,0]-drt[j,0]
               d2=drt[i,1]-drt[j,1]
               d3=drt[i,2]-drt[j,2]
               d1=d1-np.rint(d1)
               d2=d2-np.rint(d2)
               d3=d3-np.rint(d3)
               dx=d1*a1[0]+d2*a2[0]+d3*a3[0]
               dy=d1*a1[1]+d2*a2[1]+d3*a3[1]
               dz=d1*a1[2]+d2*a2[2]+d3*a3[2]
               r=np.sqrt(dx*dx+dy*dy+dz*dz)
#              print(r)
               if r < 3.:
                   blist.append(r)
    if True:
        for i in range(natot):
           if sym[i] == 'As' :
               for j in range(natot):
                   if sym[j] == 'In' :
                       d1=drt[i,0]-drt[j,0]
                       d2=drt[i,1]-drt[j,1]
                       d3=drt[i,2]-drt[j,2]
                       d1=d1-np.rint(d1)
                       d2=d2-np.rint(d2)
                       d3=d3-np.rint(d3)
                       dx1=d1*a1[0]+d2*a2[0]+d3*a3[0]
                       dy1=d1*a1[1]+d2*a2[1]+d3*a3[1]
                       dz1=d1*a1[2]+d2*a2[2]+d3*a3[2]
                       rij=np.sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
                       if rij < 3.0 :
                           for k in range(natot):
                               if sym[k] == 'In' and k > j:
                                   d1=drt[i,0]-drt[k,0]
                                   d2=drt[i,1]-drt[k,1]
                                   d3=drt[i,2]-drt[k,2]
                                   d1=d1-np.rint(d1)
                                   d2=d2-np.rint(d2)
                                   d3=d3-np.rint(d3)
                                   dx2=d1*a1[0]+d2*a2[0]+d3*a3[0]
                                   dy2=d1*a1[1]+d2*a2[1]+d3*a3[1]
                                   dz2=d1*a1[2]+d2*a2[2]+d3*a3[2]
                                   rik=np.sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
                                   if rik < 3.0 :
                                       angle=np.arccos((dx1*dx2+dy1*dy2+dz1*dz2)/(rij*rik))
                                       angle=angle*180./np.pi
                                       dlist.append(angle)
    if True:
        for i in range(natot):
           if sym[i] == 'Sb' :
               for j in range(natot):
                   if sym[j] == 'In' :
                       d1=drt[i,0]-drt[j,0]
                       d2=drt[i,1]-drt[j,1]
                       d3=drt[i,2]-drt[j,2]
                       d1=d1-np.rint(d1)
                       d2=d2-np.rint(d2)
                       d3=d3-np.rint(d3)
                       dx1=d1*a1[0]+d2*a2[0]+d3*a3[0]
                       dy1=d1*a1[1]+d2*a2[1]+d3*a3[1]
                       dz1=d1*a1[2]+d2*a2[2]+d3*a3[2]
                       rij=np.sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
                       if rij < 3.0 :
                           for k in range(natot):
                               if sym[k] == 'In' and k > j:
                                   d1=drt[i,0]-drt[k,0]
                                   d2=drt[i,1]-drt[k,1]
                                   d3=drt[i,2]-drt[k,2]
                                   d1=d1-np.rint(d1)
                                   d2=d2-np.rint(d2)
                                   d3=d3-np.rint(d3)
                                   dx2=d1*a1[0]+d2*a2[0]+d3*a3[0]
                                   dy2=d1*a1[1]+d2*a2[1]+d3*a3[1]
                                   dz2=d1*a1[2]+d2*a2[2]+d3*a3[2]
                                   rik=np.sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
                                   if rik < 3.0 :
                                       angle=np.arccos((dx1*dx2+dy1*dy2+dz1*dz2)/(rij*rik))
                                       angle=angle*180./np.pi
                                       elist.append(angle)
alist0 = np.array(alist0)
factor=alist0[0]
alist = np.array(alist)
blist = np.array(blist)
dlist = np.array(dlist)
elist = np.array(elist)
alist = alist/factor
blist = blist/factor
if False:
    n_bins = 50
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    axs[0].hist(alist, bins=n_bins)
    axs[1].hist(blist,  bins=n_bins)
    plt.show()
    plt.savefig('dist_distri.pdf')
    plt.close()


if True:
    fig, ax1 = plt.subplots()
    ax1.hist([alist, blist], bins=50, alpha=0.5, stacked=False, color=['cyan', 'Purple'], edgecolor='skyblue', lw=1)
#   ax1.hist(alist, bins=50, alpha=0.5, color='cyan', edgecolor='skyblue', lw=1)
#   ax1.hist(blist, bins=50, alpha=0.5, color='Purple', edgecolor='k', lw=2)
#ax1.set_xlim(2.65, 2.95)
    ax1.set_xlim(0.97, 1.08)
    ax1.set_ylabel("Count (arb. unit)",fontsize=18)
    ax1.set_xlabel("Bond-length ratio",fontsize=18)
#ax1.axvline(x=2.706448, color='b', label='In-As')
#ax1.axvline(x=2.706448, color='b')
    ax1.axvline(x=1.000, color='magenta', linestyle='--', linewidth=1)
    ax1.xaxis.set_major_locator(MultipleLocator(0.02))
#ax1.xaxis.set_major_formatter('{x:.0f}')
    ax1.xaxis.set_minor_locator(MultipleLocator(0.005))
    plt.tight_layout()

# Adding labels and title
#plt.xlabel('Values')
#plt.ylabel('Frequency')
#plt.title('Stacked Histogram')

    plt.legend(['In-As(reference)', 'In-As', 'In-Sb'], fontsize=18)
    plt.savefig("bddist.pdf")
    plt.show()
    plt.close()
    print('alist')
    for i in range(len(alist)):
        print(alist[i])
    print('blist')
    for i in range(len(blist)):
        print(blist[i])


if True:
    fig, ax1 = plt.subplots()
    ax1.hist([dlist, elist], bins=50, alpha=0.5, stacked=False, color=['cyan', 'Purple'], edgecolor='skyblue', lw=1)
#   ax1.hist( dlist,  bins=50, alpha=0.5, color='cyan', edgecolor='skyblue', lw=1)
#   ax1.hist( elist,  bins=50, alpha=0.5, color='Purple', edgecolor='k', lw=2)

    ax1.set_xlim(102.00, 114.00)
    ax1.set_ylabel("Count (arb. unit)",fontsize=18)
    ax1.set_xlabel("Bond angle (degree)",fontsize=18)
    ax1.axvline(x=109.4712406, color='magenta', linestyle='--', linewidth=1)
    ax1.xaxis.set_major_locator(MultipleLocator(2.00))
#ax1.xaxis.set_major_formatter('{x:.0f}')
    ax1.xaxis.set_minor_locator(MultipleLocator(0.25))
    plt.tight_layout()
    plt.legend(['In-As-In(reference)', 'In-As-In', 'In-Sb-In'], fontsize=18, loc='best')
    plt.savefig("angdist.pdf")
    plt.show()
    plt.close()
    print('dlist')
    for i in range(len(dlist)):
        print(dlist[i])
    print('elist')
    for i in range(len(elist)):
        print(elist[i])
