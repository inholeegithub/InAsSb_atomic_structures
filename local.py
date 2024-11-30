#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tensorflow as tf
import warnings
import numpy as np
import os
import sys
from scipy import stats
from pymatgen.core import Composition, Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#from pyxtal import pyxtal
from random import choice
from scipy.stats import qmc
import time
import gc
from m3gnet.models import Relaxer
for category in (UserWarning, DeprecationWarning):
    warnings.filterwarnings("ignore", category=category, module="tensorflow")
tf.get_logger().setLevel('ERROR')



def get_lvectors(a, b, c, alpha, beta, gamma):
    abcabc = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    s = str(abcabc)
    a1 = np.zeros(3)
    a2 = np.zeros(3)
    a3 = np.zeros(3)
    a1[0] = float(s.split()[0])
    a1[1] = float(s.split()[1])
    a1[2] = float(s.split()[2])
    a2[0] = float(s.split()[3])
    a2[1] = float(s.split()[4])
    a2[2] = float(s.split()[5])
    a3[0] = float(s.split()[6])
    a3[1] = float(s.split()[7])
    a3[2] = float(s.split()[8])
    del abcabc, s
    return a1, a2, a3


def get_abcabc(a1, a2, a3):
    a = np.linalg.norm(a1)
    b = np.linalg.norm(a2)
    c = np.linalg.norm(a3)
    alpha = np.arccos(np.dot(a2, a3)/b/c) * 180./np.pi
    beta = np.arccos(np.dot(a3, a1)/c/a) * 180./np.pi
    gamma = np.arccos(np.dot(a1, a2)/a/b) * 180./np.pi
    return a, b, c, alpha, beta, gamma



def check_int(s):
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()


# In[2]:
def read_poscar0():
    scale0 = 1.
    a1 = np.zeros(3)
    a2 = np.zeros(3)
    a3 = np.zeros(3)
    jline = 0
    j0 = 0
    for line in sys.stdin:
        jline = jline+1
        if jline == 1:
            ncal = 0
            if check_int(line.split()[0]) :
                ncal = int(line.split()[0])
            continue
        if jline == 2:
            scale0 = float(line.split()[0])
            continue
        if jline == 3:
            a1[0] = float(line.split()[0])*scale0
            a1[1] = float(line.split()[1])*scale0
            a1[2] = float(line.split()[2])*scale0
        if jline == 4:
            a2[0] = float(line.split()[0])*scale0
            a2[1] = float(line.split()[1])*scale0
            a2[2] = float(line.split()[2])*scale0
        if jline == 5:
            a3[0] = float(line.split()[0])*scale0
            a3[1] = float(line.split()[1])*scale0
            a3[2] = float(line.split()[2])*scale0
        if jline == 6:
            nspecies = len(line.split())
            slist = []
            for i0 in range(nspecies):
                slist.append(line.split()[i0])
        if jline == 7:
            nspecies = len(line.split())
            nlist = []
            for i0 in range(nspecies):
                nlist.append(int(line.split()[i0]))
            natot = 0
            for i0 in range(nspecies):
                natot = natot+nlist[i0]
            drt = np.zeros((natot, 3))
        if jline == 8:
            lcart = False
            cstring = line.split()[0]
            if cstring[0] == 'S' or cstring[0] == 's':
                jline = jline - 1
            if cstring[0] == 'C' or cstring[0] == 'c':
                lcart = True
            continue
        if jline > 8:
            if lcart:
                d1 = float(line.split()[0])
                d2 = float(line.split()[1])
                d3 = float(line.split()[2])
                drt[j0, 0] = d1*scale0
                drt[j0, 1] = d2*scale0
                drt[j0, 2] = d3*scale0
                drt[j0, :] = cart2drtunit(a1, a2, a3, drt[j0, :])
                drt[j0, 0] = drt[j0, 0]-round(drt[j0, 0])
                drt[j0, 1] = drt[j0, 1]-round(drt[j0, 1])
                drt[j0, 2] = drt[j0, 2]-round(drt[j0, 2])
                if drt[j0, 0] < 0.:
                    drt[j0, 0] = drt[j0, 0]+1.
                if drt[j0, 1] < 0.:
                    drt[j0, 1] = drt[j0, 1]+1.
                if drt[j0, 2] < 0.:
                    drt[j0, 2] = drt[j0, 2]+1.
            else:
                d1 = float(line.split()[0])
                d2 = float(line.split()[1])
                d3 = float(line.split()[2])
                d1 = d1-round(d1)
                d2 = d2-round(d2)
                d3 = d3-round(d3)
                if d1 < 0.:
                    d1 = d1 + 1.
                if d2 < 0.:
                    d2 = d2 + 1.
                if d3 < 0.:
                    d3 = d3 + 1.
                drt[j0, 0] = d1
                drt[j0, 1] = d2
                drt[j0, 2] = d3
            j0 = j0+1
    sys.stdin.close()
    del j0, jline, scale0, d1, d2, d3, line, natot, lcart
    return a1, a2, a3, drt, slist, nlist, ncal



def write_poscar(ncal, a, b, c, alpha, beta, gamma, slist, nlist, f, e, ypath):
    a1, a2, a3 = get_lvectors(a, b, c, alpha, beta, gamma)
    file1 = open(ypath, 'w')
    s = str(ncal) + " " + str(e) + "\n" + "1." + "\n"
    s = s + str(a1[0]) + " " + str(a1[1]) + " " + str(a1[2]) + "\n"  \
          + str(a2[0]) + " " + str(a2[1]) + " " + str(a2[2]) + "\n"   \
          + str(a3[0]) + " " + str(a3[1]) + " " + str(a3[2]) + "\n"
    s = s + " ".join(slist) + "\n" + \
        " ".join(map(str, nlist)) + "\n" + "direct" + "\n"
    file1.write(s)
    nspecies = len(slist)
    k = 0
    for i in range(nspecies):
        for j in range(nlist[i]):
            file1.write(str(f[k, 0]) + " " +
                        str(f[k, 1]) + " " + str(f[k, 2]) + "\n")
            k = k + 1
    file1.close()
    del a1, a2, a3, i, j, k, s, nspecies


def write_poscar0(ncal, a, b, c, alpha, beta, gamma, slist, nlist, f, e ):
    a1, a2, a3 = get_lvectors(a, b, c, alpha, beta, gamma)
    s = str(ncal) + " " + str(e) + "\n" + "1." + "\n"
    s = s + str(a1[0]) + " " + str(a1[1]) + " " + str(a1[2]) + "\n"  \
          + str(a2[0]) + " " + str(a2[1]) + " " + str(a2[2]) + "\n"   \
          + str(a3[0]) + " " + str(a3[1]) + " " + str(a3[2]) + "\n"
    s = s + " ".join(slist) + "\n" + \
        " ".join(map(str, nlist)) + "\n" + "direct" 
    print(s)
    nspecies = len(slist)
    k = 0
    for i in range(nspecies):
        for j in range(nlist[i]):
            s1=str(f[k, 0]) + " " + str(f[k, 1]) + " " + str(f[k, 2])
            print(s1)
            k = k + 1




def get_local_poscar(cpath0):
    lcal = 0
    a1, a2, a3, drt, slist, nlist, ncal = read_poscar0()
    if True:
        nspecies = len(slist)
        natot = 0
        for i in range(nspecies):
            for j in range(nlist[i]):
               natot = natot + 1
        sym1 = []
        for i in range(nspecies):
            for j in range(nlist[i]):
                sym1.append(slist[i])
    if False:
        factor = 0.6
        ascale = 0.
        for i in range(natot):
            tmp = atomic_radii[sym1[i]]
            ascale = ascale - 4.*np.pi/3.*(tmp**3)
        if sym1[0] == 'Hg' or sym1[0] == 'Cs':
            ascale = ascale * 0.1
            factor = 0.1
        else:
            ascale = ascale * 2.4
        drt, kode = get_trial(slist, nlist, drt, a1, a2, a3, factor)
        if kode == 1:
            a, b, c, alpha, beta, gamma, drt = get_random(ascale, slist, nlist)
            a1, a2, a3 = get_lvectors(a, b, c, alpha, beta, gamma)
    aa, bb, cc, alpha, beta, gamma = get_abcabc(a1, a2, a3)
    nspecies = len(slist)
    sym1 = []
    for i in range(nspecies):
        for j in range(nlist[i]):
            sym1.append(slist[i])
    mo = Structure(Lattice.from_parameters(
        aa, bb, cc, alpha, beta, gamma), sym1, drt)
    lines_to_append = []
    start = time.time()
    relaxer = Relaxer(relax_cell=False)
    relax_results = relaxer.relax(mo, steps=600, fmax=0.09,  verbose=False)
    final_structure = relax_results['final_structure']
    final_energy = float(relax_results['trajectory'].energies[-1])
#
    loc1 = np.zeros((natot,3))
    for j in range(natot):
        loc1[j, 0] = final_structure.frac_coords[j, 0]
        loc1[j, 1] = final_structure.frac_coords[j, 1]
        loc1[j, 2] = final_structure.frac_coords[j, 2]
    abcabc = final_structure.lattice
    s = str(abcabc)
    a1 = np.zeros(3)
    a2 = np.zeros(3)
    a3 = np.zeros(3)
    a1[0] = float(s.split()[0])
    a1[1] = float(s.split()[1])
    a1[2] = float(s.split()[2])
    a2[0] = float(s.split()[3])
    a2[1] = float(s.split()[4])
    a2[2] = float(s.split()[5])
    a3[0] = float(s.split()[6])
    a3[1] = float(s.split()[7])
    a3[2] = float(s.split()[8])
    a, b, c, alpha, beta, gamma = get_abcabc(a1, a2, a3)
    if False:
        llattice = check_lattice(a, b, c, alpha, beta, gamma, natot, loc1, a1, a2, a3)
        if not llattice:
            final_energy = 1e99
#
        lines_to_append.append('fixed cell '+str(final_energy) +
                           ' ' + str(time.time()-start) + ' s')
    if False:
        print('fixed cell', final_energy, time.time()-start, 's')
    start1 = time.time()
    relaxer = Relaxer(relax_cell=True)
    relax_results = relaxer.relax(
        final_structure, steps=600, fmax=0.09, verbose=False)
    final_structure = relax_results['final_structure']
    final_energy = float(relax_results['trajectory'].energies[-1])
#
    loc1 = np.zeros((natot,3))
    for j in range(natot):
        loc1[j, 0] = final_structure.frac_coords[j, 0]
        loc1[j, 1] = final_structure.frac_coords[j, 1]
        loc1[j, 2] = final_structure.frac_coords[j, 2]
    abcabc = final_structure.lattice
    s = str(abcabc)
    a1 = np.zeros(3)
    a2 = np.zeros(3)
    a3 = np.zeros(3)
    a1[0] = float(s.split()[0])
    a1[1] = float(s.split()[1])
    a1[2] = float(s.split()[2])
    a2[0] = float(s.split()[3])
    a2[1] = float(s.split()[4])
    a2[2] = float(s.split()[5])
    a3[0] = float(s.split()[6])
    a3[1] = float(s.split()[7])
    a3[2] = float(s.split()[8])
    a, b, c, alpha, beta, gamma = get_abcabc(a1, a2, a3)
    if False:
        llattice = check_lattice(a, b, c, alpha, beta, gamma, natot, loc1, a1, a2, a3)
        if not llattice:
            final_energy = 1e99
#
        lines_to_append.append('variable cell '+str(final_energy)+' ' +
                               str(time.time()-start1) + ' s'+' '+str(time.time()-start)+' s')
    if False:
        print('variable cell', final_energy, time.time() -
              start1, 's', time.time()-start, 's')
    symprec = 0.1
#    sga2 = SpacegroupAnalyzer(final_structure)
#    ispg = sga2.get_space_group_number()
#    sga2.get_space_group_symbol()     sga2.volume   sga2.get_space_group_info()   supercell = polonium * (2, 2, 2) get_primitive_structure()
#    print( (time.time()-start), 's', ) # ispg  , sga2.get_space_group_symbol())
    abcabc = final_structure.lattice
    s = str(abcabc)
    a1 = np.zeros(3)
    a2 = np.zeros(3)
    a3 = np.zeros(3)
    a1[0] = float(s.split()[0])
    a1[1] = float(s.split()[1])
    a1[2] = float(s.split()[2])
    a2[0] = float(s.split()[3])
    a2[1] = float(s.split()[4])
    a2[2] = float(s.split()[5])
    a3[0] = float(s.split()[6])
    a3[1] = float(s.split()[7])
    a3[2] = float(s.split()[8])
    a, b, c, alpha, beta, gamma = get_abcabc(a1, a2, a3)
    if False:
        append_multiple_lines(cpath0, lines_to_append)
    del lines_to_append
    gc.collect()
    nspecies = len(slist)
    natot = 0
    for i in range(nspecies):
        for j in range(nlist[i]):
           natot = natot + 1
    loc1_save = np.zeros((natot,3))
    a_save = a
    b_save = b
    c_save = c
    alpha_save = alpha
    beta_save = beta
    gamma_save = gamma
    e_save = final_energy
    loc1_save[:, :] = final_structure.frac_coords[:, :]
    if lcal == 0:
#       write_poscar0(ncal, a, b, c, alpha, beta, gamma, slist, nlist, final_structure.frac_coords, final_energy)
        eref = final_energy
    for iitt in range(5):
        lcal = lcal + 1
#       a1, a2, a3, drt, slist, nlist, ncal = read_poscar0()
        aa, bb, cc, alpha, beta, gamma = get_abcabc(a1, a2, a3)
#       nspecies = len(slist)
#       natot = 0
#       for i in range(nspecies):
#           for j in range(nlist[i]):
#              natot = natot + 1
        loc1 = np.zeros((natot,3))
        for k in range(natot):
            loc1[k, 0] = final_structure.frac_coords[k, 0]
            loc1[k, 1] = final_structure.frac_coords[k, 1]
            loc1[k, 2] = final_structure.frac_coords[k, 2]
            drt[k, 0]=loc1[k, 0]
            drt[k, 1]=loc1[k, 1]
            drt[k, 2]=loc1[k, 2]
        klist = np.random.permutation(natot)
        sym1 = []
        k = 0
        for i in range(nspecies):
            for j in range(nlist[i]):
                sym1.append(slist[i])
                k1 = klist[k]
                loc1[k, 0]=drt[k1, 0] 
                loc1[k, 1]=drt[k1, 1] 
                loc1[k, 2]=drt[k1, 2] 
                k = k +1
        k = 0
        for i in range(nspecies):
            tmp = np.random.random()*np.pi
            tmq = np.random.random()/10.
            for j in range(nlist[i]):
                loc1[k, 0]=loc1[k, 0] + tmq*np.sin(loc1[k, 0]*np.pi-tmp)
                loc1[k, 1]=loc1[k, 1] + tmq*np.sin(loc1[k, 1]*np.pi-tmp)
                loc1[k, 2]=loc1[k, 2] + tmq*np.sin(loc1[k, 2]*np.pi-tmp)
                k = k +1
        mo = Structure(Lattice.from_parameters(
            aa, bb, cc, alpha, beta, gamma), sym1, loc1)
        lines_to_append = []
        start = time.time()
        relaxer = Relaxer(relax_cell=False)
        relax_results = relaxer.relax(mo, steps=600, fmax=0.09,  verbose=False)
        final_structure = relax_results['final_structure']
        final_energy = float(relax_results['trajectory'].energies[-1])
        loc1 = np.zeros((natot,3))
        for j in range(natot):
            loc1[j, 0] = final_structure.frac_coords[j, 0]
            loc1[j, 1] = final_structure.frac_coords[j, 1]
            loc1[j, 2] = final_structure.frac_coords[j, 2]
        abcabc = final_structure.lattice
        s = str(abcabc)
        a1 = np.zeros(3)
        a2 = np.zeros(3)
        a3 = np.zeros(3)
        a1[0] = float(s.split()[0])
        a1[1] = float(s.split()[1])
        a1[2] = float(s.split()[2])
        a2[0] = float(s.split()[3])
        a2[1] = float(s.split()[4])
        a2[2] = float(s.split()[5])
        a3[0] = float(s.split()[6])
        a3[1] = float(s.split()[7])
        a3[2] = float(s.split()[8])
        a, b, c, alpha, beta, gamma = get_abcabc(a1, a2, a3)
        if False:
            llattice = check_lattice(a, b, c, alpha, beta, gamma, natot, loc1, a1, a2, a3)
            if not llattice:
                final_energy = 1e99
#
            lines_to_append.append('fixed cell '+str(final_energy) +
                                   ' ' + str(time.time()-start) + ' s')
        if False:
            print('fixed cell', final_energy, time.time()-start, 's')
        start1 = time.time()
        relaxer = Relaxer(relax_cell=True)
        relax_results = relaxer.relax(
            final_structure, steps=600, fmax=0.09, verbose=False)
        final_structure = relax_results['final_structure']
        final_energy = float(relax_results['trajectory'].energies[-1])
#
        loc1 = np.zeros((natot,3))
        for j in range(natot):
            loc1[j, 0] = final_structure.frac_coords[j, 0]
            loc1[j, 1] = final_structure.frac_coords[j, 1]
            loc1[j, 2] = final_structure.frac_coords[j, 2]
        abcabc = final_structure.lattice
        s = str(abcabc)
        a1 = np.zeros(3)
        a2 = np.zeros(3)
        a3 = np.zeros(3)
        a1[0] = float(s.split()[0])
        a1[1] = float(s.split()[1])
        a1[2] = float(s.split()[2])
        a2[0] = float(s.split()[3])
        a2[1] = float(s.split()[4])
        a2[2] = float(s.split()[5])
        a3[0] = float(s.split()[6])
        a3[1] = float(s.split()[7])
        a3[2] = float(s.split()[8])
        a, b, c, alpha, beta, gamma = get_abcabc(a1, a2, a3)
        if False:
            llattice = check_lattice(a, b, c, alpha, beta, gamma, natot, loc1, a1, a2, a3)
            if not llattice:
                final_energy = 1e99
#
            lines_to_append.append('variable cell '+str(final_energy)+' ' +
                                   str(time.time()-start1) + ' s'+' '+str(time.time()-start)+' s')
        if False:
            print('variable cell', final_energy, time.time() -
                  start1, 's', time.time()-start, 's')
        symprec = 0.1
#        sga2 = SpacegroupAnalyzer(final_structure)
#        ispg = sga2.get_space_group_number()
#        sga2.get_space_group_symbol()     sga2.volume   sga2.get_space_group_info()   supercell = polonium * (2, 2, 2) get_primitive_structure()
#        print( (time.time()-start), 's', ) # ispg  , sga2.get_space_group_symbol())
        abcabc = final_structure.lattice
        s = str(abcabc)
        a1 = np.zeros(3)
        a2 = np.zeros(3)
        a3 = np.zeros(3)
        a1[0] = float(s.split()[0])
        a1[1] = float(s.split()[1])
        a1[2] = float(s.split()[2])
        a2[0] = float(s.split()[3])
        a2[1] = float(s.split()[4])
        a2[2] = float(s.split()[5])
        a3[0] = float(s.split()[6])
        a3[1] = float(s.split()[7])
        a3[2] = float(s.split()[8])
        a, b, c, alpha, beta, gamma = get_abcabc(a1, a2, a3)
        if False:
            append_multiple_lines(cpath0, lines_to_append)
        del lines_to_append
        gc.collect()
        if eref > final_energy :
            eref = final_energy
#           write_poscar0(ncal, a, b, c, alpha, beta, gamma, slist, nlist, final_structure.frac_coords, final_energy)
            a_save = a
            b_save = b
            c_save = c
            alpha_save = alpha
            beta_save = beta
            gamma_save = gamma
            e_save = final_energy
            loc1_save[:, :] = final_structure.frac_coords[:, :]
    write_poscar0(ncal, a_save, b_save, c_save, alpha_save, beta_save, gamma_save, slist, nlist, loc1_save, e_save)
# In[7]:

cpath0='output'
get_local_poscar(cpath0)
