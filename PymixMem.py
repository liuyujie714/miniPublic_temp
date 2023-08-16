#!/usr/bin/python3
#coding: utf-8
# Release: 2023/08/16

import numpy as np
import MDAnalysis as mda
from MDAnalysis.transformations.rotate import rotateby
from MDAnalysis.transformations.translate import center_in_box
from MDAnalysis.lib.distances import calc_bonds
from itertools import chain
from random import shuffle
import argparse, sys

class PymixMem:
    def __init__(self, fp:str="settings.inp", fgro:str=None, dx:float=10.0, frac:float=2.2,
                 nlayer:int=2) -> None:
        self.fp = fp
        self.fgro = fgro
        self.dx = dx # the space for each lipids
        self.frac = frac # the distance scale factor between two monolayers
        if nlayer < 1 or nlayer > 2:
            raise ValueError('-n must be 1 or 2, represents monlayer and bilayer')
        self.nlayer = nlayer # the number of membrane layers, 1 or 2

    def gen_mem(self):
        """ generate a membrane structure
        """
        inp = self.read_inut()
        ntot = np.sum([inp[k][0] for k in inp.keys()])
        print(f'The number of lipids in monolayer: {ntot}')
        # [0, 0, 0, 0, 1, 1, 1, ... 2, 2, 2, ...]
        order = list(chain.from_iterable([[idx]*count for idx, (count, _) in enumerate(inp.values())]))
        pdbinfo = self.get_pdbinfo(inp)
        self.pack_mol(inp, order, pdbinfo)

    def get_pdbinfo(self, inp:dict):
        """ get coordinates of input lipids
        """
        pdbinfo = []
        for k in inp.keys():
            u = mda.Universe(k) # read lipid pdb
            # put long axis of molecule along Z axis
            n = u.atoms.n_atoms
            print(f'Reading {k}, total {n} atoms...')
            ia = ib = 0
            dist = -999999
            for i in range(n):
                for j in range(i+1, n):
                    r = calc_bonds(u.atoms.positions[i], u.atoms.positions[j],
                                   box=None)
                    if r > dist:
                        dist, ia, ib = r, i, j
            vec = u.atoms.positions[ia]-u.atoms.positions[ib]
            if (inp[k][1] != 0):
                vec = -vec
            vec /= np.linalg.norm(vec)
            zaxis = [0,0,1]
            angle = np.arccos(vec[2])*180.0/np.pi # simple dot
            raxis = np.cross(vec, zaxis)
            rotateby(angle, direction=raxis, ag=u.atoms)(u.trajectory.ts)
            # move molecules to zero
            center_in_box(ag=u.atoms, point=[0, 0, 0])(u.trajectory.ts)
            # u.atoms.write(k.rstrip(".pdb")+"_rotate.pdb")
            pdbinfo.append([u.atoms.resnames, u.atoms.names, u.atoms.positions])
        return pdbinfo

    def pack_mol(self, inp:dict, order:list, pdbinfo:list) -> None:
        """ pack lipid to build double layers membrane

            Parameters
            ----------
            inp: input settings
            order: a list of pdb file idex (0-based) with shuffle
            pdbinfo: a list contains [resnames, atomnames, positions]
        """
        # distribution number in each side in each layer
        ntot = len(order)
        nside = int(np.sqrt(ntot))
        nrest, nx = 0, 0
        while nside*nside < ntot:
            nside += 1
        # adjust the space of last column lipid
        if nside*nside - ntot > 0:
            nx = nside
            while nside*nx > ntot:
                nx -= 1
            nrest = ntot - nside*nx
        # each layer
        maxZ = -99999
        Natoms = 0
        nmol = [count for count, _ in inp.values()]
        for idx, (res, _, coords) in enumerate(pdbinfo):
            Natoms += len(res)*nmol[idx]
            val = np.max(coords)
            if val > maxZ:
                maxZ = val
        fout = open("mixmem.gro", "w")
        fmol = open("mixmem.txt", "w")
        fout.write(f'Membrane structure\n{self.nlayer*Natoms}\n')
        FMT_GRO = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"
        molID, totAtom = 1, 1
        grolines, uniq_res = [], {}
        for layer in range(self.nlayer): 
            dx = dy = self.dx
            # shuffle molecules order
            shuffle(order)
            # each molecule in this layer
            break_out = False
            cout = 0
            for i in range(nside):
                for j in range(nside):
                    if nrest > 0 and i==nx:
                        dy = self.dx*nside/nrest
                    if cout >= ntot:
                        break_out = True
                        break
                    id = order[cout] # molecule type index
                    natom = len(pdbinfo[id][0]) 
                    for k in range(natom):
                        res = pdbinfo[id][0][k]
                        name = pdbinfo[id][1][k]
                        # unit is Angstrom
                        x = pdbinfo[id][2][k][0] + i*dx
                        y = pdbinfo[id][2][k][1] + j*dy
                        z = pdbinfo[id][2][k][2]
                        # another layer
                        if layer == 1:
                            z = -z - maxZ*self.frac
                        # fout.write(FMT_GRO %(molID%100000, res, name, totAtom%100000, x/10., y/10., z/10.))
                        grolines.append(FMT_GRO %(molID%100000, res, name, totAtom%100000, x/10., y/10., z/10.))
                        totAtom += 1
                    molID+=1
                    cout += 1
                    # fmol.write('%6s %5d\n' %(pdbinfo[id][0][0], 1))
                    uniq_res[pdbinfo[id][0][0]] = natom
                if break_out:
                    break

        # sort residue order
        new_grolines = []
        skip_count = 0
        for res, nat in uniq_res.items():
            nmol_ = 0
            for idx, line in enumerate(grolines):
                if skip_count > 0:
                    skip_count -= 1
                    continue
                if res == line[5:10].strip():
                    new_grolines.append(grolines[idx:idx+nat])
                    skip_count = nat-1
                    nmol_ += 1
            fmol.write('%6s %5d\n' %(res, nmol_))
        for line in new_grolines:
            fout.writelines(line)
        fout.write("%10.5f %10.5f %10.5f\n" %(0.0, 0.0, 0.0))
        fout.close()
        fmol.close()
        print(f'\nFinal mix membrane has been saved in "mixmem.gro"')

    def read_inut(self) -> dict:
        """ read input parameters

            Returns
            -------
            a dict { 'pdb_file': (nmol, direct) }
        """
        inp = {}
        with open(self.fp, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.lstrip()
                if len(line) < 2 or line[0] in [';', '#']:
                    continue
                ll = line.split()
                if len(ll) == 3:
                    try:
                        fpdb, count, bInv = ll[0], int(ll[1]), int(ll[2])
                        inp[fpdb] = (count, bInv)
                    except:
                        raise SyntaxError(f'Can not recognize line: "{line.rstrip()}" in {self.fp}')
        return inp

    def remove_water(self, lower:float=-9999., upper:float=9999.):
        """ Remove water molecules that inside of membrane
        """
        water_res = "(resname SOL or resname TIP3 or resname HOH or resname ICE)"
        print(f'Remove water {lower} nm < z < {upper} nm')
        lower, upper = lower*10., upper*10. # nm to angstrom
        u = mda.Universe(self.fgro)
        sel = u.select_atoms(f'not (same resid as (prop z>{lower} and prop z<{upper} and {water_res} ))')
        sel.write('sol_2.gro')
        print(f'Write final structure to "sol_2.gro" with {sel.select_atoms(water_res).atoms.n_atoms} atoms of water molecules')

def parser_opt():
    parser = argparse.ArgumentParser(description='A tool for building mix monolayer and bilayer membrane')
    parser.add_argument('-f', '--file', metavar='settings.inp', type=str, default=None,
                         help='The input parameters file')
    parser.add_argument('-dx', '--dx', metavar=10.0, type=float, default=10.0, 
                        help='The space of molecules in side, default=10.0, unit is angstrom')
    parser.add_argument('-frac', '--frac', metavar=2.2, type=float, default=2.2,  
                        help='The scale factor of bilayer distance, default=2.2')
    parser.add_argument('-n', '--nlayer', metavar=2, type=int, default=2,  
                        help='The number of membrane layers, only be 1 or 2, default=2')
    parser.add_argument('-g', '--gro', metavar='sol.gro', type=str, default=None,
                         help='The input solvated box for removing water that inside of membrane')
    parser.add_argument('-low', '--low', metavar=-9999.0, type=float, default=-9999.0,
                         help='The lower (nm) for removing waters')
    parser.add_argument('-upp', '--upp', metavar=9999.0, type=float, default=9999.0,
                         help='The upper (nm) for removing waters')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit('Error! Missing input parameters')
    return parser.parse_args()

if __name__ == '__main__':
    option = parser_opt()
    ret = PymixMem(fp=option.file, fgro=option.gro, dx=option.dx, frac=option.frac, nlayer=option.nlayer)
    if option.file is not None:
        ret.gen_mem()
    if option.gro is not None:
        ret.remove_water(lower=option.low, upper=option.upp)
