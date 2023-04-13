from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from multiprocessing import Process
from MDAnalysis import *
import MDAnalysis.analysis.pca as pca
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms
import os
import matplotlib as mpl
from scipy.ndimage.filters import gaussian_filter1d

def Dist(a, b):
    r = np.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))
    return r


def Align(ref, pdb, traj, sel, file_new):
    print('Universe loading...', traj)
    u = Universe(pdb, traj)
    refe = Universe(ref)  # Universe(pdb)
    print('Alignment run ', traj)
    alignment = align.AlignTraj(u, refe, filename=file_new, select=sel)
    alignment.run()
    print('Alignment finish ', traj)


def RMSD_calc(pdb, traj, sel, out):
    print('Universe loading...', traj)
    u = Universe(pdb, traj)
    print('Rmsd calculating...', traj)
    R = rms.RMSD(u, select=sel)
    R.run()
    rmsd = R.rmsd
    np.savetxt(out, rmsd)
    print('RMSD calculation finished ', traj)


def RMSF_calc(pdb, traj, out):
    f = open(out, 'w')
    print('Universe loading...', traj)
    u = Universe(pdb, traj)
    print('Rmsf calculating...', traj)
    s = u.select_atoms('protein and name CA')
    R = rms.RMSF(s, verbose=True).run()
    N = s.resnums
    M = R.rmsf
    for n, m in zip(N, M):
        # print(n)
        f.write(str(n) + ' ' + str(m) + '\n')
    f.close()
    print('RMSF replica done!', traj)


def Gates(pdb, traj, gLoopGlyID, leuID, pheID, stride, out):
    f = open(out, 'w')  # order: no Glycine_AE Leucine_AE Phenyloalanine_AE Glycine_CE ...
    u = Universe(pdb, traj)
    print(u.atoms[:5])
    i = 1
    print('Gates calculating...', traj)
    for t in u.trajectory[::stride]:
        line = str(i) + ' '
        i += 1
        for chain1, chain2 in zip(['PROA', 'PROB'], ['PROC', 'PROD']):  # NEED TO BE ABLE TO IDENTIFY CHAINS-prepare pdb files with named chains

            # Glycine in G loop
            G1 = u.select_atoms('protein and segid ' + chain1 + f' and resid {gLoopGlyID}')
            G2 = u.select_atoms('protein and segid ' + chain2 + f' and resid {gLoopGlyID}')
            ming = 100.
            for g1 in G1:
                for g2 in G2:
                    if Dist(g1.position, g2.position) < ming: ming = Dist(g1.position, g2.position)

            line += str(ming) + ' '

            # leucine
            G1 = u.select_atoms('protein and segid ' + chain1 + f' and resid {leuID}')
            G2 = u.select_atoms('protein and segid ' + chain2 + f' and resid {leuID}')
            ming = 100.
            for g1 in G1:
                for g2 in G2:
                    if Dist(g1.position, g2.position) < ming: ming = Dist(g1.position, g2.position)

            line += str(ming) + ' '

            # phenylalanine
            G1 = u.select_atoms('protein and segid ' + chain1 + f' and resid {pheID}')
            G2 = u.select_atoms('protein and segid ' + chain2 + f' and resid {pheID}')
            ming = 100.
            for g1 in G1:
                for g2 in G2:
                    if Dist(g1.position, g2.position) < ming: ming = Dist(g1.position, g2.position)

            line += str(ming) + ' '

        line += '\n'
        print(line)
        f.write(line)




def Correlation(pdb, traj, out, stride):
    chainCorr = []
    for chain in ["A", "B", "C", "D"]:
        sel = f"segid PRO{chain}"  # "name CA"# and (((segid A C E G) and resid 1:54 353:390) or ((segid B D F H) and resid 192:286 620:677 917:994) )"
        u = Universe(pdb, traj)
        end = u.trajectory.n_frames
        beg = round(end / 10)
        #stride = int((end - beg)/100.)
        stride = stride
        s = u.select_atoms(sel)
        # for i in s.residues:
        #	print((i.resid))
        residnum = s.n_residues
        R = np.zeros(shape=(residnum, residnum))
        number_of_frames = int(end) - int(beg)
        end_of_frames = int(end)
        beg_of_frames = int(beg)
        # print(u.select_atoms(sel).residues[0].atoms.center_of_mass())
        residuas1 = [[i.atoms.center_of_mass() for i in u.select_atoms(sel).residues] for frame in
                     u.trajectory[beg_of_frames:end_of_frames:stride]]
        # COMs = []
        # for i in u.select_atoms(sel).residues:
        #	COMs.append(i)
        # for frame in u.trajectory[beg_of_frames:end_of_frames:stride]:
        #	pass
        # ^ the average position of the CA (COM) of each residue for each frame

        # residuas2 = [ [i.atoms.center_of_mass() for i in u.select_atoms(sel).residues] for frame in u.trajectory[beg_of_frames:end_of_frames:stride]]

        print(np.shape(residuas1))
        nr_res = np.shape(residuas1)[1]

        print("Starts the main loop...")
        for i in range(nr_res):
            av_res1 = np.array([0.0, 0.0, 0.0])

            for frame in range(int(number_of_frames / stride)):
                res1 = residuas1[frame][i]
                av_res1 += res1
            av_res1 /= (number_of_frames / stride)
            for j in range(nr_res):
                av_res2 = np.array([0.0, 0.0, 0.0])
                for frame in range(int(number_of_frames / stride)):
                    res2 = residuas1[frame][j]
                    av_res2 += res2
                av_res2 /= (number_of_frames / stride)

                # implenetation of dynamic cross correlation formula
                DR1 = np.array([0.0, 0.0, 0.0])
                DR2 = np.array([0.0, 0.0, 0.0])
                L = 0.0  # numerator
                M1 = 0.0  # denominator 1
                M2 = 0.0  # debinubatir 1

                for frame in range(
                        int(number_of_frames / stride)):  # Tutaj licze cross correlacje zgodnie z formula z DCC
                    res1 = residuas1[frame][i]
                    res2 = residuas1[frame][j]
                    # print(res1,res2)
                    DR1 = res1 - av_res1
                    DR2 = res2 - av_res2
                    L += DR1[0] * DR2[0] + DR1[1] * DR2[1] + DR1[2] * DR2[2]
                    M1 += DR1[0] * DR1[0] + DR1[1] * DR1[1] + DR1[2] * DR1[2]
                    M2 += DR2[0] * DR2[0] + DR2[1] * DR2[1] + DR2[2] * DR2[2]
                L /= (number_of_frames / stride)
                M1 /= (number_of_frames / stride)
                M2 /= (number_of_frames / stride)
                M1 = 1.0 / np.sqrt(M1)
                M2 = 1.0 / np.sqrt(M2)
                R[i, j] = L * M1 * M2
            print("  Progress: %s" % ((i + 1) * 100.0 / nr_res))
        chainCorr.append(R)
    print(np.shape(chainCorr))
    # need to calculate the average correlation matrix here
    averageChainCorr = np.mean(chainCorr, axis=0)
    print(np.shape(averageChainCorr))
    np.savetxt(out, averageChainCorr)


def PlotRMSF(RMSFPath, numOfChains=4, figsize=(8, 4)):
    RMSF = np.loadtxt(RMSFPath)
    RMSF = RMSF.T  # for readability - I'm more of a row man than a column man, personally
    numOfChains = numOfChains  # can be changed, but we are currently dealing with a tetramer in the case of Kir
    residueNums = list(np.array_split(RMSF[0], numOfChains)[0])  # since this contains four identical lists, only one is needed
    chainsRMSF = np.array_split(RMSF[1], numOfChains)
    fig = plt.figure(figsize=figsize)

    for chainNum in range(numOfChains):
        plt.plot(residueNums, chainsRMSF[chainNum], label=f"Chain {chainNum + 1}")

    averageChainRMSF = list(
        np.zeros(len(residueNums)))  # getting the length of one chain and creating an array of zeros

    for chainNum in range(numOfChains):
        for residNum in range(len(residueNums)):
            averageChainRMSF[residNum] += chainsRMSF[chainNum][residNum]
    averageChainRMSF = [x / numOfChains for x in averageChainRMSF]

    plt.plot(residueNums, averageChainRMSF, "k-", label="Average")

    plt.legend(loc="best")
    plt.xlabel("Residue No.")
    plt.ylabel(r"RMSF")


def CorrelationPlot(corrPath):
    m = 'RdBu'
    # Zestaw 1

    corrData = np.loadtxt(corrPath)

    fig, ax = subplots(1, 1)
    img = ax.imshow(corrData, cmap=m, vmin=-1, vmax=1)

    fig.colorbar(img)



def RMSDPlot(rmsdPath, figsize=(4, 4)):
    rmsd = np.loadtxt(rmsdPath)
    rmsd = rmsd.T  # transpose makes it easier for plotting
    time = rmsd[1]
    fig = plt.figure(figsize=figsize)
    plt.plot(time, rmsd[2], 'k-', label="all")
    plt.legend(loc="best")
    plt.xlabel("time (ps)")
    plt.ylabel(r"RMSD ($\AA$)")


def phist(vector, range, color, bin):
    h = histogram(vector, bins=bin, range=range, density=True)
    X = h[1][:-1]
    Y = gaussian_filter1d(h[0], sigma=2)
    plot(X, Y, color=color, label=' ', linewidth=3)
    fill_between(X, 0, Y, facecolor=color, alpha=0.3)

