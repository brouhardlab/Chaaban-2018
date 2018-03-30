
# coding: utf-8

# In[13]:

##########################################################################################
###############   ALL COMPARISON FILENAMES ARE IN SUBSEQUENT CELLS #######################
##########################################################################################

###############
# PDB FILE
###############
#ORDER OF CHAINS HAS TO BE BOTTOM MONOMER TO TOP MONOMER, AND IDENTICAL FOR BOTH PDB FILES

###############
# ALIGNMENT FILE
###############
#USE ALIGN CHAIN SEQUENCES IN CHIMERA TO GENERATE .ASC FILE
#EDIT SO THAT IT HAS ONLY ONE ALIGNMENT AND NO HEADERS

###############
# OUTPUT FILE
###############
#.BILD FILE WITH VECTOR INFORMATION FOR CHIMERA

###################
######TEST FILES
#pdb1 = 'PF-3JAT.pdb'
#pdb2 = 'PF-3JAS-ali-to-3JAT-bottombeta.pdb'
#aligna = '3JAT-to-3JAS-alignment-alpha.asc'
#alignb = '3JAT-to-3JAS-alignment-beta.asc'
#RUN:
#python Tubulin_Vectors.py --pdb1 PF-3JAT.pdb --pdb2 PF-3JAS-ali-to-3JAT-bottombeta.pdb --a 3JAT-to-3JAS-alignment-alpha.asc --b 3JAT-to-3JAS-alignment-beta.asc
##################

##########################################################################################
##########################################################################################

import pandas as pd
import numpy as np
import optparse
import sys

#Arrow size limits
max_cutoff = 6 #Don't show vectors longer than 6A
min_cutoff = 1 #Don't show vectors shorter than 1A

#Arrow radius
r1 = 0.15
r2 = 0.35
rho = 0.5

#Color arrows by domain
colordomain_flag = 0
ignorecterm_flag = 0 #Don't show

#If colordomain_flag = 1 or ignorecterm_flag = 1,
#at what residue does the intermediate domain and the c-terminal domain start?
intermstart = 206
ctermstart = 380

def setupParserOptions():
    parser = optparse.OptionParser()
    parser.add_option("--pdb1",dest="pdb1",type="string",metavar="FILE",
                help="PF-3JAT.pdb")
    parser.add_option("--pdb2",dest="pdb2",type="string",metavar="FILE",
        help="PF-3JAS-ali-to-3JAT-bottombeta.pdb")
    parser.add_option("--a",dest="alpha",type="string",metavar="FILE",
                help="3JAT-to-3JAS-alignment-alpha.asc")
    parser.add_option("--b",dest="beta",type="string",metavar="FILE",
        help="3JAT-to-3JAS-alignment-beta.asc")
    parser.add_option("--reverse", action="store_true",dest="reverse",default=False,
        help="Reverse direction of arrows if argument is passed")

    options,args = parser.parse_args()

    if len(args) > 1:
            parser.error("Unknown commandline options: " +str(args))

    if len(sys.argv) < 2:
            parser.print_help()
            sys.exit()

    params={}

    for i in parser.option_list:
            if isinstance(i.dest,str):
                    params[i.dest] = getattr(options,i.dest)
    return params

def draw_arrows(pdb1, pdb2, aligna, alignb, reverseflag):


    df1 = pd.read_csv(pdb1, sep='\t', header=None, skipfooter = 1)
    df2 = pd.read_csv(pdb2, sep='\t', header=None, skipfooter = 1)
    dfa = pd.read_csv(aligna, sep='\t', header=None, skipfooter = 1)
    dfb = pd.read_csv(alignb, sep='\t', header=None, skipfooter = 1)

    alignment_a1 = []
    alignment_a2 = []
    alignment_b1 = []
    alignment_b2 = []

    for index, row in dfa.iterrows():
        if index < 10:
            alignment_a1.append(int(row[0][1:3]))
            alignment_a2.append(int(row[0][-4:-2]))
        else:
            alignment_a1.append(int(row[0][1:5]))
            alignment_a2.append(int(row[0][-5:-2]))

    for index, row in dfb.iterrows():
        if index < 10:
            alignment_b1.append(int(row[0][1:3]))
            alignment_b2.append(int(row[0][-4:-2]))
        else:
            alignment_b1.append(int(row[0][1:5]))
            alignment_b2.append(int(row[0][-5:-2]))
    
    ####################################################

    coords1_1 = []
    coords1_2 = []
    coords1_3 = []
    coords1_4 = []
    coords2_1 = []
    coords2_2 = []
    coords2_3 = []
    coords2_4 = []

    chain = 0
    previous = 0

    for index, row in df1.iterrows():
        if (row[0][0:5] == "SHEET" or row[0][0:5] == "HELIX"):
            continue
        if int(row[0][23:26]) == 1 and previous != 1:
            chain += 1
        if row[0][13:15] == "CA":
            n = int(row[0][23:26])
            x = float(row[0][31:38])
            y = float(row[0][39:46])
            z = float(row[0][47:54])
            if chain == 1:
                coords1_1.append([n, x, y, z])
            if chain == 2:
                coords1_2.append([n, x, y, z])
            if chain == 3:
                coords1_3.append([n, x, y, z])
            if chain == 4:
                coords1_4.append([n, x, y, z])
        previous = int(row[0][23:26])

    chain = 0
    previous = 0

    for index, row in df2.iterrows():
        if (row[0][0:5] == "SHEET" or row[0][0:5] == "HELIX"):
            continue
        if int(row[0][23:26]) == 1 and previous != 1:
            chain += 1
        if row[0][13:15] == "CA":
            n = int(row[0][23:26])
            x = float(row[0][31:38])
            y = float(row[0][39:46])
            z = float(row[0][47:54])
            if chain == 1:
                coords2_1.append([n, x, y, z])
            if chain == 2:
                coords2_2.append([n, x, y, z])
            if chain == 3:
                coords2_3.append([n, x, y, z])
            if chain == 4:
                coords2_4.append([n, x, y, z])
        previous = int(row[0][23:26])

    ########################

    arrows = []

    for i in np.arange(1, len(alignment_a1)):

        for loc1,item in enumerate(coords1_1):
            if item[0] == alignment_a1[i]:
                break

        for loc2,item in enumerate(coords2_1):
            if item[0] == alignment_a2[i]:
                break

        coord1 = coords1_1[loc1]
        coord2 = coords2_1[loc2]

        n1, x1, y1, z1 = coord1[0], coord1[1], coord1[2], coord1[3]
        n2, x2, y2, z2 = coord2[0], coord2[1], coord2[2], coord2[3]

        dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        if dist > max_cutoff or dist < min_cutoff:
            continue
         
        if colordomain_flag == 1:
            if n1 < intermstart:
                arrows.append(".color red")
            elif n1 < ctermstart:
                arrows.append(".color blue")
            else:
                if ignorecterm_flag == 1:
                    break
                arrows.append(".color green")

        if reverseflag == 0:
            arrows.append(".arrow " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(r1) + " " + str(r2) + " " + str(rho))
        else:
            arrows.append(".arrow " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(r1) + " " + str(r2) + " " + str(rho))

    for i in np.arange(1, len(alignment_b1)):

        for loc1,item in enumerate(coords1_2):
            if item[0] == alignment_b1[i]:
                break

        for loc2,item in enumerate(coords2_2):
            if item[0] == alignment_b2[i]:
                break

        coord1 = coords1_2[loc1]
        coord2 = coords2_2[loc2]

        n1, x1, y1, z1 = coord1[0], coord1[1], coord1[2], coord1[3]
        n2, x2, y2, z2 = coord2[0], coord2[1], coord2[2], coord2[3]

        dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        if dist > max_cutoff or dist < min_cutoff:
            continue
            
        if colordomain_flag == 1:
            if n1 < intermstart:
                arrows.append(".color red")
            elif n1 < ctermstart:
                arrows.append(".color blue")
            else:
                if ignorecterm_flag == 1:
                    break
                arrows.append(".color green")

        if reverseflag == 0:
            arrows.append(".arrow " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(r1) + " " + str(r2) + " " + str(rho))
        else:
            arrows.append(".arrow " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(r1) + " " + str(r2) + " " + str(rho))

    for i in np.arange(1, len(alignment_a1)):

        for loc1,item in enumerate(coords1_3):
            if item[0] == alignment_a1[i]:
                break

        for loc2,item in enumerate(coords2_3):
            if item[0] == alignment_a2[i]:
                break

        coord1 = coords1_3[loc1]
        coord2 = coords2_3[loc2]

        n1, x1, y1, z1 = coord1[0], coord1[1], coord1[2], coord1[3]
        n2, x2, y2, z2 = coord2[0], coord2[1], coord2[2], coord2[3]

        dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        
        if dist > max_cutoff or dist < min_cutoff:
            continue
            
        if colordomain_flag == 1:
            if n1 < intermstart:
                arrows.append(".color red")
            elif n1 < ctermstart:
                arrows.append(".color blue")
            else:
                if ignorecterm_flag == 1:
                    break
                arrows.append(".color green")

        if reverseflag == 0:
            arrows.append(".arrow " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(r1) + " " + str(r2) + " " + str(rho))
        else:
            arrows.append(".arrow " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(r1) + " " + str(r2) + " " + str(rho))

    for i in np.arange(1, len(alignment_b1)):

        for loc1,item in enumerate(coords1_4):
            if item[0] == alignment_b1[i]:
                break

        for loc2,item in enumerate(coords2_4):
            if item[0] == alignment_b2[i]:
                break

        coord1 = coords1_4[loc1]
        coord2 = coords2_4[loc2]

        n1, x1, y1, z1 = coord1[0], coord1[1], coord1[2], coord1[3]
        n2, x2, y2, z2 = coord2[0], coord2[1], coord2[2], coord2[3]

        dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        if dist > max_cutoff or dist < min_cutoff:
            continue
            
        if colordomain_flag == 1:
            if n1 < intermstart:
                arrows.append(".color red")
            elif n1 < ctermstart:
                arrows.append(".color blue")
            else:
                if ignorecterm_flag == 1:
                    break
                arrows.append(".color green")

        if reverseflag == 0:
            arrows.append(".arrow " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(r1) + " " + str(r2) + " " + str(rho))
        else:
            arrows.append(".arrow " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(r1) + " " + str(r2) + " " + str(rho))

    if reverseflag == 0:
        outname = pdb1[0:-4] + '_to_' + pdb2[0:-4] + '.bild'
    else:
        outname = pdb2[0:-4] + '_to_' + pdb1[0:-4] + '.bild'
    
    exportfile = open(outname, 'w')
    for item in arrows:
        exportfile.write("%s\n" % item)
    exportfile.close()

    print('\n\nSuccessfully written:\n' + outname + '\n')
    
def mainloop(params):

    pdb1 = params["pdb1"]
    pdb2 = params["pdb2"]
    aligna = params["alpha"]
    alignb = params["beta"]
    if params["reverse"] == True:
        reverseflag = 1
    else:
        reverseflag = 0
        
    draw_arrows(pdb1, pdb2, aligna, alignb, reverseflag)
    
if __name__ == "__main__":
    params = setupParserOptions()
    mainloop(params)


# In[ ]:



