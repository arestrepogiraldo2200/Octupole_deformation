# Libraries
import numpy as np
import pandas as pd
import subprocess

# ===================================================================================================================
# Function to invoke Schur program and perform products (λ1,μ1)x(λ2,μ2) 
def SU3prod(l1,m1,l2,m2):
    
    # Formatting parameters to Schur program
    if (l1+m1 >= 10):
        arg1 = "!"+str(l1+m1)
    else:
        arg1 = str(l1+m1)
    
    if (l2+m2 >= 10):
        arg2 = "!"+str(l2+m2)
    else: 
        arg2 = str(l2+m2)
        
    prod = f"prod {{{arg1} {m1}}},{{{arg2} {m2}}}\n"
    
    # Input file to deliver to Schur
    with open('schurinput.txt','w') as f:
        f.write("rep\n")
        f.write("gr1 su3\n")
        f.write(prod)
        f.write("^C")
    f.close()

    # Run Schur
    STR = subprocess.check_output("schur < schurinput.txt", shell=True, text=True)
    os.system("rm schurinput.txt")
    
    # Extract the irreps
    index = 0
    for i in STR:
        if i == "     ^":
            break
        index += 1
        
    irreps = ""
    for i in STR[23:index]:
        irreps += i
        
    # Separate each irrep
    irrepsSep = irreps.split(" + ")
        
    irrepsArr = []
    for i in irrepsSep:
        
        openkey = i.rfind("{")
        closekey = i.rfind("}")
        
        # Extract multiplicites
        mult = 1
        if openkey != 0:
            mult = int(i[0:openkey])
        
        l_m = i[openkey+1:closekey].split(" ")
            
        # Extract (λ1,μ1) values
        if (len(l_m)==3):
            irrepsArr.append([mult, int(l_m[0])-int(l_m[1]), int(l_m[1])])
            
        if (len(l_m)==2):
            if (l_m[1]==""):
                irrepsArr.append([mult, int(l_m[0]), 0])
            else:
                irrepsArr.append([mult, int(l_m[0])-int(l_m[1]), int(l_m[1])])

        if (len(l_m)==1):
            if (len(l_m[0])==2):
                irrepsArr.append([mult, int(l_m[0][0])-int(l_m[0][1]), int(l_m[0][1])])
            if (len(l_m[0])==1):
                irrepsArr.append([mult, int(l_m[0][0]), 0])
                
    values_df = pd.DataFrame(np.array(irrepsArr), columns = ["mult","lam","mu"])

    return values_df
    
# ===================================================================================================================    
# Function to invoke su3lib and compute Wigner coefficients
def SU3WignerCoeff(l1,m1,L1,l2,m2,L2,l3,m3,L3):
    
    # Input file to deliver to su3lib
    with open('su3wignercoeffs.txt','w') as f:
        f.write(f"{l1} {m1}\n")
        f.write(f"{L1}\n")
        f.write(f"{l2} {m2}\n")
        f.write(f"{L2}\n")        
        f.write(f"{l3} {m3}\n")
        f.write(f"{L3}\n")
    f.close()
    
    # Run SU3_SO3_WignerCoeffs
    
    
    lines = subprocess.check_output("/home/alejandro/Downloads/Thesis_Programs/su3libCorrection/su3lib/tools/SU3_SO3_WignerCoeffs < su3wignercoeffs.txt", shell=True, text=True)
    os.system("rm su3wignercoeffs.txt")

    # Values and labels stored
    fl = list(filter(lambda x: x != '', lines[8].split(" ")))
    flv = [fl[2], fl[3], fl[6], fl[7], fl[10:]]
    values = [[i for row in flv for i in row]]
    values = values + [list(filter(lambda x: x != '', l.split(" "))) for l in lines[9:]]
    
    # Column names and multiplicities
    cols = ["k1", "L1", "k2", "L2", "k3", "L3"]
    rhos = []
    for i in range(1,len(values[0])-5):
        cols.append("rho{0}".format(i))
        rhos.append("rho{0}".format(i))

    # Values stored
    values_df = pd.DataFrame(np.array(values), columns = cols)
    values_df[["k1", "L1", "k2", "L2", "k3", "L3"]] = values_df[["k1", "L1", "k2", "L2", "k3", "L3"]].astype(int)
    values_df[rhos] = values_df[rhos].astype(float)

    return values_df
        
# ===================================================================================================================    
# Function to invoke su3lib and compute 9-(λ,μ)
def SU39lammu(l1,m1,l2,m2,l12,m12,
              l3,m3,l4,m4,l34,m34,
              l13,m13,l24,m24,l,m):
    
    # Input file to deliver to su3lib
    with open('su39lammu.txt','w') as f:
        f.write(f"{l1} {m1}\n")
        f.write(f"{l2} {m2}\n")
        f.write(f"{l12} {m12}\n")
        f.write(f"{l3} {m3}\n")
        f.write(f"{l4} {m4}\n")
        f.write(f"{l34} {m34}\n")
        f.write(f"{l13} {m13}\n")
        f.write(f"{l24} {m24}\n")
        f.write(f"{l} {m}\n")
    f.close()
    
    # Run SU3_U9Coeffs
    lines = subprocess.check_output("/home/alejandro/Downloads/Thesis_Programs/su3libCorrection/su3lib/tools/SU3_U9Coeffs < su39lammu.txt", shell=True, text=True)
    os.system("rm su39lammu.txt")
        
    # Column names of the multiplicities
    cols = list(filter(lambda x: x != '', lines[20].split(" ")))
    cols.append("value")
    cols[5] = cols[5][:5]
    
    # Multiplicities values stored
    values = [list(filter(lambda x: x != '', l.split(" "))) for l in lines[21:]]
    values_df = pd.DataFrame(np.array(values), columns = cols)
    values_df["value"] = values_df["value"].astype(float)
    values_df[["rho1324","rho24", "rho13", "rho1234", "rho34", "rho12"]] = values_df[["rho1324","rho24", "rho13", "rho1234", "rho34", "rho12"]].astype(int)
    
    return values_df
    
# ===================================================================================================================        
# Function to invoke abljint and obtain the SU(3) irreps contained in U(N) 
def UNtoU3(eta, UNirrep):
    
    # Input file to deliver to su3lib
    with open('UNtoU3.txt','w') as f:
        f.write("1\n")
        f.write(f"{eta}\n")
        
        irr = ""
        for i in range(0, (eta+1)*(eta+2)//2):
            try:
                irr += str(UNirrep[i]) + " "
            except:
                irr += "0 "    
        
        f.write(f"{irr}\n")
        f.write("^C")
    f.close()

    # Run UNtoU3 algorithm
    lines = subprocess.check_output("/home/alejandro/Downloads/Thesis_Programs/Fortran_H.O/abljint < UNtoU3.txt", shell=True, text=True)
    os.system("rm UNtoU3.txt")
    
    # Extract the lines
    beginreps = lines.index("0            ( LM, MU)       NUM      DIM       C2       C3")
    endreps = lines[1:].index("1 DO YOU WANT DETAILED OUTPUT? ... YES (0) OR NO (1)")

    # Values and labels stored
    vv = [list(filter(lambda x: x != '', l.split(" "))) for l in lines[beginreps+1:endreps+1]]
    values = []
    for i in vv:
        values.append([int(i[1].rstrip(i[1][-1])), int(i[2].rstrip(i[2][-1])), int(i[3]), int(i[4]), int(i[5]), int(i[6])])

    # Values stored in dataframe
    values_df = pd.DataFrame(np.array(values), columns = ["lam","mu","num", "dim", "C2", "C3"])

    return values_df    
