#Converts a yaml file into FDS SPEC and REAC inputs.
#yaml file must have thermodynamic data inlcuding NASA7 or NASA9 polynomial and diameter and well-depth.
#run from command line as:
#python cantera2fds.py yamlfile.yaml background_species> fdsfile.fds
# where yamlfile.yaml is the name of the yaml file to process and background spcecies is the species name to set as BACKGROUND

import cantera as ct
import numpy as np
import sys
import math

gas = ct.Solution(sys.argv[1])
bg = sys.argv[2]
k_b = 1.380649E-23
t_r = 298.15
r0 = 8.314472

n_species = len(gas.species())

for i in range(n_species):
   if (bg==gas.species(i).name):
      bg_index = i
      break

for i2 in range(n_species):
   if (i2==0):
      i=bg_index
   elif (i2<=bg_index):
      i=i2-1
   elif (i2>bg_index):
      i=i2
   name = gas.species(i).name
   element_list = list(gas.species(i).composition.keys())
   atoms_list = list(gas.species(i).composition.values())
   sigma = gas.species(i).transport.diameter*1E10
   epsok = gas.species(i).transport.well_depth / k_b
   n_elem = len(atoms_list)
   formula = ''
   for j in range(n_elem):
      formula += element_list[j]+str(atoms_list[j])

   gas_list = list(gas.species(i).thermo.input_data.values())
   gasinit=gas.species(i).name+':1'
   gas.TPX=293.15,101325,gasinit

   Pr = gas.cp_mass * gas.viscosity / gas.thermal_conductivity
   poly = gas_list[0]
   temp_bands= gas_list[1]
   outstr="&SPEC ID='"+name+"',"
   if (i2==0):
      outstr+="BACKGROUND=T,"
   print(outstr)
   print("      PR_GAS=",int(Pr*1000)/1000,",")
   outstr="      FORMULA='"+formula+"',"
   print(outstr)
   print("      SIGMALJ=",int(sigma*1000)/1000,",")
   print("      EPSILONKLJ=",epsok,",")
   outstr="      POLYNOMIAL='"+poly+"',"
   print(outstr)
   band_l = len(temp_bands)
   outstr = ''
   for j in range(band_l):
      outstr+=str(temp_bands[j])+','
   print("      POLYNOMIAL_TEMP=",outstr)
   poly_len=9
   if (poly=='NASA7'):
      poly_len=7
   for j in range(band_l-1):
     if (j==0 and temp_bands[j]>t_r):
        t_r2 = t_r
        t_r = temp_bands[j]
        if (poly=='NASA7'):
           c_p = gas_list[2][j][0]+gas_list[2][j][1]*t_r+\
           gas_list[2][j][2]*t_r**2+gas_list[2][j][3]*t_r**3+gas_list[2][j][4]*t_r**4
           h_f = gas_list[2][j][0]*t_r+0.5*gas_list[2][j][1]*t_r**2+\
           1./3.*gas_list[2][j][2]*t_r**3+0.25*gas_list[2][j][3]*t_r**4+0.2*gas_list[2][j][4]*t_r**5+gas_list[2][j][5]
           h_f = h_f - c_p * (t_r - t_r2)
           h_f = int(h_f * r0 * 100)/100000.
           t_r = t_r2
        else:
           c_p = gas_list[2][j][0]/t_r**2+gas_list[2][j][1]/t_r+gas_list[2][j][2]+gas_list[2][j][3]*t_r+\
           gas_list[2][j][4]*t_r**2+gas_list[2][j][5]*t_r**3+gas_list[2][j][6]*t_r**4
           h_f = -gas_list[2][j][0]/t_r+gas_list[2][j][1]*math.log(t_r)+gas_list[2][j][2]*t_r+0.5*gas_list[2][j][3]*t_r**2+\
           1./3.*gas_list[2][j][4]*t_r**3+0.25*gas_list[2][j][5]*t_r**4+0.2*gas_list[2][j][6]*t_r**5+gas_list[2][j][7]
           h_f = h_f - c_p * (t_r - t_r2)
           h_f = int(h_f * r0 * 100)/100000.
           t_r = t_r2
     if (temp_bands[j]<t_r):
        if (poly=='NASA7'):
           h_f = gas_list[2][j][0]*t_r+0.5*gas_list[2][j][1]*t_r**2+\
           1./3.*gas_list[2][j][2]*t_r**3+0.25*gas_list[2][j][3]*t_r**4+0.2*gas_list[2][j][4]*t_r**5+gas_list[2][j][5]
           h_f = int(h_f * r0 * 100)/100000.
        else:
           h_f = -gas_list[2][j][0]/t_r+gas_list[2][j][1]*math.log(t_r)+gas_list[2][j][2]*t_r+0.5*gas_list[2][j][3]*t_r**2+\
           1./3.*gas_list[2][j][4]*t_r**3.+0.25*gas_list[2][j][5]*t_r**4+0.2*gas_list[2][j][6]*t_r**5+gas_list[2][j][7]
           h_f = int(h_f * r0 * 100)/100000.
     outstr = ''
     for k in range(poly_len):
        outstr+=str(gas_list[2][j][k])+','
     namestr='      POLYNOMIAL_COEFF(1:'+str(k+1)+','+str(j+1)+")="
     print(namestr,outstr)
   print("      ENTHALPY_OF_FORMATION=",h_f,"/")

numreac = len(gas.reactions())

A=[]
Ea=[]
b=[]
rlist=[]
plist=[]
three=[]
efflist=[]
explist=[]
for i in range(numreac):
	rlist.append(list(gas.reaction(i).reactants.items()))
	explist.append(0)
	for j in range(len(list(gas.reaction(i).reactants.items()))):
		explist[i]=explist[i]+list(gas.reaction(i).reactants.items())[j][1]
	plist.append(list(gas.reaction(i).products.items()))
	try:
		rate=gas.reaction(i).rate.input_data['low-P-rate-constant']
		A.append(rate['A'])
		Ea.append(rate['Ea'])
		b.append(rate['b'])

	except (KeyError):
		rate=gas.reaction(i).rate.input_data['rate-constant']
		A.append(rate['A'])
		Ea.append(rate['Ea'])
		b.append(rate['b'])

	if gas.reaction(i).reaction_type=='three-body' or gas.reaction(i).reaction_type=='three-body-Arrhenius':
		three.append(True)
		efflist.append(list(gas.reaction(i).efficiencies.items()))
	else:
		three.append(False)
		efflist.append([])

for i in range(len(rlist)):
	print(f"&REAC ID='R{i+1}',")
	if (gas.reaction(i).reversible):
		print("     REVERSE=T,")
	print("     RADIATIVE_FRACTION=0,")
	if (three[i]):
		print("     THIRD_BODY=T,")
		print("     A=",int(A[i]*10000*1000**(explist[i]))/1E4,",")
		if (len(efflist[i])>0):
			effs=str(np.array(efflist[i])[:,0])
			effn=str(np.array(efflist[i])[:,1])
			effs=effs.replace('[','').replace(']','').replace("' '","','").replace("\n",",")
			effn=effn.replace('[','').replace(']','').replace(" ",",").replace("'","").replace("\n",",")
			print("     THIRD_EFF_ID=",effs,",")
			print("     THIRD_EFF=",effn,",")
	else:
		print("     A=",int(A[i]*10000*1000**(explist[i]-1))/1E4,",")
	print("     E=",int(Ea[i]*1E6)/1E9,",")
	if (b[i]!=0): print("     N_T=",b[i],",")
	rlist2=str(np.array(rlist[i])[:,0])
	plist2=str(np.array(plist[i])[:,0])
	rlist2=rlist2.replace('[','').replace(']','').replace("' '","','")
	plist2=plist2.replace('[','').replace(']','').replace("' '","','")
	print("     SPEC_ID_NU=",rlist2,",",plist2,",")
	rlist2=str(np.array(rlist[i])[:,1])
	plist2=str(np.array(plist[i])[:,1])
	rlist2=rlist2.replace('[','-').replace(']','').replace(" ",",-").replace("'","")
	plist2=plist2.replace('[','').replace(']','').replace(" ",",").replace("'","")
	print("     NU=",rlist2,",",plist2,"/")



#gas.species(1).name
#gas.species(1).composition
#{'element':val}
#list(gas.species().composition)
#len elements
#gas.speces().thermo.coeffs
#gas.species().thermo.min_temp  max_temp
#gas.species().transport.diameter * 1E10
#gas.species().transport.well_depth / k_b


#x=list(gas.species().thermo.input_data)
#x[0]=poly
#x[1]=temp range
#x[2][0:N] = poly data for range

