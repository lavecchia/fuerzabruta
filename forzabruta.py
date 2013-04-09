#!/usr/bin/env python
# forzabruta racional v. 1.0.0
# algoritmo para generar conformaciones a partir de la variacion de parametros. El directorio debe contener:
# - archivo de geometria (pdb) Ejemplo. estructura.pdb
# - archivo de parametros (.in). Ejemplo. estructura.in.
#
# En este ultimo se listaran los diedros a modificar D atom1 atom2 atom3 atom4. Si se desea que un diedro varie entre
# valores en particular hagalo luego de esta forma: D atom1 atom2 atom3 atom4 = val1 val2 val3 val4 ...

import string
import os
import datetime
import sys
from pymol import cmd
import pymol
from optparse import OptionParser
import subprocess


# extraido de http://doeidoei.wordpress.com/2009/02/11/pymol-api-simple-example/
import __main__ 
__main__.pymol_argv = [ 'pymol', '-qcr'] # Quiet and no GUI
pymol.finish_launching()

#CONSTANTES
ANGLESTEP = 20
MINMAXNUMBER = 5 #numero maximo de minimos por angulo diedro
FORCEFIELD = "pm6 1scf" #campo de fuerza utilizado en MOPAC
MOPACPATH = "/opt/mopac/MOPAC2012.exe" #ruta del MOPAC


#VARIABLES
atomIDs = [] #almacena los Id de los atomos a modificar los parametros
paramod = [] #almacen parametros modificables 
matrizval = [] #almacena por cada columna los valores respectivos para cada archivo de salida



#FUNCIONES

#funcion que devuelve perfil de energia
def calc_energyprofile(atom1, atom2, atom3, atom4, structurefilename, forcefield, anglestep, mopacpath):
    energyprofilelist=[]
    pymol.cmd.load(os.path.abspath(os.curdir)+ "//" + structurefilename)
    for value in range (-180,180,anglestep):
        cmd.set_dihedral("ID " + str(atom1),"ID " + str(atom2),"ID " + str(atom3),"ID " + str(atom4), value)
        tmpstructurefilename = structurefilename.replace(".pdb", "_tmp.pdb")
        moptmpstructurefilename = structurefilename.replace(".pdb", "_tmp.mop")
        outtmpstructurefilename = structurefilename.replace(".pdb", "_tmp.out")
        cmd.save(tmpstructurefilename)

        #devnull = open('/dev/null', 'w') #para ocultar stdout
        command = "babel -ipdb " + tmpstructurefilename + " -omop " + moptmpstructurefilename + " -xk '" + forcefield + "' 2> /dev/null; " + mopacpath + " " + moptmpstructurefilename
        #process = subprocess.Popen(command , shell=True, bufsize=2048, stdout=devnull, stderr=devnull)
        os.system(command)
        #process.wait()
        #devnull.close()
        outtmpstructurefile = open(outtmpstructurefilename, "r")
        for line in outtmpstructurefile:
            if 'FINAL HEAT' in line:
                energyprofilelist.append([value,float(line.split()[5])])
        outtmpstructurefile.close()
    return energyprofilelist

#calcula pendiente en cada punto tomando el cociente incremental entre i-1 y i+1
def calc_slope(energyprofilelist):
    slopelist=[]
    for i in range(0,len(energyprofilelist)):
        if i == 0: #primer punto
            slope=(energyprofilelist[1][1]-energyprofilelist[-1][1])/(energyprofilelist[1][0]-(energyprofilelist[-1][0]-360))
        elif i == len(energyprofilelist) - 1: #ultimo punto
            slope=(energyprofilelist[1][1]-energyprofilelist[i-1][1])/(energyprofilelist[1][0]+360-energyprofilelist[i-1][0])       
        else:
            slope=(energyprofilelist[i+1][1]-energyprofilelist[i-1][1])/(energyprofilelist[i+1][0]-energyprofilelist[i-1][0])
        slopelist.append([energyprofilelist[i][0],energyprofilelist[i][1],slope])
    return slopelist
    
#funcion que devuelve en base al perfil de energia los minimos
# slopelist lista [elemento1:diedro, elemento2:energia, elemento3:pendiente]
def find_min(slopelist):
    minlist = []
    for i in range(0,len(slopelist)):
        if (i == 0) and (slopelist[-1][2]<0) and (slopelist[0][2]>0):
            #estima minimo y asigna el valor de energia
            minvar = (slopelist[-1][0]-360)*abs(slopelist[0][2])/(abs(slopelist[0][2])+abs(slopelist[-1][2]))+(slopelist[0][0])*abs(slopelist[-1][2])/(abs(slopelist[-1][2])+abs(slopelist[0][2]))
            if minvar < -180:
                minvar = minvar + 360
            minlist.append([minvar, min(slopelist[0][1],slopelist[-1][1])])
        elif (slopelist[i-1][2]<0) and (slopelist[i][2]>0):
            #estima minimo y asigna el valor de energia
            minvar = slopelist[i][0]*abs(slopelist[i-1][2])/(abs(slopelist[i-1][2])+abs(slopelist[i][2]))+slopelist[i-1][0]*abs(slopelist[i][2])/(abs(slopelist[i][2])+abs(slopelist[i-1][2]))
            minlist.append([minvar, min(slopelist[i-1][1],slopelist[i][1])])
    return minlist
    
def set_min(atom1, atom2, atom3, atom4, structurefilename, forcefield, anglestep, minmaxnumber, mopacpath):
    minselect = []
    minlist = find_min(calc_slope(calc_energyprofile(atom1, atom2, atom3, atom4, structurefilename, forcefield, anglestep, mopacpath)))
    minlist = sorted(minlist,key=lambda x: x[1])
    #si el numero de minimos supera el numero maximo determinado en minmaxnumber
    
    for i in range(0,len(minlist)):
        if i>minmaxnumber:
            print "Warning: el numero de minimos hallados para " + str(atom1) + "-" + str(atom2) + "-" + str(atom3) + "-" + str(atom4) + " supera el numero de minimos establecidos como maximos. Acotando..."
            break
        minselect.append(minlist[i][0])
    return minselect
    

#PARSEO
usage = "%prog [opciones]"
parser = OptionParser(usage=usage, version="%prog 1.0")
parser.add_option('-i','--infile', action='store', type ='string', dest='infile', metavar='nombre', help='especifica el nombre de la entrada .pdb\n')
#parser.add_option('-p','--param', action='store', type ='string', dest='parametros', metavar='parametros', help='especifica archivo con parametros a modificar. por defecto nombreentrada.in" \n')
(options, args) = parser.parse_args()

# Especifica nombre de entrada    
if options.infile:
    structurefile = options.infile
    structurefilename = structurefile
else:
    print "Debe especificar un archivo de entrada. Ej. forzabruta -i ejemplo.pdb"
    sys.exit()


infilenamediv = structurefile.split(".")
infilename = infilenamediv[0]+ ".in"

    
#CODIGO
#1. Lectura de archivo con parametros a modificar. Genera Matriz con valores de diedros variables

pymol.cmd.load(os.path.abspath(os.curdir)+ "//" + structurefile)

infile = open(infilename,"r") #abre archivo con parametros a modificar

nroparam = 0
for line in infile:
    paramod.append([])
    linediv = line.split("=") #corta la linea
    line1div = linediv[0].split()
    if line1div[0] == "D":
        atom1 = line1div[1]
        atom2 = line1div[2]
        atom3 = line1div[3]
        atom4 = line1div[4]
    if len(linediv)>1:
        line2div = linediv[1].split()
        nrovalores=len(line2div) #calcula cuantos valores va a tomar este parametro
        if nrovalores > 0: #si existe algun valor...
            for valor in range(0,nrovalores):
                paramod[nroparam].append(line2div[valor])                
    else: #sino calcula mediante el perfil energetico que valores tomar
        paramod[nroparam] = set_min(atom1, atom2, atom3, atom4, structurefilename, FORCEFIELD, ANGLESTEP, MINMAXNUMBER, MOPACPATH)
        #paramod[nroparam]=variacion #asigna variacion por defecto
    atomIDs.append("ID " + atom1 + ", ID " + atom2 + ", ID " + atom3 + ", ID " + atom4) #guarda en atomIDs los atomos a modificar en cada linea
    nroparam = nroparam + 1
infile.close()
pymol.cmd.load(os.path.abspath(os.curdir)+ "//" + structurefile)

#2. Genera matriz matrizval nxm donde n es el numero de archivos a generar y m los parametros a modificar
nrosalidas=1
for n in range(0,len(paramod)):
    nrosalidas = nrosalidas*len(paramod[n])

if len(paramod)>0:
    print "Parametros a modificar: " + str(nroparam)
    print "Salidas a generar: " + str(nrosalidas)
    for n in range(0,nrosalidas):
        matrizval.append([])
    for npam in range(0,len(paramod)):
        nsal = 0
        for p in range(0,nrosalidas/len(paramod[npam])):
            for valor in paramod[npam]:
                matrizval[nsal].append(valor)
                nsal = nsal + 1
        matrizval.sort() #importante para el programa


    #3. Genera Archivos de Salidas
    infilenamediv=infilename.split(".")
    for nsal in range(0,nrosalidas):
        salida = [] #lista que contiene las lineas para un archivo de salida especifico
        npam = 0
        extname = infilenamediv[0]
        for line in atomIDs: 
            linediv=line.split(",")
            val = float(matrizval[int(nsal)][int(npam)])
            cmd.set_dihedral(linediv[0],linediv[1],linediv[2],linediv[3], val)
            extname = extname + "_D" + str(npam) + "_" + str(int(val)) #nombre de salida con valores de parametros            
            npam = npam + 1
        cmd.set('pdb_conect_all')
        cmd.save(extname+".pdb")
        print "Se creo el archivo " + extname+".pdb"
else:
    print "No se encuentran parametros para modificar"


# Get out!
pymol.cmd.quit()





    
