#!/bin/python

import sys
import glob 
sys.path.append('/home/prlarsen/usr/genericio/python')
import genericio as gio

num1 = int(sys.argv[1])
num2 = int(sys.argv[2])

file_list = glob.glob('/eagle/LastJourney/heitmann/OuterRim/M000/L4225/HACC000/analysis/Particles/STEP247/m000.mpicosmo.247#*')
ids = [152616279983,39674159997]

for i in range(num1,num2):
    print(i)
    file1 = file_list[i]
    id_list = gio.gio_read(file1,'id')
    if ids[1] in id_list:
        print(file1)
    if ids[0] in id_list:
        print(file1)

halo_list = np.array([152616279983,  39674159997,  37787265719, 123808104886,
       117622284790, 265355933680, 241767785788, 192580203025,
        92348040658, 257603310439, 237472094696,  42526244172,
        33633148914, 414730779187, 326109211825, 499440052482,
       492653215375, 449172439588, 366970255881, 349133183051,
       292503301288, 612707547911, 786957140741, 598859535685,
       740122031447, 136538713701, 265110944511, 236063269900,
        84802265492, 135957498852])

file_list = glob.glob('/eagle/LastJourney/heitmann/OuterRim/M000/L4225/HACC000/analysis/Halos/HaloCatalog/03_31_2018.OR.253.fofproperties#*')
for i in range(len(file_list)):
    id_list = gio.gio_read(file_list[i],'fof_halo_tag')
    mask = np.isin(halo_list,id_list)
    
