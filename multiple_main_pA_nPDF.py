#!/usr/bin/python3.8

import os
import sys
from time import sleep


options = ['SoftQCD:all=on ', 'SoftQCD:inelastic=on ', 'HardQCD:allwdi=on ', 'HardQCD:all=on ', #0-3
           'HardQCD:onlydi=on ', 'PromptPhoton:all=on ', 'PhaseSpace:pTHatMin=0. ', 'PDF:pSet=LHAPDF6:cteq6l1 ',#4-7
           'MultipartonInteractions:Kfactor=0.5 ']#8
names = ['soft', 'inelastic', 'allwdi', 'hard', 'onlydi', 'promt', 'pthat', 'nPDFSetA7', 'Kfac']
system = ['pp1', 'np3', 'pn3', 'nn3', 'pAl1', 'pAu1', 'dAu1', 'HeAu1']
nPSFs = ['nNNPDF10_nlo_as_0118_N1 ', 'nCTEQ15_3_2 ', 'nCTEQ15_27_13 ', 'nCTEQ15_197_79 ', #0-3
         'nNNPDF10_nlo_as_0118_D2 ', 'nNNPDF10_nlo_as_0118_Al27 ', 'nNNPDF10_nlo_as_0118_Au197 ', #4-6
         'EPPS16nlo_CT14nlo_Al27 ', 'EPPS16nlo_CT14nlo_Au197 ', '0 '] #7-9
pdf_name = ['nCTEQ15_1_1', 'nCTEQ15_3_2', 'nCTEQ15_27_13', 'nCTEQ15_197_79',
            'nNNPDF10_nlo_as_0118_D2', 'nNNPDF10_nlo_as_0118_Al27', 'nNNPDF10_nlo_as_0118_Au197',
            'EPPS16nlo_CT14nlo_Al27', 'EPPS16nlo_CT14nlo_Au197', '']
triggers = ['', 'trig']
projectile = [2212, 2112, 2212, 2112, 2212, 2212, 1000010020, 1000020030]
target = [2212, 2212, 2112, 2112, 1000130270, 1000791970, 1000791970, 1000791970]
opt = [1, 7, 8]
syst = 0
trigs = 0
pdf  = [2, 3, 7, 8, 3, 8, 3, 8]
pdfA = [0, 0, 0, 0, 4, 4, 1, 1]
option1 = ''
name = ''
for i in opt:
    option1 = option1 + options[i]
    name = name + names[i]
name = name + 'nPDFs_opt_pnPDF'
nEvents = 4e7

nPDFsetsA = ''
nPDFsetsB = ''
for i in pdfA:
    nPDFsetsA = nPDFsetsA + nPSFs[i]

for i in pdf:
    nPDFsetsB = nPDFsetsB + nPSFs[i]

file = 'main113_pA_nPDF'
file_out = 'SmallSystems/'+system[syst]+'/'+system[syst]+name+triggers[trigs]
os.system('make ' + file)

n_threads = 25
start = 941020
for i in range(0, n_threads):
    option = option1
    print(option)
    os.system(f'./{file} {start + i * 10000} {i} {n_threads} {nEvents} {file_out} {projectile[syst]} '
              f'{target[syst]} {trigs} {len(pdf)} {nPDFsetsB} {nPDFsetsA} {len(opt)} {option} &')

try:
    while True:
        sleep(10)
        print("yolo")
        ss = os.popen('ps -fu | grep '+file).read()
        if ss.find('./'+file) == -1:
            break
    os.system('rm ' + file_out + '.root')
except:
    print("first try "*10)
if n_threads > 1:
    os.system('hadd ' + file_out + '.root ' + file_out + '_*' ' && rm ' + file_out + '_*')

print("\a"*10)