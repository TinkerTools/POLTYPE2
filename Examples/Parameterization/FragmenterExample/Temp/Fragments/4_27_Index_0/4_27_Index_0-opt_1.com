%RWF=/scratch/liuchw/Gau-4_27_Index_0-farH_in17/,1267GB
%Nosave
%Chk=4_27_Index_0-opt_1.chk
%Mem=20GB
%Nproc=5
#P opt=(ModRedundant,maxcycles=400,Loose) MP2/6-31G* MaxDisk=1267GB 

4_27_Index_0 Gaussian OPT Calculation on node165.bme.utexas.edu

0 1
 C    7.862500   -0.649200    1.540900
 N    6.609300   -0.198000    0.919400
 C    5.417000   -0.912000    1.425500
 C    4.207900   -0.750700    0.463500
 C    6.490300    1.264400    0.984300
 H    5.638500   -1.984900    1.468800
 H    5.161100   -0.544900    2.425100
 O    4.378600   -0.195400   -0.626700
 H    6.572900    1.598000    2.027600
 H    5.517900    1.623500    0.630900
 H    8.111400   -1.635000    1.129800
 H    7.742600   -0.778300    2.624200
 H    8.640600    0.102000    1.340400
 H    3.214000   -1.125600    0.749000
 H    7.292000    1.663400    0.345600


