%RWF=/scratch/liuchw/Gau-32_40_Index_0-farH_in17/,1267GB
%Nosave
%Chk=32_40_Index_0-opt_1.chk
%Mem=20GB
%Nproc=5
#P opt=(ModRedundant,maxcycles=400,Loose) MP2/6-31G* MaxDisk=1267GB 

32_40_Index_0 Gaussian OPT Calculation on node165.bme.utexas.edu

1 1
 H   10.925500    1.350500   -0.390400
 C   10.063400    2.007500   -0.249000
 N    8.877000    1.153900    0.067600
 H    8.771100    0.509800   -0.732200
 H    9.852300    2.551900   -1.173100
 H   10.231500    2.697000    0.582500
 C    7.579700    1.964000    0.134000
 C    9.073600    0.297700    1.305400
 H    7.817300    2.944100    0.561500
 H    7.252900    2.084900   -0.904200
 H    9.994300   -0.278300    1.166400
 H    9.195000    0.989200    2.147000
 H    6.802300    1.480900    0.744100
 H    8.231700   -0.387300    1.484200


