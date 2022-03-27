function state = convertvarstostateQ8(PreambleTransM,PreambleTransCV,DeltaNps,Pacb)

format long



PreambleTransMref=[0:1:29]; %30 estados
PreambleTransCVref=[0:0.2:1]; %5 estados
DeltaNpsref=[1:1:3];%1 crecio,2 disminuyo, 3 igual : 3 estados
Pacbref=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]; %16 estados



%Pacb=floor(Pacb*100/5)*0.05;

%total=8*8*10*6*6;

total=30*5*3*16;

%find(abs(PreambleTransMref-PreambleTransM)<0.001)-1
%find(abs(PreambleTransCV-PreambleTransCVref)<0.001)
%find(abs(DeltaNps-DeltaNpsref)<0.001)
%find(abs(Pacb-Pacbref)<0.001)

state=(find(abs(PreambleTransMref-PreambleTransM)<0.001)-1)*16*3*5+(find(abs(PreambleTransCV-PreambleTransCVref)<0.001)-1)*16*3+(find(abs(DeltaNps-DeltaNpsref)<0.001)-1)*16+find(abs(Pacb-Pacbref)<0.001);

