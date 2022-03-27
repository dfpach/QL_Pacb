function state = convertvarstostateQ4(RARsent,PreambleTrans,Pacb)

format long

RARsentref=[0:1:15];

PreambleTransref=[0:1:54];

Pacbref=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1];



Pacb=floor(Pacb*100/5)*0.05;

%total=8*8*10*6*6;

total=16*55*20;

%find(abs(RARsentref-RARsent)<0.001)
%find(abs(PreambleTransref-PreambleTrans)<0.001)
%find(abs(Pacb-Pacbref)<0.001)

state=(find(abs(RARsentref-RARsent)<0.001)-1)*20*55+(find(abs(PreambleTransref-PreambleTrans)<0.001)-1)*20+find(abs(Pacb-Pacbref)<0.001);

