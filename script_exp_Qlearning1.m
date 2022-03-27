Q = zeros((16*55*20),20);

epsilon=0.99;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov1celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.9;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov2celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.9;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov3celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.9;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov4celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end


epsilon=0.9;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov5celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.9;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov6celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.9;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov7celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.8;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov8celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.8;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov9celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.8;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov10celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.8;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov11celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.8;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov12celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.8;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov13celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.8;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov14celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end

epsilon=0.8;

for i=1:144

rng('shuffle')
archivotraficoH2H = textread('/Users/diego/Documents/ECBRACH/datoscelda5161/nov15celda5161.txt','%s');

traficoH2H(i)=str2double(cell2mat(archivotraficoH2H(i)));

traficoH2H(i)=traficoH2H(i)*4.1581;

[averagePerRAO, avPreamStatsPerRAO,Q,vectorPacb] = LTEA_M_H_ACB_QL(1e4.*betarnd(3,4,30000,1), unifrnd(0,10*60*1000,floor(traficoH2H(i)),1),Q,epsilon);

end