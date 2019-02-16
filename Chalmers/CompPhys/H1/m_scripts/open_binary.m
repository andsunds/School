clc;clf;


FILENAME = '../data/INIDATA_temp-700_pres-1.bin';
fID=fopen(FILENAME,'rb');
data1=fread(fID,[3,inf],'real*8').';
fclose(fID);

AMU = 1.0364e-4;
m_Al = 27*AMU;
kB= 8.61733e-5;
N_atoms=4^4;

T=sum(sum(data1.^2,2),1) / (3*m_Al*N_atoms*kB)
data2=load('../data/phase-space_temp-500_pres-1.tsv');
T=sum(sum(data2(:,4:end).^2,2),1) / (3*m_Al*N_atoms*kB)