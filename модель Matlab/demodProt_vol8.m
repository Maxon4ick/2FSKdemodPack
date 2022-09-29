%Считывание сигнала
fid = 'FSK2_8000_50Bd_7_5st5.wav';
f = fopen(fid,'rb');
fseek(f,44, 'bof');
s = fread(f,[1,inf],'short');
fclose(f);
%Выделение огибающей
sig = hilbert(s);
out2 = bandpass(sig,[1100,1300],8e3);
c2 = abs(out2);
out3 = bandpass(sig,[700,900],8e3);
c3 = abs(out3);
out = c2-c3;
cik = fix(length(out)/1200);
ro = zeros(1,1200*cik);
for i = 1:1:cik
nach = (i-1)*1200+1;
kon = i*1200;
strob = out(nach:kon);
delta = sum(strob)/1200;
strob = strob-delta;
ro(nach:kon) = strob;
end
tr = zeros(1,length(ro));
 for i = 1:1:length(ro)
 if(ro(i)<0)
     tr(i) = -1000;
 else
     tr(i) = 1000;
 end
 end
 %корреляция
 %Формирование опроного сигнала синхронизации
a1 = -ones(160,1)';
a3 = ones(240,1)';
a2 = zeros(800,1)';
a4 = zeros(801,1)';
op = [a1 a2 a3];
p1 = [a1 a4 a3];
po = [op op op p1];
en = length(tr)-4800;
strobOp = zeros(4801,1);
synx = zeros(en,1);
for i = 1:1:en
    strobOp = tr(i:i+4799);
    synx(i)= Cor(strobOp,po);
end
noc = en/1200;
k= zeros(1,noc);
point = zeros(1,noc);
for i = 1:1:noc
nach = (i-1)*1200+1;
kon = i*1200;
strobCo = synx(nach:kon);
k(i) = max(strobCo);
end
for m = 1:1:noc
nach = (m-1)*1200+1;
kon = m*1200;
strobMa = synx(nach:kon);
for i = 1 :1:1200
   if ((strobMa(i) == k(m))&& (strobMa(i) >10e5))
       point(m) = (m-1)*1200+ i;
   end
end
end
size = nnz(point);
trPoint = nonzeros(point);
su = 0;
sigDecod = zeros(1,5*size);
%Декодирование
for i = 1:1:size
strobKon =tr(trPoint(i):trPoint(i)+1200);
for m = 1:1:5
n = m*160+1;
b = (m+1)*160;
sn = (i-1)*5;
su = sum(strobKon(n:b));
if (su >0)
    sigDecod(sn+m) = 1000;
end
end
end
