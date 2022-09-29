function r = Cor(x,y) 
En = length(x);
r = 0;
for i = 1:1:En
   r = r+x(i)*y(i);
end