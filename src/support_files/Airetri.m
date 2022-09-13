function [sec]=Airetri(a,b,c)
% Calcul de la surface d'un triangle

aa = max([a, b, c]);
cc = min([a, b, c]);
if aa ~= b & cc ~= b
bb = b;
elseif aa ~= a & cc ~= a
bb = a;
else
bb = c;
end

if b + c ~= a
sec = 0.25 * ((aa + bb + cc) * (-aa + bb + cc) * (aa - bb + cc) * (aa + bb - cc)) ^ (1 / 2);
else
sec = 0;
end