function Color=Translate_Color(Value,Min,Max)
colormap(jet(64));
cmap=colormap;
Scale=transpose(Min:(Max-Min)/63:Max);
if Value<Min
    Color=[0,0,0.5625];
elseif Value>Max
    Color=[0.5,0,0];
else
   Color=interp1(Scale,cmap,Value);
end