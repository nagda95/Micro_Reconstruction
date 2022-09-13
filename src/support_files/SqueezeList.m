function Listen=SqueezeList(Liste)
Liste=sort(Liste);
Listen=zeros(size(Liste,1),1);
Listen(1,1)=Liste(1,1);
num=1;
for i=2:size(Liste,1)
    if Liste(i,1)==Liste(i-1,1)
        continue
    else
        num=num+1;
        Listen(num,1)=Liste(i,1);
    end
end
Listen=Listen(1:num,1);
    