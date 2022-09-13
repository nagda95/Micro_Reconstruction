function [Sigma]=Sp_Coeffs_PS(A,n)
t=1;
for l=1:n
    for m=0:l
        if m==0
            Sigma(t,t)=A(l);
            t=t+1;
            continue
        end        
        Sigma(t,t)=A(l)/2;
        t=t+1;
        Sigma(t,t)=A(l)/2;
        t=t+1;
    end
end