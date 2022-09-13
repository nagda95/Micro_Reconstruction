function [Y,first_coeff]=Yharmonic_coeffs(n,phi,theta)
    dirs=[phi theta];
    Y=zeros(1,(n+1)^2-1);
    basisType='complex';
    Y_N = getSH(n, dirs, basisType);
    x=3;
    for i=1:n
        Y(1,i^2)=real(Y_N(1,x));
        for j=1:i
            Y(1,i^2+2*(j-1)+1)=2*real(Y_N(1,x+j));
            Y(1,i^2+2*(j-1)+2)=2*imag(Y_N(1,x+j));
        end
        x=x+2*(i+1);
    end
    first_coeff=real(Y_N(1,1));
end
