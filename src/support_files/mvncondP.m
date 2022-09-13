function [mubar, Sigmabar]=mvncondP(mu,Sigma,values)
    ind=numel(values);
    n=size(Sigma,1);
    mu1=mu((1:ind),1);
    mu2=mu(((ind+1):n),1);
    Sigma11=Sigma(1:ind,1:ind);
    Sigma12=Sigma(1:ind,(ind+1):n);
    Sigma21=Sigma((ind+1):n,1:ind);
    Sigma22=Sigma((ind+1):n,(ind+1):n);

    s21s11inv=Sigma21/Sigma11;

    mubar    = mu2 + s21s11inv*(values - mu1);
    Sigmabar = Sigma22 - s21s11inv*Sigma12;
end