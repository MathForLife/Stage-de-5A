function Hn=Heavyside_eta(z,eta,n)
    if n==1
        Hn=zeros(size(z));
        
        cond=abs(z)<=eta;
        
        Hn(z>eta)=1;
        Hn(z<-eta)=0;
        Hn(cond)=0.5*(1+z(cond)/eta+sin(pi*z(cond)/eta)/pi);
    elseif n==2
        Hn=0.5*(1+2/pi*atan(z/eta));
    end
end