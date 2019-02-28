function dn=delta_eta(z,eta,n)
    if n==1
        dn=zeros(size(z));
        
        cond=abs(z)<=eta;
        
        dn(cond)=0.5/eta*(1+cos(pi*z(cond)/eta));
    elseif n==2
        dn=(eta*pi*(1+(z/eta).^2)).^-1;
    end
end
    