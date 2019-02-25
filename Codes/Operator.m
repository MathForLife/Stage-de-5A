function result=Operator(g,u,T)
    
    sz=size(g);
    if T
        g=reshape(g,[],sz(end));
        g=g*u;
        sz(end)=[];
        result=reshape(g,sz);
    else
        result=sum(sum(u.*g,1),2);
        result=reshape(result,sz(end),1);
    end