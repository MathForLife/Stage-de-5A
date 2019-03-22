function result=Operator(Hist,Var,Transpose)
    % Hist est la matrice de dimension [m*n,nbins] reresentant l'operateur A ou B
    % Si Transpose==True alors la variable est une des variables duales q2
    % ou q3
    % Au contraire, si Transpose==False alors la variable consideree est la
    % variable primale u
    if Transpose
        result=Hist'*Var;
    else
        result=Hist*Var;
    end