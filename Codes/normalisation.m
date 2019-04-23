function u=normalisation(vect,a,b)
% Condition : a<=0 et b>=0
Min=min(vect);
Max=max(vect);
   
u=(a+b)*(vect-Min)/(Max-Min);
end