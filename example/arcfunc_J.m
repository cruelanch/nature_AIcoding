function j = arcfunc_J(x)
if x>=0 && x<=0.3646
    j=1.09542*x^2+0.214217*x+2.33727*sqrt(x);
elseif x<=0.9999
    j=-0.706692*log(0.386013*(1-x))+1.75017*x;
else
    j=100;
end
end

