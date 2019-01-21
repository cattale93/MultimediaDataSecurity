function [av,z]=bf_average_omega(t,w)

t=abs(t);

nonzero=find(t>0);

elements=zeros(size(t));

elements(nonzero)=t(nonzero).^(-1i*w/log(10));

% elements(nonzero)=exp(-j*2*pi*n*log10(t(nonzero)));


z=sum(elements)/length(t);

av=abs(z);

end

