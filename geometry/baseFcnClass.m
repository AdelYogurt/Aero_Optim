function C=baseFcnClass(U,N1,N2)
% default of class function
%
NP=calNormPar(N1,N2);
C=U.^N1.*(1-U).^N2./NP;
end

function nomlz_par=calNormPar(N1,N2)
% calculate normailize class function parameter by N1, N2
%
nomlz_par=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
nomlz_par((N1 == 0) & (N2 == 0))=1;
end
