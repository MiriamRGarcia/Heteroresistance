function dzdt = Odes_MVCV(tt, zz, mug, muk, AA)

% Sizes of the problem:
nr      = size(muk, 1);
Nexp    = size(muk, 2);

% ODEs for the mean:
xx      = zz(1:Nexp*nr);
dxdt    = AA*xx;

% ODEs for variance:
varx    = zz((Nexp*nr + 1):Nexp*(nr + 1));
xx_aux  = reshape(xx, nr, Nexp);
xT      = sum(xx_aux, 1);
munT    = (sum(xx_aux.*(mug - muk), 1)./xT).';

dvarxdt = 2*diag(munT)*varx + sum(xx_aux.*(mug + muk), 1).';

% ODEs for covariance:
covx    = zz((Nexp*(nr + 1) + 1):Nexp*(nr + 2));

dcovxdt = diag(munT)*covx;

dzdt = [dxdt;
        dvarxdt;
        dcovxdt];
    
end