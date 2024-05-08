function dcvdt = Odes_CV(tt, cv, tmod, munT)

munT  = interp1(tmod, munT, tt);

dcvdt = munT*cv;


end