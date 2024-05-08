function dcvdt = Odes_SensCV(tt, ss, tmod, munT, sens_munT)

np        = size(sens_munT, 2);

munT      = interp1(tmod, munT, tt);
sens_munT = interp1(tmod, sens_munT, tt).';

cv  = ss(1);
scv = ss(2:(np + 1));

dcvdt = [munT*cv;
         cv*sens_munT + munT*scv];

end