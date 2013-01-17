function [var]=PhaseVar(realPhase,estPhase)
[M N] = size(realPhase);



dif = angle( exp( 1i*(estPhase-realPhase) ) );
piston = sum(dif(:))/(M*N);

estPhase = angle( exp( 1i*(estPhase-piston) ) );

dif = angle( exp( 1i*(estPhase-realPhase) ) );
var = (sum(sum(dif.^2))/(M*N));



end