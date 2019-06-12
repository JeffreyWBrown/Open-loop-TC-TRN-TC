%question 3 of chapter 5
%integrate and fire model


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function integrate using 2nd order scheme of Abbott and Dayan (5.48)
%
%V(t + dt) = timeint(V(t),V'(t), dt, tauv)
%
%V is function value
%V' is first time derivative
%dt is time step
%tauv is timescale of circuit
%Vinf is time value of inf

function [Vout] = timeint(Vinf,x, dt, tau)
  Vout = Vinf + (x - Vinf).*exp(-dt./tau);
end
