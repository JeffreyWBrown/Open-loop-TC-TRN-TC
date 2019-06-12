function [ spike ] = spike_check(X,SpikeTh)
%This function will return a value of "True" if a spike is measured
%
%X is the last three voltages
%SpikeTh is the spike threshold
%
%Spikes will be defined when the middle value of the set Xin, is greater
%than SpikeTh along with the requirement that
%X(2) - X(1) > 0 and X(2) - X(3) < 0, where 
%X(3) is voltage at current time,
%X(2) is voltage at one time step backward, and
%X(1) is voltage at 2 steps back


spike = (X(2) >= SpikeTh) & (X(1) - X(2) < 0) & (X(2) - X(3) > 0);
    

