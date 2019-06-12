classdef Synapse < handle
    %Synapse will create a synapse that will
    %return conductances
    
    properties
        Model = 'Saturating Synapse';
        E_s = 0    ;        ; %potentials of synapse in mV
        Pmax = 1.%.50;         %max probability of openning
        gs = 1.2*10^3;      % nominal conductance in ms/mm^2
        tau_s = 5.6;          %time in ms (GABAa)
        
        z = 0;              %gating variable initially at 0
        Ps = 0;             %initial open probability 
        Gs = 0;             %calculated conductance.
        
    end
    
    methods
        function obj = Synapse(input)
            if(nargin > 0)
                obj.Model = input;
            end
        end
        
        function integrate(ss,dt)
            
            %This will integrate the synapse in time and determine
            %conductance
            Pinf = exp(1)*ss.Pmax*ss.z/(exp(1)*ss.Pmax*ss.z + 1);
            tau_i = ss.tau_s/(exp(1)*ss.Pmax*ss.z + 1.);
            
            %update Ps
            newval = timeint(Pinf,ss.Ps,dt,tau_i);
            ss.Ps = newval;
            
            %update z
            newval = ss.z*exp(-dt/ss.tau_s);
            ss.z = newval;
            
            %calculated conductance
            ss.Gs = ss.Ps*ss.gs;
        end
        
        function SpikeReceived(ss)
            %This will update z variable to open channel.
            ss.z = 1.;
        end
            
    end
    
end

