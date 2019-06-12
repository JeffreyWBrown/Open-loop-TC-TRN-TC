classdef Synapse_Depression < handle
    %Synapse demonstrating depression from Tsodyks and Markarm PNAS 1997
    %  Receiveing spike changes the amount of effective synapse from
    % E to E + U_se*R, where R is the amount of recovered synpases
    % E neurons decay to Inactive state I (I = 1 - R - E) and inactiave
    % reocer to R at rate tau_rec.
    properties
        Model = 'Depressing Synapse';
        E_s = 0            ; %potentials of synapse in mV
        U_se = .67          ;%  utilization of synaptic efficiency    %max probability of openning
        A_se = 500;      % nominal current in pA (used for testing)
        tau_inact = 3;       %time of inactivation
        tau_rec = 800;  % time scale of roveryy
        tau_mem = 50;   
        gs = (.25*10^-4)*10^6% %nominal conductance (used for simulation)
        
        R = 1;              %prob of recoverved synapses
        E = 0;             %prob of effective synapses 
        I = 0;             % prob of inactiave synapses
        Gs = 0;             %calculated conductance
        
    end
    
    methods
        function obj = Synapse(input)
            if(nargin > 0)
                obj.Model = input;
            end
        end
        
        function integrate(ss,dt)
            
            %This will integrate effictive values of synpases
            newval = timeint(0,ss.E,dt,ss.tau_inact);
            ss.E = newval;
            
            %update Recovered synapses
            newval  = timeint(1-ss.E,ss.R,dt,ss.tau_rec);
            ss.R = newval;
            
            ss.I = 1 - ss.R - ss.E;
            
            ss.Gs = ss.E*ss.gs;
           
        end
        
        function SpikeReceived(ss)
            %Utlilize U_se of recovered channels.
            ss.E = ss.E + ss.U_se*ss.R;
            ss.R = ss.R - ss.U_se*ss.R;
        end
            
    end
    
end

