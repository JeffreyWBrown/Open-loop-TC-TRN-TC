%Initial Parameters
SimTime =1500;          %Simulation time in ms
dt = .1 ;               %Time step in ms
Nt = SimTime/dt;        %Number of steps in simulation.
SpkTh = 00;             %Spike Threshold           
SimNum = 1000;          %Number of simulations at fixed frequencies
Vinit = -66;            %Initial Membrane Potential
EqlbTime = 200;         %Time in ms to equilibrate and not gather statistics

%Define Synapses
%Excitatory external synapse TC_A
SynInE_A = Synapse_Depression;  
SynInE_A.E_s = 0;
SynInE_A.gs = 32;
SynInE_A.U_se = 0.76;
SynInE_A.tau_inact = 2.64;
SynInE_A.tau_rec = 125.0;

% %Excitatory external synapse TC_B
SynInE_B = Synapse_Depression; 
SynInE_B.E_s = 0;
SynInE_B.gs = 32;
SynInE_B.U_se = 0.76;
SynInE_B.tau_inact = 2.64;
SynInE_B.tau_rec = 125.0;

% %Excitatory external synapse TC_C
SynInE_C = Synapse_Depression;  
SynInE_C.E_s = 0;
SynInE_C.gs = 32;
SynInE_C.U_se = 0.76;
SynInE_C.tau_inact = 2.64;
SynInE_C.tau_rec = 125.0;

%Inhibitory synapse from TRN_A to TC_A 
SynTRN_A_TC_A = Synapse_Depression;
SynTRN_A_TC_A.E_s = -80;
SynTRN_A_TC_A.U_se = 0.62;
SynTRN_A_TC_A.tau_inact = 16.62;
SynTRN_A_TC_A.tau_rec = 167.29;

%Inhibitory synapse from TRN_A to TC_B 
SynTRN_A_TC_B = Synapse_Depression;
SynTRN_A_TC_B.E_s = -80;
SynTRN_A_TC_B.U_se = 0.62;
SynTRN_A_TC_B.tau_inact = 16.62;
SynTRN_A_TC_B.tau_rec = 167.29;

%Inhibitory synapse from TRN_B to TC_B 
SynTRN_B_TC_B = Synapse_Depression;
SynTRN_B_TC_B.E_s = -80;
SynTRN_B_TC_B.U_se = 0.62;
SynTRN_B_TC_B.tau_inact = 16.62;
SynTRN_B_TC_B.tau_rec = 167.29;

%Inhibitory synapse from TRN_B to TC_C 
SynTRN_B_TC_C = Synapse_Depression;
SynTRN_B_TC_C.E_s = -80;
SynTRN_B_TC_C.U_se = 0.62;
SynTRN_B_TC_C.tau_inact = 16.62;
SynTRN_B_TC_C.tau_rec = 167.29;

%Inhibitory synapse from TRN_C to TC_C 
SynTRN_C_TC_C = Synapse_Depression;
SynTRN_C_TC_C.E_s = -80;
SynTRN_C_TC_C.U_se = 0.62;
SynTRN_C_TC_C.tau_inact = 16.62;
SynTRN_C_TC_C.tau_rec = 167.29;

%Inhibitory synapse from TRN_A to TRN_B 
SynTRN_A_TRN_B = Synapse_Depression;
SynTRN_A_TRN_B.E_s = -75;
SynTRN_A_TRN_B.U_se = 0.62;
SynTRN_A_TRN_B.tau_inact = 15;
SynTRN_A_TRN_B.tau_rec = 225;

%Inhibitory synapse from TRN_A to TRN_C 
SynTRN_A_TRN_C = Synapse_Depression;
SynTRN_A_TRN_C.E_s = -75;
SynTRN_A_TRN_C.U_se = 0.62;
SynTRN_A_TRN_C.tau_inact = 15;
SynTRN_A_TRN_C.tau_rec = 225;

%Inhibitory synapse from TRN_B to TRN_A 
SynTRN_B_TRN_A = Synapse_Depression;
SynTRN_B_TRN_A.E_s = -75;
SynTRN_B_TRN_A.U_se = 0.62;
SynTRN_B_TRN_A.tau_inact = 15;
SynTRN_B_TRN_A.tau_rec = 225;

%Inhibitory synapse from TRN_B to TRN_C 
SynTRN_B_TRN_C = Synapse_Depression;
SynTRN_B_TRN_C.E_s = -75;
SynTRN_B_TRN_C.U_se = 0.62;
SynTRN_B_TRN_C.tau_inact = 15;
SynTRN_B_TRN_C.tau_rec = 225;

%Inhibitory synapse from TRN_C to TRN_A 
SynTRN_C_TRN_A = Synapse_Depression;
SynTRN_C_TRN_A.E_s = -75;
SynTRN_C_TRN_A.U_se = 0.62;
SynTRN_C_TRN_A.tau_inact = 15;
SynTRN_C_TRN_A.tau_rec = 225;

%Inhibitory synapse from TRN_C to TRN_B 
SynTRN_C_TRN_B = Synapse_Depression;
SynTRN_C_TRN_B.E_s = -75;
SynTRN_C_TRN_B.U_se = 0.62;
SynTRN_C_TRN_B.tau_inact = 15;
SynTRN_C_TRN_B.tau_rec = 225;

%Excitatory Synapse from TC_A to TRN_A 
SynTC_A_TRN_A = Synapse_Depression;
SynTC_A_TRN_A.E_s = 0;
SynTC_A_TRN_A.gs =  150;
SynTC_A_TRN_A.U_se = 0.76;
SynTC_A_TRN_A.tau_inact = 2.64;
SynTC_A_TRN_A.tau_rec = 500;

%Excitatory Synapse from TC_B to TRN_B 
SynTC_B_TRN_B = Synapse_Depression;
SynTC_B_TRN_B.E_s = 0;
SynTC_B_TRN_B.gs =  150;
SynTC_B_TRN_B.U_se = 0.76;
SynTC_B_TRN_B.tau_inact = 2.64;
SynTC_B_TRN_B.tau_rec = 500;

%Excitatory Synapse from TC_C to TRN_C
SynTC_C_TRN_C = Synapse_Depression;
SynTC_C_TRN_C.E_s = 0;
SynTC_C_TRN_C.gs =  150;
SynTC_C_TRN_C.U_se = 0.76;
SynTC_C_TRN_C.tau_inact = 2.64;
SynTC_C_TRN_C.tau_rec = 500;

%Excitatory Synapse (Depressing) TC_A to C_A
SynTC_A_C_A = Synapse_Depression;
SynTC_A_C_A.E_s = 0;
SynTC_A_C_A.gs = 50;
SynTC_A_C_A.U_se = 0.8113;
SynTC_A_C_A.tau_inact = 3;
SynTC_A_C_A.tau_rec = 160;

%Excitatory Synapse (Depressing) TC_B to C_B
SynTC_B_C_B = Synapse_Depression;
SynTC_B_C_B.E_s = 0;
SynTC_B_C_B.gs = 50;
SynTC_B_C_B.U_se = 0.8113;
SynTC_B_C_B.tau_inact = 3;
SynTC_B_C_B.tau_rec = 160;

%Excitatory Synapse (Depressing) TC_C to C_C
SynTC_C_C_C = Synapse_Depression;
SynTC_C_C_C.E_s = 0;
SynTC_C_C_C.gs = 50;
SynTC_C_C_C.U_se = 0.8113;
SynTC_C_C_C.tau_inact = 3;
SynTC_C_C_C.tau_rec = 160;

%CreateNeurons
TC_A = Neuron_Thalamic_Relay_deleuze2012; 

TRN_A = Neuron_TRN;  

C_A = Neuron_Regular_Spiking; 

TC_B = Neuron_Thalamic_Relay_deleuze2012; 

TRN_B = Neuron_TRN;  

C_B = Neuron_Regular_Spiking;  

TC_C = Neuron_Thalamic_Relay_deleuze2012; 

TRN_C = Neuron_TRN; 

C_C = Neuron_Regular_Spiking; 

%Define and Loop Over Full Homogeneous Synaptic Parameter Space

RRg_values=[0 50 100 150 200 250 300 350 400 450];
CC_values=[0 0.06 0.12 0.18 0.24 0.3 0.36];
Op_values=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

RRg_length=length(RRg_values);
CC_length=length(CC_values);
Op_length=length(Op_values);

Spike_C_A_Times = zeros(RRg_length,CC_length,Op_length,SimTime,SimNum);
Spike_C_B_Times = zeros(RRg_length,CC_length,Op_length,SimTime,SimNum);
Spike_C_C_Times = zeros(RRg_length,CC_length,Op_length,SimTime,SimNum);
RRg_index=0;
CC_index=0;
Op_index=0;

for m = 1:RRg_length
    RRg_index = RRg_index + 1;
    RRg = RRg_values(m);
    SynTRN_A_TRN_B.gs =  0.9956*RRg;
    SynTRN_A_TRN_C.gs =  0.9824*RRg;
    SynTRN_B_TRN_A.gs =  0.9956*RRg;
    SynTRN_B_TRN_C.gs =  0.9956*RRg;
    SynTRN_C_TRN_A.gs =  0.9824*RRg;
    SynTRN_C_TRN_B.gs =  0.9956*RRg;
    
    for o = 1:CC_length
       CC_index = CC_index + 1;
       CC = CC_values(o);
        
       for p = 1:Op_length
           Op_index = Op_index + 1;
            Op = Op_values(p);
            SynTRN_A_TC_A.gs =  (1 - Op)*80;
            SynTRN_B_TC_B.gs =  (1 - Op)*80;
            SynTRN_C_TC_C.gs =  (1 - Op)*80;
            SynTRN_A_TC_B.gs =  Op*80;
            SynTRN_B_TC_C.gs =  Op*80;
           
%Initialize Membrane Potentials
V_C_A = zeros(Nt,1);     
V_C_A(1) = Vinit;
V_TC_A = zeros(Nt,1);
V_TC_A(1) = Vinit;
V_TRN_A = zeros(Nt,1);
V_TRN_A(1) = Vinit;

V_C_B = zeros(Nt,1);     
V_C_B(1) = Vinit;
V_TC_B = zeros(Nt,1);
V_TC_B(1) = Vinit;
V_TRN_B = zeros(Nt,1);
V_TRN_B(1) = Vinit;

V_C_C = zeros(Nt,1);     
V_C_C(1) = Vinit;
V_TC_C = zeros(Nt,1);
V_TC_C(1) = Vinit;
V_TRN_C = zeros(Nt,1);
V_TRN_C(1) = Vinit;
    
inputE_A = zeros(Nt,1); %External Synapse to TC_A
inputE_B = zeros(Nt,1); %External Synapse to TC_B
inputE_C = zeros(Nt,1) ;%External Synapse to TC_C
time = zeros(Nt,1);   %Time Vector for Histogram

%Generate Inputs

run=0;
for l = 1:SimNum
TC_A_rate=0.01;
TC_B_rate=0.04;
TC_C_rate=0.04;

run=run+1;
status=[m o p run]
   
        spikein_A =0;             
        spikeout_A = 0;          
        n1Aspikeout = 0;         
        spikein_B =0;            
        spikeout_B = 0;           
        n1Bspikeout = 0;        
        spikein_C =0;             
        spikeout_C = 0;          
        n1Cspikeout = 0;         

     inputE_A(:) = -60;
     inputE_B(:) = -60;
     inputE_C(:) = -60;


%Reset Potentials
   
    TC_A.V = Vinit;
    C_A.V = Vinit;
    TRN_A.V = Vinit;

    TC_B.V = Vinit;
    C_B.V = Vinit;
    TRN_B.V = Vinit;

    TC_C.V = Vinit;
    C_C.V = Vinit;
    TRN_C.V = Vinit;
    
%Reset Synapses to Closed
   
    SynInE_A.E = 0.0;
    
    SynInE_B.E = 0.0;
    
    SynInE_C.E = 0.0;
    
    SynTC_A_C_A.E = 0.0;

    SynTC_B_C_B.E = 0.0;
 
    SynTC_C_C_C.E = 0.0;
    
    SynTRN_A_TC_A.E = 0.0;
    
    SynTRN_A_TC_B.E = 0.0;
     
    SynTRN_B_TC_B.E = 0.0;
    
    SynTRN_B_TC_C.E = 0.0;
    
    SynTRN_C_TC_C.E = 0.0;

    SynTRN_A_TRN_B.E = 0.0;
    
    SynTRN_A_TRN_C.E = 0.0;

    SynTRN_B_TRN_A.E = 0.0;
    
    SynTRN_B_TRN_C.E = 0.0;
    
    SynTRN_C_TRN_A.E = 0.0;
    
    SynTRN_C_TRN_B.E = 0.0;

    SynTC_A_TRN_A.E = 0.0; 

    SynTC_B_TRN_B.E = 0.0; 
        
    SynTC_C_TRN_C.E = 0.0; 
    
    %Input Spike Generation
     spkE_Acnt = 1; 
     spikein_A = 0;  
     spikeC_A = 0;  
     spikeout_A = 0;
     
     spkE_Bcnt = 1; 
     spikein_B = 0; 
     spikeC_B = 0;
     spikeout_B = 0;
     
     spkE_Ccnt = 1; 
     spikein_C = 0;   
     spikeC_C = 0;  
     spikeout_C = 0;
     
%Reset Inputs to -60 mV and Generate External Input to TC layer
     inputE_A(:) = -60;
     inputE_B(:) = -60;
     inputE_C(:) = -60;
     
     onsetE_A = EqlbTime + 100*rand(1);
     onsetE_B = EqlbTime + 100*rand(1);
     onsetE_C = EqlbTime + 100*rand(1);
     
    [nspikeE_A,tempspikeE_A] = poisson(TC_A_rate,onsetE_A,400);
    tempspikeE_A=[tempspikeE_A 400:5:1495];
    nspikeE_A=length(tempspikeE_A);
     
    [nspikeE_B,tempspikeE_B] = poisson(TC_B_rate,onsetE_B,SimTime);
     
    [nspikeE_C,tempspikeE_C] = poisson(TC_C_rate,onsetE_C,SimTime);
    
     for i = 1:Nt 
         time(i) = (i-1)*dt;
        
         if spkE_Acnt <= nspikeE_A   
             if tempspikeE_A(spkE_Acnt) < (i+1)*dt
                 inputE_A(i+1) = 30;
                 spkE_Acnt = spkE_Acnt + 1;
                 spikein_A = spikein_A + 1;
             end
         end
         
          if spkE_Bcnt <= nspikeE_B  
             if tempspikeE_B(spkE_Bcnt) < (i+1)*dt
                 inputE_B(i+1) = 30;
                 spkE_Bcnt = spkE_Bcnt + 1;
                 spikein_B = spikein_B + 1;
             end
          end
         
           if spkE_Ccnt <= nspikeE_C  
             if tempspikeE_C(spkE_Ccnt) < (i+1)*dt
                 inputE_C(i+1) = 30;
                 spkE_Ccnt = spkE_Ccnt + 1;
                 spikein_C = spikein_C + 1;
             end
           end
     end

%Integration     

for i = 2:Nt
    
%TRN-TRN Gap Junctions
     TRN_ig=22;
     
    GJTRN_A_TRN_B_Gs=0.929/(1/(TRN_ig*CC) - 1/TRN_ig);
    GJTRN_A_TRN_B_Es=TRN_A.V;
    
    GJTRN_A_TRN_C_Gs=0.745/(1/(TRN_ig*CC) - 1/TRN_ig);
    GJTRN_A_TRN_C_Es=TRN_A.V;
    
    GJTRN_B_TRN_A_Gs=0.929/(1/(TRN_ig*CC) - 1/TRN_ig);
    GJTRN_B_TRN_A_Es=TRN_B.V;
    
    GJTRN_B_TRN_C_Gs=0.929/(1/(TRN_ig*CC) - 1/TRN_ig);
    GJTRN_B_TRN_C_Es=TRN_B.V;
    
    GJTRN_C_TRN_A_Gs=0.745/(1/(TRN_ig*CC) - 1/TRN_ig);
    GJTRN_C_TRN_A_Es=TRN_C.V;
    
    GJTRN_C_TRN_B_Gs=0.929/(1/(TRN_ig*CC) - 1/TRN_ig);
    GJTRN_C_TRN_B_Es=TRN_C.V;
   
%Chemical Synapses
    
     EsGsTC_A = ...
         SynInE_A.Gs*SynInE_A.E_s + ...   
         SynTRN_A_TC_A.Gs*SynTRN_A_TC_A.E_s;
        
     EsGsTC_B = ...
         SynInE_B.Gs*SynInE_B.E_s + ...   
         SynTRN_A_TC_B.Gs*SynTRN_A_TC_B.E_s + ...  
         SynTRN_B_TC_B.Gs*SynTRN_B_TC_B.E_s;

     EsGsTC_C = ...
         SynInE_C.Gs*SynInE_C.E_s + ...  
         SynTRN_B_TC_C.Gs*SynTRN_B_TC_C.E_s + ... 
         SynTRN_C_TC_C.Gs*SynTRN_C_TC_C.E_s;

     EsGsC_A = SynTC_A_C_A.Gs*SynTC_A_C_A.E_s;    
        
     EsGsC_B = SynTC_B_C_B.Gs*SynTC_B_C_B.E_s;   
        
     EsGsC_C = SynTC_C_C_C.Gs*SynTC_C_C_C.E_s;     
 
    EsGsTRN_A = SynTC_A_TRN_A.Gs*SynTC_A_TRN_A.E_s + ...  
                 SynTRN_B_TRN_A.Gs*SynTRN_B_TRN_A.E_s + ...
                 SynTRN_C_TRN_A.Gs*SynTRN_C_TRN_A.E_s + ...
                 GJTRN_B_TRN_A_Gs*GJTRN_B_TRN_A_Es + ...
                 GJTRN_C_TRN_A_Gs*GJTRN_C_TRN_A_Es;
        
     EsGsTRN_B = SynTC_B_TRN_B.Gs*SynTC_B_TRN_B.E_s + ...  
                 SynTRN_A_TRN_B.Gs*SynTRN_A_TRN_B.E_s + ...
                 SynTRN_C_TRN_B.Gs*SynTRN_C_TRN_B.E_s + ...
                 GJTRN_A_TRN_B_Gs*GJTRN_A_TRN_B_Es + ...
                 GJTRN_C_TRN_B_Gs*GJTRN_C_TRN_B_Es;
             

     EsGsTRN_C = SynTC_C_TRN_C.Gs*SynTC_C_TRN_C.E_s + ... 
                 SynTRN_A_TRN_C.Gs*SynTRN_A_TRN_C.E_s + ...
                 SynTRN_B_TRN_C.Gs*SynTRN_B_TRN_C.E_s + ...
                 GJTRN_A_TRN_C_Gs*GJTRN_A_TRN_C_Es + ...
                 GJTRN_B_TRN_C_Gs*GJTRN_B_TRN_C_Es;
    
     GsTC_A = SynInE_A.Gs + SynTRN_A_TC_A.Gs; 
     GsC_A = SynTC_A_C_A.Gs;
     GsTRN_A = SynTC_A_TRN_A.Gs + SynTRN_B_TRN_A.Gs + SynTRN_C_TRN_A.Gs + ...
         GJTRN_B_TRN_A_Gs + GJTRN_C_TRN_A_Gs;
     
     GsTC_B = SynInE_B.Gs + SynTRN_A_TC_B.Gs + SynTRN_B_TC_B.Gs;
     GsC_B = SynTC_B_C_B.Gs;
     GsTRN_B = SynTC_B_TRN_B.Gs + SynTRN_A_TRN_B.Gs + SynTRN_C_TRN_B.Gs + ...
         GJTRN_A_TRN_B_Gs + GJTRN_C_TRN_B_Gs;


     GsTC_C = SynInE_C.Gs + SynTRN_B_TC_C.Gs + SynTRN_C_TC_C.Gs;
     GsC_C = SynTC_C_C_C.Gs;
     GsTRN_C = SynTC_C_TRN_C.Gs + SynTRN_A_TRN_C.Gs + SynTRN_B_TRN_C.Gs + ...
         GJTRN_A_TRN_C_Gs + GJTRN_B_TRN_C_Gs;
     
     TC_A.integrate(dt,EsGsTC_A,GsTC_A);
     TRN_A.integrate(dt,EsGsTRN_A,GsTRN_A);
     C_A.integrate(dt,EsGsC_A,GsC_A);
     
     TC_B.integrate(dt,EsGsTC_B,GsTC_B);
     TRN_B.integrate(dt,EsGsTRN_B,GsTRN_B);
     C_B.integrate(dt,EsGsC_B,GsC_B);
    
     TC_C.integrate(dt,EsGsTC_C,GsTC_C);
     TRN_C.integrate(dt,EsGsTRN_C,GsTRN_C);
     C_C.integrate(dt,EsGsC_C,GsC_C);
    
    SynInE_A.integrate(dt);
    SynTRN_A_TC_B.integrate(dt);
    SynTRN_A_TC_A.integrate(dt);
    SynTRN_A_TRN_B.integrate(dt);
    SynTRN_A_TRN_C.integrate(dt);
    SynTC_A_TRN_A.integrate(dt);
    SynTC_A_C_A.integrate(dt);    
    SynInE_B.integrate(dt);
    SynTRN_B_TC_C.integrate(dt);
    SynTRN_B_TC_B.integrate(dt);
    SynTRN_B_TRN_A.integrate(dt);
    SynTRN_B_TRN_C.integrate(dt);
    SynTC_B_TRN_B.integrate(dt); 
    SynTC_B_C_B.integrate(dt);
    SynInE_C.integrate(dt);
    SynTRN_C_TC_C.integrate(dt);
    SynTRN_C_TRN_A.integrate(dt);
    SynTRN_C_TRN_B.integrate(dt);
    SynTC_C_TRN_C.integrate(dt); 
    SynTC_C_C_C.integrate(dt);

     V_TC_A(i) = TC_A.V;
     V_TRN_A(i) = TRN_A.V;
     V_C_A(i) = C_A.V;
     
     V_TC_B(i) = TC_B.V;
     V_TRN_B(i) = TRN_B.V;
     V_C_B(i) = C_B.V;

     V_TC_C(i) = TC_C.V;
     V_TRN_C(i) = TRN_C.V;
     V_C_C(i) = C_C.V;
    
    
%Process Spiking

 delay=10;
 if i > 10
     
     if spike_check(inputE_A(i-delay:i),SpkTh) 
         SynInE_A.SpikeReceived;
     end
     
     if spike_check(inputE_B(i-delay:i),SpkTh)  
         SynInE_B.SpikeReceived;   
     end
     
     if spike_check(inputE_C(i-delay:i),SpkTh) 
         SynInE_C.SpikeReceived;   
     end

     if spike_check(V_TC_A(i-delay:i),SpkTh)    
         SynTC_A_C_A.SpikeReceived;
         SynTC_A_TRN_A.SpikeReceived;  
         spikeout_A = spikeout_A + 1;
     end

     if spike_check(V_TC_B(i-delay:i),SpkTh)  
         SynTC_B_C_B.SpikeReceived;
         SynTC_B_TRN_B.SpikeReceived;  
         spikeout_B = spikeout_B + 1;
     end

     if spike_check(V_TC_C(i-delay:i),SpkTh) 
         SynTC_C_C_C.SpikeReceived;
         SynTC_C_TRN_C.SpikeReceived;  
         spikeout_C = spikeout_C + 1;
     end
     
    if spike_check(V_TRN_A(i-delay:i),SpkTh)   
        SynTRN_A_TC_B.SpikeReceived;
        SynTRN_A_TC_A.SpikeReceived;
        SynTRN_A_TRN_B.SpikeReceived;
        SynTRN_A_TRN_C.SpikeReceived;
    end

    if spike_check(V_TRN_B(i-delay:i),SpkTh)   
        SynTRN_B_TC_B.SpikeReceived;
        SynTRN_B_TC_C.SpikeReceived;
        SynTRN_B_TRN_A.SpikeReceived;
        SynTRN_B_TRN_C.SpikeReceived;
    end
    
    if spike_check(V_TRN_C(i-delay:i),SpkTh)   
        SynTRN_C_TC_C.SpikeReceived;
        SynTRN_C_TRN_A.SpikeReceived;
        SynTRN_C_TRN_B.SpikeReceived;
    end 
     
    if spike_check(V_C_A(i-delay:i),SpkTh)   
        spikeC_A = spikeC_A +1;
        Spike_C_A_Times(m,o,p,round(i*dt),l) = 1;
    end

    if spike_check(V_C_B(i-delay:i),SpkTh) 
        spikeC_B = spikeC_B +1;
        Spike_C_B_Times(m,o,p,round(i*dt),l) = 1;
    end   

    if spike_check(V_C_C(i-delay:i),SpkTh)   
        spikeC_C = spikeC_C +1;
          Spike_C_C_Times(m,o,p,round(i*dt),l) = 1;
    end
      end
 end %time stepping

end
       end
    end
end
%% 

%Generate L4 Spike Histograms
n=10; %Bin Size

Spike_C_A_hold=zeros(SimTime,SimNum);
Spike_C_A_hishold=zeros(SimTime,1);
Spike_C_A_hisraw=zeros(RRg_length,CC_length,Op_length,SimTime);
Spike_C_A_his=zeros(RRg_length,CC_length,Op_length,(SimTime/n));


for i=1:RRg_length
    for j=1:CC_length
        for k=1:Op_length           
            Spike_C_A_hold=squeeze(Spike_C_A_Times(i,j,k,:,:));
            Spike_C_A_hisraw(i,j,k,:)=Spike_C_A_hold*ones(SimNum,1);
            Spike_C_A_hishold=squeeze(Spike_C_A_hisraw(i,j,k,:));
            x=reshape(Spike_C_A_hishold,n,[]);
            Spike_C_A_hishold=sum(x);
            Spike_C_A_his(i,j,k,:)=Spike_C_A_hishold;
        end
    end
end

Spike_C_B_hold=zeros(SimTime,SimNum);
Spike_C_B_hishold=zeros(SimTime,1);
Spike_C_B_hisraw=zeros(RRg_length,CC_length,Op_length,SimTime);
Spike_C_B_his=zeros(RRg_length,CC_length,Op_length,(SimTime/n));

for i=1:RRg_length
    for j=1:CC_length
        for k=1:Op_length
            Spike_C_B_hold=squeeze(Spike_C_B_Times(i,j,k,:,:));
            Spike_C_B_hisraw(i,j,k,:)=Spike_C_B_hold*ones(SimNum,1);
            Spike_C_B_hishold=squeeze(Spike_C_B_hisraw(i,j,k,:));
            x=reshape(Spike_C_B_hishold,n,[]);
            Spike_C_B_hishold=sum(x);
            Spike_C_B_his(i,j,k,:)=Spike_C_B_hishold;
        end
    end
end

Spike_C_C_hold=zeros(SimTime,SimNum);
Spike_C_C_hishold=zeros(SimTime,1);
Spike_C_C_hisraw=zeros(RRg_length,CC_length,Op_length,SimTime);
Spike_C_C_his=zeros(RRg_length,CC_length,Op_length,(SimTime/n));

for i=1:RRg_length
    for j=1:CC_length
        for k=1:Op_length
            Spike_C_C_hold=squeeze(Spike_C_C_Times(i,j,k,:,:));
            Spike_C_C_hisraw(i,j,k,:)=Spike_C_C_hold*ones(SimNum,1);
            Spike_C_C_hishold=squeeze(Spike_C_C_hisraw(i,j,k,:));
            x=reshape(Spike_C_C_hishold,n,[]);
            Spike_C_C_hishold=sum(x);
            Spike_C_C_his(i,j,k,:)=Spike_C_C_hishold;
        end
    end
end

clearvars -except CC_length Op_length RRg_length SimNum SimTime ... 
Spike_C_A_his Spike_C_B_his Spike_C_C_his n;

%%
%Calcuate Oscillation, Propagation, Optimization
for i=1:RRg_length
    for j=1:CC_length
        for k=1:Op_length
            status=[i j k];
  Spike_C_A_plot=squeeze(Spike_C_A_his(i,j,k,:));
  Spike_C_B_plot=squeeze(Spike_C_B_his(i,j,k,:));
  Spike_C_C_plot=squeeze(Spike_C_C_his(i,j,k,:));
  
%Propagation
 
prop(i,j,k)=max(Spike_C_C_plot(58:73));

[M I_A]=max(Spike_C_A_plot(40:47));
peak_A=39+I_A;
[M I_B]=max(Spike_C_B_plot(47:60));
peak_B=46+I_B;
period_his=peak_B - peak_A;
prop_period(i,j,k)=period_his*n;  %%%Period 

            
%Oscillation
Spike_C_C_dtr=detrend(Spike_C_C_plot(55:100));
Spike_C_C_dtr(Spike_C_C_dtr<0)=0;
Spike_C_C_acorr=xcorr(Spike_C_C_dtr);
osci(i,j,k)=max(Spike_C_C_acorr(56:60))/max(Spike_C_C_acorr);
[M I]=max(Spike_C_C_acorr(54:60));
osci_period(i,j,k)=10*(7+I);

        end
    end
end
  
%Optimization

max_osci=0.7343;
max_prop=1057;
 
for i=1:RRg_length
    for j=1:CC_length
        for k=1:Op_length     

opti(i,j,k)=sqrt((osci(i,j,k)/max_osci)^2 + (prop(i,j,k)/max_prop)^2)-abs(prop(i,j,k)/max_prop - osci(i,j,k)/max_osci);
        end
    end
end
