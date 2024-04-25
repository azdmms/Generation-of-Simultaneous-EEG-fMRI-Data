clc,clear all ; close all;

%% 1.11 : algorythm for making ERP-Matrix
num_trains = 250 ;
num_trains2 = 100 ;
timeStepS = 0.001;          % 1 msec is short emough to have just one spike in that period

%%% making spike rate distribution 

mean  = 10;    % mean
median = 8;     % median 
mu = log(median); % calculating mu 
sigma = sqrt(2*log(mean/median));% calculating sigma
pd = makedist('Lognormal','mu',mu,'sigma',sigma);%Create a Lognormal Distribution Object Using Specified Parameters
%rng(47) ; % For reproducibility
spikesPerS = random(pd,num_trains*110,1);
load 'leadfield.mat'
gain=gain(:,1:3:15003);
%% BSC PROJECT plot part 
figure();
subplot(2,1,1)
nbins = 50 ;
h1 = histogram(spikesPerS,nbins,'Normalization','probability')
xlabel('spike rate');ylabel('probabilty');title( "spike rate histogram")
subplot(2,1,2)
h2 = histogram(log(spikesPerS),nbins,'Normalization','probability')
xlabel('spike rate');ylabel('probabilty');title( "log spike rate histogram") 
%%
for k=1:1
    %spikesPerS = randi([10, 100], 1,num_trains) ;            % 50 spikes per second, on average it is like r(t)
    durationS = 50.000;          % 1 sec simulation
    num_activclust=5;
    [ERP_matrix0,exitatory0]=gen_sources(pd,num_activclust,durationS,40,40,10);
    [ERP_matrix,exitatory]=gen_sources(pd,num_trains,durationS,10,10,3);
    for i=1:num_trains
        ERP_matrix(i,:)=(ERP_matrix(i,:)*26+ERP_matrix0(1+floor((i-1)/50),:)*100)/126;
        exitatory(i,:)=(exitatory(i,:)*13+exitatory0(1+floor((i-1)/50),:)*50)/63;
    end
    %noisy distributed sources

    times = 0:timeStepS:10;
    t_peak = 5*1e-3 ; 
    d = times.*(exp(-1*times/t_peak)) ;
    d_convolv = [ zeros(1,3) ,d]; 
    d = [times.*(exp(-1*times/t_peak)) ,zeros(1,3)] ;
    d2 = d- d_convolv;
    ERP_matrix2 = zeros( num_trains2 , length(d)+durationS*1000);
    exitatory2 = zeros( num_trains2 ,length(d)+durationS*1000 ); 
    for train = 1: num_trains2
        s11 = zeros( 1, length(d)+durationS*1000);
        s12 = zeros( 1 , length(d)+durationS*1000);
        s13 = zeros( 1 , length(d_convolv)+durationS*1000);
        s14 = zeros( 1 , length(d_convolv)+durationS*1000);
        for i = 1:50
            tslot=1;
            n=floor(durationS/tslot);
            spikesPerS = random(pd,n,1);
            spikes = spike_generator(timeStepS, spikesPerS, durationS,tslot);
            s11 = s11 + conv(spikes, d ) ;

        end
        for i=51:100 
            tslot=1;
            n=floor(durationS/tslot);
            spikesPerS = random(pd,n,1);
            spikes = spike_generator(timeStepS, spikesPerS, durationS,tslot);
            s12 = s12 + conv(spikes, -1*d ) ;
        end
        %s14=s11;
        for i=101:110
            tslot=1;
            n=floor(durationS/tslot);
            spikesPerS = random(pd,n,1);
            spikes = spike_generator(timeStepS, spikesPerS, durationS,tslot);
            s13 = s13 + conv(spikes,d2 ) ;
            s14 = s14 + conv(spikes,d ) ;
        end

        ERP_matrix2(train , : )= (s11 + s12 + s13) /120 ;
        exitatory2(train, : ) = (s11+s14)/100 ; 
    end
    ERP_matrix=[ERP_matrix;ERP_matrix2];
    exitatory=[exitatory;exitatory2];


    sina = 1

    % 2.0.0 pre-processing for fmri data
    TR_forward = 1;            % time resouloution 
    fs = 1/timeStepS;
    indx = floor(size(ERP_matrix,2)/(TR_forward*fs));          % each TR * fs should be one sample
    ERP__resol_mat1=zeros(num_trains+num_trains2 , indx);
    fin111=zeros(num_trains , indx+20);
    for train = 1 : 1 : size(ERP_matrix,1)
        counter = 0 ;
        for i= 1 :  TR_forward *fs : (indx-1)*TR_forward*fs
            counter = counter + 1 ;
            e = ERP_matrix(train , [i : (i + TR_forward *fs -1)]);
            ERP__resol_mat1(train,counter) = sum((e.*e))/(TR_forward *fs);
        end
        fin2 = (ERP__resol_mat1(train, :));
        fin1_indx =find(fin2~=0);                                   % Nonzero Elements
        fin1 = fin2(fin1_indx);
        min_fin1= min(fin1([1:length(fin1_indx)-2]));
        max_fin1= max(fin2);
        fin1=(fin2-min(fin2)) /max_fin1;
        fin111(train,:)=[ zeros(1,20) fin1 ]; 
    end
    %%%
    ERP_p1_resol = zeros(num_trains , indx);
    fin222=zeros(num_trains , indx+20);
    for train = 1 : 1 : size(exitatory,1)
        counter = 0 ;
        for i= 1 :  TR_forward *fs : (indx-1)*TR_forward*fs
            counter = counter + 1 ;
            e = exitatory(train , [i : (i + TR_forward *fs -1)]);
            ERP_p1_resol(train,counter) = sum((e))/(TR_forward *fs);
        end
        fin2 = (ERP_p1_resol(train, :));
        fin1_indx =find(fin2~=0);                                   % Nonzero Elements
        fin1 = fin2(fin1_indx);
        min_fin1= min(fin1([1:length(fin1_indx)-2]));
        max_fin1= max(fin2);
        fin1=(fin2-min(fin2)) /max_fin1;
        fin222(train,:)=[ zeros(1,20) fin1 ]; 
    end
    sina = 2
    %%%%
%     ERP_p3_resol = zeros(num_trains , indx);
%     fin333=zeros(num_trains , indx+20);
%     for train = 1 : 1 : size(ERP_matrix,1)
%         counter = 0 ;
%         for i= 1 :  TR_forward *fs : (indx-1)*TR_forward*fs
%             counter = counter + 1 ;
%             e = ERP_matrix(train , [i : (i + TR_forward *fs -1)]);
%             ERP_p3_resol(train,counter) = sum((e))/(TR_forward *fs);
%         end
%         fin2 = (ERP_p3_resol(train, :));
%         fin1_indx =find(fin2~=0);                                   % Nonzero Elements
%         fin1 = fin2(fin1_indx);
%         min_fin1= min(fin1([1:length(fin1_indx)-2]));
%         max_fin1= max(fin2);
%         fin1=(fin2-min(fin2)) /max_fin1;
%         fin333(train,:)=[ zeros(1,20) fin1 ]; 
%     end
    % %% 2.1.1 :initializations for straight forward atitude for generating fmri data 
    T_alg = 110;
    fmri_data = fmri_generator(T_alg,TR_forward,fin111);
    sina =3
    % 2.1.2 : plotting fmri data in straight forward atitiude 
    figure();
    for train = 1 : 1 :49 %size(fmri_data,1)
        subplot(7,7,train)
        plot(fmri_data(train,:))
        title(" FMRI data ","fontsize",5)
    end
    figure();
    plot(fmri_data(1,:))
    %% 2.2.1 : initializations for better atitude (balloon model) for generating fmri data
    T =110;                                                 % should be the same with T of part 2.1.1
    bold = zeros(num_trains,(T)*10);     %(T-1)*10                      % 
    for train = 1 : 1 : size(ERP__resol_mat1,1)
        E0 = 0.4; 
        RT = 2;
        DUR = 1;
        dCBF = 0.7;
        tao0 = 2; %3.5
        taov = 40; 
        alpha = 0.4;
        fin =fin222(train,:);%*dCBF+1; *1.6
        %%% making the kernel
        tcharge=RT+DUR;
        decharge_time= RT;
        taw_charge= tcharge/3;
        taw_decharge = decharge_time/5;
        fink = (dCBF * [[0:RT]/RT ones(1,DUR) [RT:-1:0]/RT]);
        domain = 1-exp(-1*tcharge/taw_charge);
        fink = [1-exp(-1*[0:RT+DUR]./taw_charge) domain*exp(-[1:RT]./taw_decharge)];%fink = [1-exp(-1*[0:0.1:tcharge]./taw_charge) domain*exp(-[0.1:0.1:decharge_time]./taw_decharge)];
        %%%
        fin =  conv(fin, fink)+1;
    %     tspan = [1,T];
    %     y0=[0;0];
    %     epsi=0.7;  taus=1; tauf=2;
    %     ode = @(t,y) cbf_ode(t,y,epsi,taus,tauf,fin0);
    %     [t,y] = ode113(ode, tspan, y0);
    %     tt = [1:1:T];
    %     fin = 1+interp1(t,y(:,1),tt) ;
        %%%% making m(t)
        m = fin111(train,:);
        m2=conv(m, fink)/2+1; 
        Mu = struct;
        Mu.m = m2;
        Mu.fin = fin;
        Mu.E0 = E0;
        Mu.tao0 = tao0;
        Mu.taov = taov;
        Mu.alpha = alpha;
        tspan = [1,85];           
        y0 = [1; 1; 1];
        ode1 = @(t,y) balloon_ode1(t,y,Mu);
        [t1,y1] = ode45(ode1, tspan, y0);
        tt = [0.1:0.1:T];
        hbr1 = interp1(t1,y1(:,1),tt) ;
        v1 = interp1(t1,y1(:,2),tt);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        b1 = (7*E0*(1-hbr1)+2*(1-hbr1./v1)+(2*E0-.2)*(1-v1));
%         b1 = 2*(3.4*(1-hbr)-1*(1-v));
        %b1 = 7*E0*((1-hbr)-(1-v));

        bold(train,:)=b1;
    end
     save(['bold',num2str(1)],'bold')
     %%
    figure(100)
    for train = 1 : 1 :4 %size(fmri_data,1)
        voxel_number = train;
        subplot(2,2,voxel_number)
        a=bold(voxel_number,:);
        plot (tt,a,'b')
        hold on 
        plot(tt,fmri_data(voxel_number,:),'r')
        %legend({'balloon','hrf-conv'})
        xlabel("t(s)")
        ylabel("BOLD(%)")
    end
%% plotting color map 

%%

    % produce eeg signals:
    xx=randperm(length(loc2),5);
    centers=loc2(xx,:);
    clusts=[];
    locs=[];
    for i=1:5
        [cluster,loc4]=choose_cluster(loc2,centers(i),50);
        clusts=[clusts,cluster];
        locs=[locs;loc4];
    end
    scatts=randperm(length(loc2),num_trains2);
    scatlocs=loc2(scatts,:);
    locs=[locs;scatlocs];
    gain2=gain(:,[clusts,scatts]);
    eeg=gain2*ERP_matrix(:,1:50000);
    bold2=bold(:,200:800);
    save(['eeg',num2str(k)],'eeg')
    save(['locs',num2str(k)],'locs')
    save(['bold',num2str(k)],'bold2')
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikes = spike_generator(timeStepS, spikesPerS, durationS, tslot)
    times = [0:timeStepS:durationS];
    n=floor(durationS/tslot);
    m=length(times);
    spikes = zeros(1,m );
    tslot=tslot*1000;
    prob_spike_occurance_in_timestep = timeStepS * spikesPerS ;
    vt = rand(size(times));
    for slot=1:n
       spikes( (slot-1)*tslot+1:slot*tslot) = vt((slot-1)*tslot+1:slot*tslot) < (prob_spike_occurance_in_timestep(slot)) ;
    end   
     spikes( n*tslot+1:m) = vt(n*tslot+1:m) < (prob_spike_occurance_in_timestep(slot)) ;
end
%%%
function rasterPlot(spikes, timeStepS)

figure(1);

times = [0:timeStepS:timeStepS * (length(spikes) - 1)];
axes('position', [0.1, 0.1, 0.8, 0.8]);
axis([0, length(spikes) - 1, 0, 1]);
trains = size(spikes, 1); 
ticMargin = 0.01;                                      % gap between spike trains (full scale is 1)
ticHeight = (1.0 - (trains + 1) * ticMargin) / trains;

for train = 1:trains
    spikeTimes = find(spikes(train, :) == 1);
    yOffset = ticMargin + (train - 1) * (ticMargin + ticHeight);
    for i = 1:length(spikeTimes)
        line([spikeTimes(i), spikeTimes(i)], [yOffset, yOffset + ticHeight], 'linewidth',2);
    end
end

xlabel('Time ')
title('Raster plot of spikes with 50 active sources ','fontsize',12);
end
%%%
 function fmri_data_matrix = fmri_generator(T_alg,TR,spike_train)
 % T_alg = duration of simulation
 % TR : Time resouloution
 hrf = load('hrf2.mat');
 tt = [0.1:0.1:T_alg];%[1.1:0.1:T_alg];
 hrf=hrf.hrf2;
 hrf = resample(hrf , 1,1000);
 %hrf = conv(hrf,ones(1,max(1,round(duration/trold))));
 % resample to desired TR
 %hrf = interp1((0:length(hrf)-1)*trold,hrf,0:TR:(length(hrf)-1)*trold,'pchip');
 % make the peak equal to one
 %hrf = hrf / max(hrf);
 hrf = transpose(hrf);
t=0:99;
p=8.6; q=0.547;
y=(t/(p*q)).^p .* exp(p-t/q);
p2=15; q2=1;
r=1/6;
y2=(t/(p2*q2)).^p2 .* exp(p2-t/q2);
hrf2=y-r*y2;
 %%% fmri data 
 fmri_data_matrix = zeros( size(spike_train,1) , length(tt));
 for train = 1: 1 : size(spike_train,1) 
     zz= conv(spike_train(train, :), hrf2 ) ;
     t12=[1:1:length(zz)];
     fmri_data_matrix(train , : )=interp1(t12,zz,tt);
 end
 end
 %%%
 function dydt = balloon_ode1(t,y,Mu)
% Based on:
% Mildner et al (2001) A qualitative test of the balloon model for BOLD-based MR signal changes at 3T. Magnetic resonance in medicine 46 (5) 891-9
% y(1): q (HbR)
% y(2): v (CBV)
% y(3): p (HbO) + q (HbR)
%
% This is the function to be called by ode45

t1 = floor(t);
t2 = ceil(t);
if(t1==t2)
    fin_t = Mu.fin(t);
    m_t = Mu.m(t);
else
    fin_t = interp1([t1 t2], [Mu.fin(t1) Mu.fin(t2)], t);
    m_t = interp1([t1 t2], [Mu.m(t1) Mu.m(t2)], t);
end

E0 = Mu.E0;
tao0 = Mu.tao0;
taov = Mu.taov;
alpha = Mu.alpha;

q = y(1);
v = y(2);
h = y(3);
dydt = [0;0;0]; 
f_out=v.^(1/alpha);
E_t = 1-(1-E0).^(1/m_t);
dydt(1) = 1/tao0 * (m_t-fin_t*q/v);% + 1/taov * (fin_t - v.^(1/alpha))*q/v;
dydt(2) =  1/taov * (fin_t - f_out);
%f_out=v.^(1/alpha)+taov*dydt(2);
%dydt(1) =  m_t/tao0 * (E_t/E0 - q/v)+ 1/taov * (fin_t - v.^(1/alpha))*q/v; %1/tao0 * (m_t-f_out*q/v);
dydt(3) = 1/taov * (fin_t - v.^(1/alpha)*h/v);   

%dydt(3) = fin_t/tao0 * ((1-E_t)/E0 - h/v) + 1/taov * (fin_t - v.^(1/alpha))*h/v;
 end
 
function [ERP_matrix,exitatory]=gen_sources(pd,num_trains,durationS,num_exci,num_inhi,num_mix)
timeStepS = 0.001;
times = 0:timeStepS:10;	% a vector with each time step

% % figure()
%trains = 22000; %size(spikes, 1);
t_peak = 5*1e-3 ; 
d = times.*(exp(-1*times/t_peak)) ;
d_convolv = [ zeros(1,3) ,d]; 
d = [times.*(exp(-1*times/t_peak)) ,zeros(1,3)] ;
d2 = d- d_convolv;
ERP_matrix = zeros( num_trains , length(d)+durationS*1000);
exitatory = zeros( num_trains ,length(d)+durationS*1000 );  %length(conv(spikes(1, :), d ))
for train = 1: num_trains
    %spikesPerS = random(pd,1,1);
    tslot=1;
    n=floor(durationS/tslot);
    spikesPerS = random(pd,n,1);
    s11 = zeros( 1, length(d)+durationS*1000);
    s12 = zeros( 1 , length(d)+durationS*1000);
    s13 = zeros( 1 , length(d_convolv)+durationS*1000);
    s14 = zeros( 1 , length(d_convolv)+durationS*1000);
    for i = 1:num_exci
        spikes = spike_generator(timeStepS, spikesPerS, durationS,tslot);
        s11 = s11 + conv(spikes, d ) ;
    end
    for i=1:num_inhi 
        spikes = spike_generator(timeStepS, spikesPerS, durationS,tslot);
        s12 = s12 + conv(spikes, -1*d ) ;
    end
    %s14=s11;
    for i=1:num_mix
        spikes = spike_generator(timeStepS, spikesPerS, durationS,tslot);
        s13 = s13 + conv(spikes,d2 ) ;
        s14 = s14 + conv(spikes,d ) ;
    end

    ERP_matrix(train , : )= (s11 + s12 + s13) /(num_exci+num_inhi+num_mix*2) ;
    exitatory(train, : ) = (s11+s14)/(num_exci+num_mix) ; 
end

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function dydt = balloon_ode2(t,y,Mu)
% Based on:
% Mildner et al (2001) A qualitative test of the balloon model for BOLD-based MR signal changes at 3T. Magnetic resonance in medicine 46 (5) 891-9
% y(1): q (HbR)
% y(2): v (CBV)
% y(3): p (HbO) + q (HbR)
%
% This is the function to be called by ode45

t1 = floor(t);
t2 = ceil(t);
if(t1==t2)
    m_t = Mu.m(t);
else
    m_t = interp1([t1 t2], [Mu.m(t1) Mu.m(t2)], t);
end

E0 = Mu.E0;
tao0 = Mu.tao0;
taov = Mu.taov;
alpha = Mu.alpha;

q = y(1);
v = y(2);
h = y(3);
dydt = [0;0;0]; 
E_t = 1-(1-E0).^(1/fin_t);
dydt(1) = m_t/tao0 * (E_t/E0 - q/v) + 1/taov * (m_t - v.^(1/alpha))*q/v;
dydt(2) =  1/taov * (m_t - v.^(1/alpha));
dydt(3) = 1/taov * (fin_t - v.^(1/alpha)*h/v);           
%dydt(3) = fin_t/tao0 * ((1-E_t)/E0 - h/v) + 1/taov * (fin_t - v.^(1/alpha))*h/v;
 end

