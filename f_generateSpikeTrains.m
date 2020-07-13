function [t,tu_new,p_new] = f_generateSpikeTrains(t_vec,tro,freq,spikes,w,plot_flag)

%%% INPUT ARGUMENTS %%%

%   t_vec: time vector for kenral (ms)
%   tro: time at response onset/ begining of non-stocastic spiking (ms)
%   freq: frequency of spiking during onset period (Hz)
%   spikes: number of spikes during onset period
%   w: breakdown of falling an rising time constants where w = sqrt(tau1^2
%   + tau2^2) (ms)
%   plot_flag: binary flag for plotting spike train (1 = plot train, 0 = no plot)

%%% OUTPUT ARGUMENTS %%%

%   t: spike times (ms)
%   tu_new: time vector of underlying rate function, shifted by tro
%   p_new: vector of underlying rate function response, shifted by tro

%%%%%%%%%%%%%%%%%%%%%%%%

tu = t_vec; 
tf = tu(end); % final time
b = freq/1e3; 
A = spikes;

tau2 = w/sqrt(5); % falling time constant
tau1 = 2*tau2; % rising time constant

for i = 1:numel(tu)
    if tu(i) < 0
        B(i) = 0;
    else
        B(i) = (1/(tau1 - tau2))*(exp(-tu(i)/tau1) - exp(-tu(i)/tau2));
    end
    p(i) = b + A*B(i);
end

Q = trapz(tu,B); % make sure B has unit area

% shift in time according to tro
tu_shift = tu + tro;
tu_idle = tu(tu < tro);
tu_new = horzcat(tu_idle,tu_shift(tu_shift < tf));

p_idle = b*ones(1,numel(tu_idle));
p_new = horzcat(p_idle,p);
p_new = p_new(1:numel(tu_new));

% Generate the spike trains from the underlying rate function(s) above

tf = 1000; % ms, time to stop simulation
t0 = 0;
n  = 0;
while t0<tf
    n = n + 1;
    if t0 == 0
        lam = b; % essentially the baseline of p_new
    else
        lam = p_new(find(tu_new > t0,1));
    end
    t0 = t0 - (1/lam*log(1-rand));
    if t0 < tf
        t(n) = t0; % output spike times (ms)
    end
end

if (plot_flag == 1)
    stem(t,ones(numel(t),1),'k','Marker','none','LineWidth',1)
    xlim([0 1000])
    set(gca,'YTickLabel',[]) % remove y tick labels
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'YTick',[],'XTick',[])
end
end