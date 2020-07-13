% Taylor Hansen
% Mitch Thomas
% BIOEN 6005 Project

clearvars;
close all;

%% Step 1: Set parameters for underlying rate function and generate a spike train

tu        = linspace(0,1000,10000); % time vector for underlying rate function (ms)
tro       = 400; % time at response onset (ms)
freq      = 20; % desired frequency of spiking during onset period
spikes    = 20; % deisred number of spikes during onset period
w         = 100; % response width (determines falling and rising time constants) (ms)
plot_flag = 0; % binary flag for plotting spike train (1 = plot)

[ts,tu,p] = f_generateSpikeTrains(tu,tro,freq,spikes,w,plot_flag); % note: times are returned in ms
ts = ts/1000; % convert spike times to s

% plot the underlying rate function
figure;
subplot(4,1,1)
plot(tu,p*1e3,'k')
% set(gca,'YTickLabel',[]) % remove y tick labels
set(gca,'XTickLabel',[]) % remove x tick labels
% set(gca,'YTick',[],'XTick',[])
set(gca,'Xtick',[])
xlim([0 1000])
% xlabel('Time, ms')
ylabel('\rho, Hz')
title('Underlying Rate Function')

%% Step 2: From spike train generated above, convolve with kernel to estimate underlying firing rate

sig = 50; % Kernel width, ms
tk = linspace(-1000,1000,10000); % t for kernel, ms

[K,sup] = generateKernels(tk,sig);

% plot the generated spike train
subplot(4,1,2)
stem(ts,ones(numel(ts),1),'k','Marker','none','LineWidth',1)
xlim([0 1])
set(gca,'YTickLabel',[]) % remove y tick labels
set(gca,'XTickLabel',[]) % remove x tick labels
set(gca,'YTick',[],'XTick',[])
title('Spike Train')

% figure;
subplot(4,1,3)
% set variables to work with code from Initial_development.m
i = 1;
T{i} = ts;
t = tk;

% overlay a simple triangular kernel on top of spike train to ensure proper timing
if numel(T{i}) > 1
    sup_t_idx = t >= sup.triangle(1) & t <= sup.triangle(2);
    sup_t = t(sup_t_idx);

    frate{i} = 0;
    for j = 1:numel(T{i})
        % get vector for current kernel
        currK = K.triangle(sup_t_idx);

        % find closest value in t to current spike time
        spikeloc = find(t > T{i}(j)*1000,1);
        % place it appropriately in time
        idx1 = find(sup_t_idx,1);
        idx3  = find(sup_t_idx,1,'last');
        idx2   = floor((idx1 + idx3)/2);

        fullK = zeros(1,numel(t));
        fullK(idx1 + (spikeloc-idx2) : idx3 + (spikeloc-idx2)) = currK;

        % sum kernels to get estimate of firing rate (in KHz)
        frate{i} = frate{i} + fullK(1:numel(t)); % prevent overflow

        plot((sup_t/1000)+T{i}(j),K.triangle(sup_t_idx),'k');
        hold on
    end
    hold off
    xlim([0 1])
    set(gca,'YTickLabel',[]) % remove y tick labels
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'YTick',[],'XTick',[])
    title('Triangular kernels')

    subplot(4,1,4)
    % plot estimate of firing rate
    plot(t/1000,frate{i}*1000,'k')
    xlim([0 1])
    xlabel('Time, s')
    ylabel('\lambda, Hz')
    title('Estimated Firing Rate')
    tmpylim = get(gca,'YLim');
    if max(frate{i}*1000) > 0.95*tmpylim(end)
        ylim([0 tmpylim(end)*1.1])
    end
end

%% Step 3: Calculate Integrated Square Error (ISE) between estimated and true firing rate

ISE = sum((frate{i} - p).^2);
ISE = ISE/numel(T{i}); % may want to divide by total number of spikes

% may want to use (max(frate{i}) - max(p)).^2 as error statistic
Peak_err = (max(frate{i}) - max(p)).^2;