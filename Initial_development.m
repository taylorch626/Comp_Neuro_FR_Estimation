% Taylor Hansen
% Mitch Thomas
% BIOEN 6005 Project

clearvars;
close all;

%% Generate figure for varying average firing rates (triangular kernel)

clear K sup T frate

sig = 50; % Kernel width, ms
t = linspace(-1000,1000,10000); % t for kernel, ms

[K,sup] = generateKernels(t,sig);

% Examples of simple single spike trains (NOT generated in NEURON)

% 4/7/19 - Just use Poisson-generated ones from PS3.4
lambda = [3 10 30 100]; % Poisson firing rates, Hz
figure;
for i = 1:numel(lambda)
    T{i} = spiketime(lambda(i));
    T{i} = cumsum(T{i});
    if i <= 2
        subplot(8,2,i)
    else
        subplot(8,2,i+6)
    end
    stem(T{i}(2:end),ones(1,numel(T{i}(2:end))),'k','Marker','none','LineWidth',1)
    xlim([0 1])
    set(gca,'YTickLabel',[]) % remove y tick labels
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'YTick',[],'XTick',[])
    title(sprintf('Spike train for average \\lambda = %d Hz',lambda(i)))
    
    % overlay a simple kernel on top of spike train to ensure proper timing
    if numel(T{i}) > 1
        if i <= 2
            subplot(8,2,i+2)
        else
            subplot(8,2,i+8)
        end
        sup_t_idx = t >= sup.triangle(1) & t <= sup.triangle(2);
        sup_t = t(sup_t_idx);
        
        frate{i} = 0;
        for j = 2:numel(T{i})
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
            
            % test if integration is truly 1
            Q = trapz(sup_t+T{i}(j)*1000,K.triangle(sup_t_idx));
        end
        xlim([0 1])
%         xlabel('Time, s')
        set(gca,'YTickLabel',[]) % remove y tick labels
        set(gca,'XTickLabel',[]) % remove x tick labels
        set(gca,'YTick',[],'XTick',[])
        title('Triangular kernels')
        
        % plot estimate of firing rate
        if i <= 2
            subplot(8,2,i+4)
        else
            subplot(8,2,i+10)
        end
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
    hold off
end

%% Generate figure for varying sigma (for triangular kernel)

clear K sup T frate

sig = [20, 50, 100]; % Kernel width, ms
% sig = [10, 20, 300];
t = linspace(-1000,1000,10000); % t for kernel, ms

% Examples of simple single spike trains (NOT generated in NEURON)

% 4/7/19 - Just use Poisson-generated ones from PS3.4
lambda = 30; % Poisson firing rates, Hz
T = spiketime(lambda);
T = cumsum(T);

figure;
for h = 1:numel(sig)
    [K,sup] = generateKernels(t,sig(h));

    subplot(3,3,h)
    stem(T(2:end),ones(1,numel(T(2:end))),'k','Marker','none','LineWidth',1)
    xlim([0 1])
    set(gca,'YTickLabel',[]) % remove y tick labels
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'YTick',[],'XTick',[])
    title(sprintf('\\sigma = %d ms',sig(h)))

    % overlay a simple kernel on top of spike train to ensure proper timing
    if numel(T) > 1
        subplot(3,3,h+3)
        sup_t_idx = t >= sup.triangle(1) & t <= sup.triangle(2);
        sup_t = t(sup_t_idx);

        frate{h} = 0;
        for j = 2:numel(T)
            % get vector for current kernel
            currK = K.triangle(sup_t_idx);

            % find closest value in t to current spike time
            spikeloc = find(t > T(j)*1000,1);
            % place it appropriately in time
            idx1 = find(sup_t_idx,1);
            idx3  = find(sup_t_idx,1,'last');
            idx2   = floor((idx1 + idx3)/2);

            fullK = zeros(1,numel(t));
            fullK(idx1 + (spikeloc-idx2) : idx3 + (spikeloc-idx2)) = currK;

            % sum kernels to get estimate of firing rate (in KHz)
            frate{h} = frate{h} + fullK(1:numel(t)); % prevent overflow

            plot((sup_t/1000)+T(j),K.triangle(sup_t_idx),'k');
            hold on

            % test if integration is truly 1
            Q = trapz(sup_t+T(j)*1000,K.triangle(sup_t_idx));
        end
        xlim([0 1])
        set(gca,'YTickLabel',[]) % remove y tick labels
        set(gca,'XTickLabel',[]) % remove x tick labels
        set(gca,'YTick',[],'XTick',[])
        title('Triangular kernels')
        
        if h == 1
            tmp = get(gca,'YLim');
        else
            ylim(tmp)
        end

        % plot estimate of firing rate
        subplot(3,3,h+6)
        plot(t/1000,frate{h}*1000,'k')
        xlim([0 1])
        xlabel('Time, s')
        ylabel('\lambda, Hz')
        title('Estimated Firing Rate')
        tmpylim = get(gca,'YLim');
        if max(frate{h}*1000) > 0.95*tmpylim(end)
            ylim([0 tmpylim(end)*1.1])
        end
    end
    hold off
end

%% Generate figure for varying kernel shape

clear K sup T frate sup_t_idx sup_t
clearvars

sig = 50; % Kernel width, ms
t = linspace(-1000,1000,10000); % t for kernel, ms

[K,sup] = generateKernels(t,sig);

% Examples of simple single spike trains (NOT generated in NEURON)

% 4/7/19 - Just use Poisson-generated ones from PS3.4
lambda = 30; % Poisson firing rates, Hz
T = spiketime(lambda);
T = cumsum(T);

figure;
% plot example spike train
subplot(6,1,1)
stem(T(2:end),ones(1,numel(T(2:end))),'k','Marker','none','LineWidth',1)
xlim([0 1])
set(gca,'YTickLabel',[]) % remove y tick labels
set(gca,'XTickLabel',[]) % remove x tick labels
set(gca,'YTick',[],'XTick',[])
title(sprintf('Average %d Hz Spike Train',lambda))

% plot kernels for various shapes
if numel(T) > 1

    sup_t_idx.boxcar   = t >= sup.boxcar(1) & t <= sup.boxcar(2);
    sup_t_idx.triangle = t >= sup.triangle(1) & t <= sup.triangle(2);
    sup_t_idx.epan     = t >= sup.epan(1) & t <= sup.epan(2);
    sup_t_idx.gauss    = t >= sup.gauss(1) & t <= sup.gauss(2);
    
    sup_t.boxcar   = t(sup_t_idx.boxcar);
    sup_t.triangle = t(sup_t_idx.triangle);
    sup_t.epan     = t(sup_t_idx.epan);
    sup_t.gauss    = t(sup_t_idx.gauss);

    frate.boxcar   = 0;
    frate.triangle = 0;
    frate.epan     = 0;
    frate.gauss    = 0;
    for j = 2:numel(T)
        % get vector for current kernel
        currK.boxcar   = K.boxcar*ones(1,numel(sup_t.boxcar));
        currK.triangle = K.triangle(sup_t_idx.triangle);
        currK.epan     = K.epan(sup_t_idx.epan);
        currK.gauss    = K.gauss(sup_t_idx.gauss);

        % find closest value in t to current spike time
        spikeloc = find(t > T(j)*1000,1);
        % place it appropriately in time
        idx1.boxcar   = find(sup_t_idx.boxcar,1);
        idx1.triangle = find(sup_t_idx.triangle,1);
        idx1.epan     = find(sup_t_idx.epan,1);
        idx1.gauss    = find(sup_t_idx.gauss,1);
        
        idx3.boxcar   = find(sup_t_idx.boxcar,1,'last');
        idx3.triangle = find(sup_t_idx.triangle,1,'last');
        idx3.epan     = find(sup_t_idx.epan,1,'last');
        idx3.gauss    = find(sup_t_idx.gauss,1,'last');
        
        idx2.boxcar   = floor((idx1.boxcar + idx3.boxcar)/2);
        idx2.triangle = floor((idx1.triangle + idx3.triangle)/2);
        idx2.epan     = floor((idx1.epan + idx3.epan)/2);
        idx2.gauss    = floor((idx1.gauss + idx3.gauss)/2);

        fullK.boxcar = zeros(1,numel(t));
        fullK.boxcar(idx1.boxcar + (spikeloc-idx2.boxcar) : idx3.boxcar + (spikeloc-idx2.boxcar)) = currK.boxcar;        
        fullK.triangle = zeros(1,numel(t));
        fullK.triangle(idx1.triangle + (spikeloc-idx2.triangle) : idx3.triangle + (spikeloc-idx2.triangle)) = currK.triangle;
        fullK.epan = zeros(1,numel(t));
        fullK.epan(idx1.epan + (spikeloc-idx2.epan) : idx3.epan + (spikeloc-idx2.epan)) = currK.epan;
        fullK.gauss = zeros(1,numel(t));
        fullK.gauss(idx1.gauss + (spikeloc-idx2.gauss) : idx3.gauss + (spikeloc-idx2.gauss)) = currK.gauss;        

        % sum kernels to get estimate of firing rate (in KHz)
        frate.boxcar   = frate.boxcar + fullK.boxcar(1:numel(t)); % prevent overflow
        frate.triangle = frate.triangle + fullK.triangle(1:numel(t)); % prevent overflow
        frate.epan     = frate.epan + fullK.epan(1:numel(t)); % prevent overflow
        frate.gauss    = frate.gauss + fullK.gauss(1:numel(t)); % prevent overflow
    end
    
    % get consistent y lim for all kernels
    ymax = max([max(K.boxcar*ones(1,numel(sup_t.boxcar))),...
                max(K.triangle(sup_t_idx.triangle)),...
                max(K.epan(sup_t_idx.epan)),...
                max(K.gauss(sup_t_idx.gauss))]);
        
    % plot kernel overlays
    subplot(6,1,2) % Boxcar
    for j = 2:numel(T)
        dt = abs(sup_t.boxcar(2) - sup_t.boxcar(1))/2;
        plot([((sup_t.boxcar(1)-dt)/1000)+T(j), (sup_t.boxcar/1000)+T(j), ((sup_t.boxcar(end)+dt)/1000)+T(j)],...
             [0, K.boxcar*ones(1,numel(sup_t.boxcar)), 0],'k')
        hold on
    end
    hold off
    xlim([0 1])
    ylim([0 ymax*1.1])
    set(gca,'YTickLabel',[]) % remove y tick labels
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'YTick',[],'XTick',[])
    title('Boxcar kernels')
    
    subplot(6,1,3) % Triangle
    for j = 2:numel(T)
        plot((sup_t.triangle/1000)+T(j),K.triangle(sup_t_idx.triangle),'k--');
        hold on
    end
    hold off
    xlim([0 1])
    ylim([0 ymax*1.1])
    set(gca,'YTickLabel',[]) % remove y tick labels
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'YTick',[],'XTick',[])
    title('Triangular kernels')
    
    subplot(6,1,4) % Epanechnikov
    for j = 2:numel(T)
        plot((sup_t.epan/1000)+T(j),K.epan(sup_t_idx.epan),'k:');
        hold on
    end
    hold off
    xlim([0 1])
    ylim([0 ymax*1.1])
    set(gca,'YTickLabel',[]) % remove y tick labels
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'YTick',[],'XTick',[])
    title('Epanechnikov kernels')
    
    subplot(6,1,5) % Gaussian
    for j = 2:numel(T)
        plot((sup_t.gauss/1000)+T(j),K.gauss(sup_t_idx.gauss),'k-.');
        hold on
    end
    hold off
    xlim([0 1])
    ylim([0 ymax*1.1])
    set(gca,'YTickLabel',[]) % remove y tick labels
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'YTick',[],'XTick',[])
    title('Gaussian kernels')

    % plot estimate of firing rate
    subplot(6,1,6)
    plot(t/1000,frate.boxcar*1000,'k')
    hold on
    plot(t/1000,frate.triangle*1000,'k--')
    hold on
    plot(t/1000,frate.epan*1000,'k:')
    hold on
    plot(t/1000,frate.gauss*1000,'k-.')
    hold off
    xlim([0 1])
    xlabel('Time, s')
    ylabel('\lambda, Hz')
    title('Estimated Firing Rate')
    tmpylim = get(gca,'YLim');
    ratemax = max([max(frate.boxcar*1000),...
                   max(frate.triangle*1000),...
                   max(frate.epan*1000),...
                   max(frate.gauss*1000)]);
    if ratemax > 0.95*tmpylim(end)
        ylim([0 tmpylim(end)*1.1])
    end
end

%% Functions called above

function [K,sup] = generateKernels(t,sig)
    % Define some standard kernel types K(t,sig)
    K.boxcar = 1/(2*sig*sqrt(3));
    K.triangle = (1/(6*sig^2))*((sig*sqrt(6)) - abs(t));
    K.epan = (3/(4*sig*sqrt(5)))*(1 - (t.^2/(5*sig^2)));
    K.gauss = (1/(sig*sqrt(2*pi)))*exp(-t.^2/(2*sig^2));

    % Define support vectors for kernel types
    % i.e. domain over which kernel functions assume non-zero values
    sup.boxcar = [-sig*sqrt(3), sig*sqrt(3)];
    sup.triangle = [-sig*sqrt(6), sig*sqrt(6)];
    sup.epan = [-sig*sqrt(5), sig*sqrt(5)];
    sup.gauss = [-inf, inf];
end

function T = spiketime(lambda) % pulled from TH PS3.4
    
    % Invert CDF to find T
    T = 0;
    flag = 0;
    while flag == 0
        % Pick a uniform random number
        alpha = rand(1);
        beta = 1-alpha;
        Tnew = (-log(beta))/lambda;
        if sum([T Tnew]) < 1
            T = [T Tnew];
        else
            flag = 1;
        end
    end
end