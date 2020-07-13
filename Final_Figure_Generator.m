% Taylor Hansen
% Mitch Thomas
% BIOEN 6005 Project

%%% Generate Final Figures for Paper

clearvars;
close all;

%% Figure 1: Varying Underlying Firing Rate

clearvars;
% close all;

tu_init   = linspace(0,1000,10000); % time vector for underlying rate function (ms)
tro       = 400; % time at response onset (ms)
freq      = [20, 30, 40]; % desired frequency of spiking during onset period
spikes    = 40; % desired number of spikes during onset period
w         = [80, 110, 140]; % response width (determines falling and rising time constants) (ms)
plot_flag = 0; % binary flag for plotting spike train (1 = plot)

%%% Generate and plot the underlying firing rates
figure;
set(gcf,'Name','Varying Underlying Firing Rate');
maxY = 0;
for i = 1:numel(freq)
    try
        [ts{i},tu{i},p{i}] = f_generateSpikeTrains(tu_init,tro,freq(i),spikes,w(i),plot_flag); % note: times are returned in ms
    catch ME
        warning('Spike train generation insufficient...rerunning code')
        close all;
        Final_Figure_Generator; % rerun code if it crashes during spike train generation
    end
    ts{i} = ts{i}/1000; % convert spike times to s
    ax{i} = subplot(4,3,i);
    plot(tu{i},p{i}*1e3,'k');
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'Xtick',[])
    xlim([0 1000])
    ylabel('\rho, Hz')
    if max(get(ax{i},'YLim')) > maxY
        maxY = max(get(ax{i},'YLim'));
    end
    for j = 1:i
        set(ax{j},'YLim',[0 maxY*1.1]) % set all ylims to be equivalent
    end
end

%%% From spike trains generated above, convolve with kernel to estimate underlying firing rate

sig = 50; % Kernel width, ms
tk = linspace(-1000,1000,10000); % t for kernel, ms

[K,sup] = generateKernels(tk,sig);

% plot the generated spike train
for i = 1:numel(freq)
    ax{3+i} = subplot(4,3,3+i);
    stem(ts{i},ones(numel(ts{i}),1),'k','Marker','none','LineWidth',0.75,'ShowBaseLine','off')
    xlim([0 1])
end

% figure;
for h = 1:numel(freq)
    ax{6+h} = subplot(4,3,6+h);
    % set variables to work with code from Initial_development.m
    i = 1;
    T{h}{i} = ts{h};
    t = tk;

    % overlay a simple triangular kernel on top of spike train to ensure proper timing
    if numel(T{h}{i}) > 1
        sup_t_idx = t >= sup.triangle(1) & t <= sup.triangle(2);
        sup_t = t(sup_t_idx);

        frate{i} = 0;
        for j = 1:numel(T{h}{i})
            % get vector for current kernel
            currK = K.triangle(sup_t_idx);

            % find closest value in t to current spike time
            spikeloc = find(t > T{h}{i}(j)*1000,1);
            % place it appropriately in time
            idx1 = find(sup_t_idx,1);
            idx3  = find(sup_t_idx,1,'last');
            idx2   = floor((idx1 + idx3)/2);

            fullK = zeros(1,numel(t));
            fullK(idx1 + (spikeloc-idx2) : idx3 + (spikeloc-idx2)) = currK;

            % sum kernels to get estimate of firing rate (in KHz)
            frate{i} = frate{i} + fullK(1:numel(t)); % prevent overflow

            plot((sup_t/1000)+T{h}{i}(j),K.triangle(sup_t_idx),'k');
            hold on
        end
        hold off
        xlim([0 1])
        set(gca,'YTickLabel',[]) % remove y tick labels
        set(gca,'XTickLabel',[]) % remove x tick labels
        set(gca,'YTick',[],'XTick',[])
        title('Triangular kernels')

        ax{9+h} = subplot(4,3,9+h);
        % plot estimate of firing rate
        plot(t/1000,frate{i}*1000,'k')
        xlim([0 1])
        set(ax{9+h},'XTick',[0 0.5 1],'XTickLabel',{'0','0.5','1'})
        xlabel('Time, s')
        ylabel('\lambda, Hz')
        if max(get(ax{9+h},'YLim')) > maxY
            maxY = max(get(ax{9+h},'YLim'));
        end
        for xx = 1:h
            set(ax{xx},'YLim',[0 maxY*1.1]) % set all ylims to be equivalent
            set(ax{9+xx},'YLim',[0 maxY*1.1]) % set all ylims to be equivalent
        end
    end
end

% Adjust positions and layout of figure
for i = 1:12
    set(ax{i},'Fontname','Times')
    set(ax{i},'Fontsize',10)
    if i == 1 || i == 2 || i == 3 || i == 10 || i == 11 || i == 12
        set(ax{i},'Box','off')
        tmp = get(ax{i},'Ylabel');
        if i ~= 1 && i ~= 10
            tmp.String = '';
            set(ax{i},'Ylabel',tmp)
            set(ax{i},'Ytick',[])
        end
        if i == 10
            tmp10 = get(ax{i},'Position');
            tmp10(1,2) = tmp10(1,2)*4.15;
            set(ax{i},'Position',tmp10);
        elseif i > 10
            tmp = get(ax{i},'Position');
            tmp(1,2) = tmp10(1,2);
            set(ax{i},'Position',tmp);
        end
    elseif i == 4 || i == 5 || i == 6
        set(ax{i},'Visible','off')
        if i == 4
            tmp4 = get(ax{i},'Position');
            tmp4(1,4) = tmp4(1,4)/3;
            tmp4(1,2) = tmp4(1,2)*1.27;
            set(ax{i},'Position',tmp4);
        else
            tmp = get(ax{i},'Position');
            tmp(1,4) = tmp4(1,4);
            tmp(1,2) = tmp4(1,2);
            set(ax{i},'Position',tmp);
        end
    elseif i == 7 || i == 8 || i == 9
        set(ax{i},'Visible','off')
        if i == 7
            tmp7 = get(ax{i},'Position');
            tmp7(1,4) = tmp7(1,4)/3;
            tmp7(1,2) = tmp7(1,2)*1.93;
            set(ax{i},'Position',tmp7);
        else
            tmp = get(ax{i},'Position');
            tmp(1,4) = tmp7(1,4);
            tmp(1,2) = tmp7(1,2);
            set(ax{i},'Position',tmp);
        end
    end
end


%% Figure 2: Varying Kernel Width

clearvars;
% close all;

tu        = linspace(0,1000,10000); % time vector for underlying rate function (ms)
tro       = 400; % time at response onset (ms)
freq      = 20; % desired frequency of spiking during onset period
spikes    = 20; % desired number of spikes during onset period
w         = 100; % response width (determines falling and rising time constants) (ms)
plot_flag = 0; % binary flag for plotting spike train (1 = plot)

try
    [ts,tu,p] = f_generateSpikeTrains(tu,tro,freq,spikes,w,plot_flag); % note: times are returned in ms
catch ME
    warning('Spike train generation insufficient...rerunning code')
    close all;
    Final_Figure_Generator; % rerun code if it crashes during spike train generation
end
ts = ts/1000; % convert spike times to s

sig = [20, 50, 100]; % Kernel width, ms

% plot the underlying rate function
figure;
set(gcf,'Name','Varying Kernel Width');
for i = 1:numel(sig)
    if i == 2 % avoid generating redundant plots
        ax{i} = subplot(4,3,i);
        plot(tu,p*1e3,'k')
        set(gca,'Xtick',[])
        xlim([0 1000])
        set(gca,'Xtick',[0 500 1000],'Xticklabel',{'0','0.5','1'})
        xlabel('Time, s')
        ylabel('\rho, Hz')
        maxY = max(get(ax{i},'YLim'));
        ylim([0 1.1*maxY])
    end
end

%%% From spike train generated above, convolve with kernel to estimate underlying firing rate

tk = linspace(-1000,1000,10000); % t for kernel, ms

for h = 1:numel(sig)

    [K,sup] = generateKernels(tk,sig(h));

    % plot the generated spike train
    if h == 2 % avoid generating redundant plots
        ax{3+h} = subplot(4,3,3+h);
        stem(ts,ones(numel(ts),1),'k','Marker','none','LineWidth',0.75,'ShowBaseLine','off')
        xlim([0 1])
        set(gca,'YTickLabel',[]) % remove y tick labels
        set(gca,'XTickLabel',[]) % remove x tick labels
        set(gca,'YTick',[],'XTick',[])
        title('Example Spike Train')
    end

    % figure;
    ax{6+h} = subplot(4,3,6+h);
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
        title(sprintf('Triangular kernels, \\sigma = %d ms',sig(h)))
        if h == 1
            tmp = get(gca,'YLim');
        else
            ylim(tmp)
        end

        ax{9+h} = subplot(4,3,9+h);
        % plot estimate of firing rate
        plot(t/1000,frate{i}*1000,'k')
        xlim([0 1])
        set(gca,'Xtick',[0 0.5 1],'Xticklabel',{'0','0.5','1'})
        xlabel('Time, s')
        ylabel('\lambda, Hz')
        if max(get(ax{9+h},'YLim')) > maxY
            maxY = max(get(ax{9+h},'YLim'));
        end
        for j = 1:h
            set(ax{j},'YLim',[0 maxY*1.1]) % set all ylims to be equivalent
            set(ax{9+j},'YLim',[0 maxY*1.1]) % set all ylims to be equivalent
        end
    end
end

% Adjust positions and layout of figure
for i = 1:12
    set(ax{i},'Fontname','Times')
    set(ax{i},'Fontsize',10)
    if i == 1 || i == 2 || i == 3 || i == 10 || i == 11 || i == 12
        set(ax{i},'Box','off')
        tmp = get(ax{i},'Ylabel');
        if i ~= 2 && i ~= 10
            tmp.String = '';
            set(ax{i},'Ylabel',tmp)
            set(ax{i},'Ytick',[])
        end
        if i == 2
            tmp2 = get(ax{i},'Position');
            tmp2(1,4) = tmp2(1,4)*1.1;
            tmp2(1,2) = tmp2(1,2) + 0.04;
            set(ax{i},'Position',tmp2);
        end
        if i == 10
            tmp10 = get(ax{i},'Position');
            tmp10(1,2) = tmp10(1,2)*2.75;
            tmp10(1,4) = tmp2(1,4);
            set(ax{i},'Position',tmp10);
        elseif i > 10
            tmp = get(ax{i},'Position');
            tmp(1,2) = tmp10(1,2);
            tmp(1,4) = tmp2(1,4);
            set(ax{i},'Position',tmp);
        end
    elseif i == 5
        set(ax{i},'Visible','off')
        tmp = get(ax{i},'Position');
        tmp(1,4) = tmp(1,4)/3;
        tmp(1,2) = tmp(1,2)*1.18;
        set(ax{i},'Position',tmp);
    elseif i == 7 || i == 8 || i == 9
        set(ax{i},'Visible','off')
        if i == 7
            clear tmp7
            tmp7 = get(ax{i},'Position');
            tmp7(1,4) = tmp7(1,4)*4/7;
            tmp7(1,2) = tmp7(1,2)*1.5;
            set(ax{i},'Position',tmp7);
        else
            tmp = get(ax{i},'Position');
            tmp(1,4) = tmp7(1,4);
            tmp(1,2) = tmp7(1,2);
            set(ax{i},'Position',tmp);
        end
    end
end


%% Figure 3: Varying Kernel Shape

clearvars;
% close all;

%%% Generate and plot the underlying firing rates

tu        = linspace(0,1000,10000); % time vector for underlying rate function (ms)
tro       = 400; % time at response onset (ms)
freq      = 20; % desired frequency of spiking during onset period
spikes    = 20; % desired number of spikes during onset period
w         = 100; % response width (determines falling and rising time constants) (ms)
plot_flag = 0; % binary flag for plotting spike train (1 = plot)

try
    [ts,tu,p] = f_generateSpikeTrains(tu,tro,freq,spikes,w,plot_flag); % note: times are returned in ms
catch ME
    warning('Spike train generation insufficient...rerunning code')
    close all;
    Final_Figure_Generator; % rerun code if it crashes during spike train generation
end
ts = ts/1000; % convert spike times to s

figure;
set(gcf,'Name','Varying Kernel Shape');
ax{1} = subplot(7,1,1);
plot(tu,p*1e3,'k')
set(gca,'XTickLabel',[]) % remove x tick labels
set(gca,'Xtick',[])
xlim([0 1000])
ylabel('\rho, Hz')
maxY = max(get(ax{1},'YLim'));
ylim([0 1.1*maxY])

%%% From spike trains generated above, convolve with kernel to estimate underlying firing rate

sig = 50; % Kernel width, ms
tk = linspace(-1000,1000,10000); % t for kernel, ms
shape = {'boxcar','triangle','epan','gauss'};
types = {'k','k--','k:','k-.'};

% plot the generated spike train
ax{2} = subplot(7,1,2);
stem(ts,ones(numel(ts),1),'k','Marker','none','LineWidth',0.75,'ShowBaseLine','off')
xlim([0 1])
set(gca,'YTickLabel',[]) % remove y tick labels
set(gca,'XTickLabel',[]) % remove x tick labels
set(gca,'YTick',[],'XTick',[])
% title('Example Spike Train')

% set variables to work with code from Initial_development.m
i = 1;
T{i} = ts;
t = tk;

maxYshape = 0;
for h = 1:numel(shape)
    ax{h+2} = subplot(7,1,2+h);

    [K{h},sup{h}] = generateKernels(tk,sig,shape{h});

    % overlay a simple triangular kernel on top of spike train to ensure proper timing
    if numel(T{i}) > 1
        sup_t_idx = t >= sup{h}(1) & t <= sup{h}(2);
        sup_t = t(sup_t_idx);

        frate{h} = 0;
        for j = 1:numel(T{i})
            % get vector for current kernel
            if h == 1 % different if boxcar case
                currK = K{h}*ones(1,numel(sup_t));
            else
                currK = K{h}(sup_t_idx);
            end

            % find closest value in t to current spike time
            spikeloc = find(t > T{i}(j)*1000,1);
            % place it appropriately in time
            idx1 = find(sup_t_idx,1);
            idx3  = find(sup_t_idx,1,'last');
            idx2   = floor((idx1 + idx3)/2);

            fullK = zeros(1,numel(t));
            fullK(idx1 + (spikeloc-idx2) : idx3 + (spikeloc-idx2)) = currK;

            % sum kernels to get estimate of firing rate (in KHz)
            frate{h} = frate{h} + fullK(1:numel(t)); % prevent overflow

            if h == 1
                dt = abs(sup_t(2) - sup_t(1))/2;
                plot([((sup_t(1)-dt)/1000)+T{i}(j), (sup_t/1000)+T{i}(j), ((sup_t(end)+dt)/1000)+T{i}(j)],...
                     [0, K{h}*ones(1,numel(sup_t)), 0],types{h})
            else
                plot((sup_t/1000)+T{i}(j),K{h}(sup_t_idx),types{h});
            end
            hold on
        end
        hold off
        xlim([0 1])
        if max(get(ax{2+h},'YLim')) > maxYshape
            maxYshape = max(get(ax{2+h},'YLim'));
        end
        for xx = 1:h
            set(ax{2+xx},'YLim',[0 1.1*maxYshape])
        end
        set(gca,'YTickLabel',[]) % remove y tick labels
        set(gca,'XTickLabel',[]) % remove x tick labels
        set(gca,'YTick',[],'XTick',[])
    end
end

ax{7} = subplot(7,1,7);
% plot estimate of firing rate
for h = 1:size(frate,2)
    plot(t/1000,frate{h}*1000,types{h})
    xlim([0 1])
    set(gca,'Xtick',[0 0.5 1],'Xticklabel',{'0','0.5','1'})
    xlabel('Time, s')
    ylabel('\lambda, Hz')
    hold on
end
hold off
% legend(shape)
if max(get(ax{7},'YLim')) > maxY
    maxY = max(get(ax{7},'YLim'));
    set(ax{7},'YLim',[0 1.1*maxY])
    set(ax{1},'YLim',[0 1.1*maxY])
else
    set(ax{7},'YLim',[0 1.1*maxY])
end

% Adjust positions and layout of figure
for i = 1:7
    set(ax{i},'Fontname','Times')
    set(ax{i},'Fontsize',10)
    if i == 1
        set(ax{i},'Box','off');
        tmp = get(ax{i},'Position');
        tmp(1,4) = 2*tmp(1,4);
        tmp(1,2) = tmp(1,2) - 0.03;
        set(ax{i},'Position',tmp);
    elseif i > 1 && i < 7
        set(ax{i},'Visible','off')
        if i == 2
            tmp = get(ax{i},'Position');
            tmp(1,4) = tmp(1,4)/2;
            tmp(1,2) = tmp(1,2) + 0.02;
            set(ax{i},'Position',tmp);
        elseif i == 3
            tmp = get(ax{i},'Position');
            tmp(1,4) = tmp(1,4)*2/3;
            tmp(1,2) = tmp(1,2) + 0.05;
            set(ax{i},'Position',tmp);
        elseif i == 4
            tmp = get(ax{i},'Position');
            tmp(1,4) = tmp(1,4)*2/3;
            tmp(1,2) = tmp(1,2) + 0.07;
            set(ax{i},'Position',tmp);
        elseif i == 5
            tmp = get(ax{i},'Position');
            tmp(1,4) = tmp(1,4)*2/3;
            tmp(1,2) = tmp(1,2) + 0.1;
            set(ax{i},'Position',tmp);
        elseif i == 6
            tmp = get(ax{i},'Position');
            tmp(1,4) = tmp(1,4)*2/3;
            tmp(1,2) = tmp(1,2) + 0.12;
            set(ax{i},'Position',tmp);
        end
    elseif i == 7
        set(ax{i},'Box','off');
        tmp = get(ax{i},'Position');
        tmp(1,4) = 2*tmp(1,4);
        tmp(1,2) = tmp(1,2) + 0.03;
        set(ax{i},'Position',tmp);
    end
end


%% Figure 4: Optimizing Kernel Parameters

clearvars;
% close all;

tu_init   = linspace(0,1000,10000); % time vector for underlying rate function (ms)
tro       = 400; % time at response onset (ms)
freq      = [20, 30, 40]; % desired frequency of spiking during onset period
spikes    = 40; % desired number of spikes during onset period
w         = [80, 110, 140]; % response width (determines falling and rising time constants) (ms)
plot_flag = 0; % binary flag for plotting spike train (1 = plot)

for i = 1:numel(freq)
    try
        [err_opt(i), best_shape{i}, sig_opt(i)] = optimizeWidth(tu_init, tro, freq(i),spikes,w(i),plot_flag); %#ok<*SAGROW>
    catch ME
        warning('Spike train generation insufficient...rerunning code')
        close all;
        Final_Figure_Generator; % rerun code if it crashes during spike train generation
    end
end

% Plug sig_opt back in to generate optimal estimates for given firing rates
shape = {'boxcar','triangle','epan','gauss'};

%%% Generate and plot the underlying firing rates
figure;
set(gcf,'Name','Optimized Kernel Parameters');
shade = (192/255)*ones(1,3);
maxY = 0;
for i = 1:numel(freq)
    try
        [ts{i},tu{i},p{i}] = f_generateSpikeTrains(tu_init,tro,freq(i),spikes,w(i),plot_flag); % note: times are returned in ms
    catch ME
        warning('Spike train generation insufficient...rerunning code')
        close all;
        Final_Figure_Generator; % rerun code if it crashes during spike train generation
    end
    ts{i} = ts{i}/1000; % convert spike times to s
    ax{i} = subplot(1,3,i);
    area(tu{i},p{i}*1e3,'EdgeColor',shade,'FaceColor',shade);
    set(gca,'XTickLabel',[]) % remove x tick labels
    set(gca,'Xtick',[])
    xlim([0 1000])
    ylabel('\rho, Hz')
    if max(get(ax{i},'YLim')) > maxY
        maxY = max(get(ax{i},'YLim'));
    end
    for j = 1:i
        set(ax{j},'YLim',[0 maxY*1.1]) % set all ylims to be equivalent
    end
end

%%% From spike trains generated above, convolve with kernel to estimate underlying firing rate

sig = sig_opt; % use optimal widths found above
tk = linspace(-1000,1000,10000); % t for kernel, ms

% figure;
for h = 1:numel(freq)
    % set variables to work with code from Initial_development.m
    i = 1;
    T{h}{i} = ts{h};
    t = tk;
    
    [K{h},sup{h}] = generateKernels(tk,sig(h),best_shape{h});

    % overlay a simple triangular kernel on top of spike train to ensure proper timing
    if numel(T{h}{i}) > 1
        sup_t_idx = t >= sup{h}(1) & t <= sup{h}(2);
        sup_t = t(sup_t_idx);

        frate{i} = 0;
        for j = 1:numel(T{h}{i})
            % get vector for current kernel
            switch best_shape{h}
                case 'boxcar'
                    currK = K{h}*ones(1,numel(sup_t));
                otherwise
                    currK = K{h}(sup_t_idx);
            end

            % find closest value in t to current spike time
            spikeloc = find(t > T{h}{i}(j)*1000,1);
            % place it appropriately in time
            idx1 = find(sup_t_idx,1);
            idx3  = find(sup_t_idx,1,'last');
            idx2   = floor((idx1 + idx3)/2);

            fullK = zeros(1,numel(t));
            fullK(idx1 + (spikeloc-idx2) : idx3 + (spikeloc-idx2)) = currK;

            % sum kernels to get estimate of firing rate (in KHz)
            frate{i} = frate{i} + fullK(1:numel(t)); % prevent overflow
        end
        % plot estimate of firing rate
        subplot(1,3,h)
        hold on
        plot(ax{h},t,frate{i}*1000,'k')
        hold off
        xlabel('Time, s')
        ylabel('\lambda, Hz')
        if max(get(ax{h},'YLim')) > maxY
            maxY = max(get(ax{h},'YLim'));
        end
        for xx = 1:h
            set(ax{xx},'YLim',[0 maxY*1.1]) % set all ylims to be equivalent
        end
    end
end

% Adjust positions and layout of figure
for i = 1:3
    set(ax{i},'Fontname','Times')
    set(ax{i},'Fontsize',10)
    tmp = get(ax{i},'Ylabel');
    if i > 1
        tmp.String = '';
        set(ax{i},'Ylabel',tmp)
    else
        tmp.String = 'Hz';
        set(ax{i},'Ylabel',tmp)
        legend(ax{i},'Underlying Rate','Estimated Rate','Location','northeast')
        legend(ax{i},'boxoff')
    end
    clear tmp
    set(ax{i},'Box','off')
    if i > 1
        set(ax{i},'Ytick',[])
    end
    set(ax{i},'Xtick',[0 500 1000],'Xticklabel',{'0','0.5','1'})
end