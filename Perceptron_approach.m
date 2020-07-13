% Taylor Hansen
% Mitch Thomas
% BIOEN 6005 Project

clearvars;
close all;

%% Framework for optimizing kernel parameters

% Generate all training examples

% Randomly initialize the weight for kernel width
% Should be a random value in a range

% Shuffle the training examples

% for example in all examples
%     Generate an estimate of the firing rate for current kernel width for each kernel shape
%     Check against ground truth
%     Check error direction against previous error from last example
%     Update weights in appropriate direction
% end

%% Generate training examples
Nex = 1000; % number of training examples to generate
tspan = [0 1000]; % time over which to run simulations (ms)

tin       = linspace(tspan(1),tspan(2),10000); % time vector for underlying rate function (ms)
tro       = 50; % time at response onset (ms)
freq      = 20; % desired frequency of spiking during onset period
spikes    = 20; % desired number of spikes during onset period
w         = 100; % response width (determines falling and rising time constants) (ms)
plot_flag = 0; % binary flag for plotting spike train (1 = plot)

strains = cell(Nex,1);
[strains{1},tu,p] = f_generateSpikeTrains(tin,tro,freq,spikes,w,plot_flag); % note: times are returned in ms
strains{1} = strains{1}/1000; % convert spike times to s
for i = 2:Nex
    [strains{i},~,~] = f_generateSpikeTrains(tin,tro,freq,spikes,w,plot_flag); % note: times are returned in ms
    strains{i} = strains{i}/1000; % convert spike times to s
end

%% Train perceptron

ep  = 1; % learning rate
Nst = 1; % number of spike trains to average

% Define indices for kernel shape
% 1 = Boxcar
% 2 = Triangle
% 3 = Epanechnikov
% 4 = Gaussian
shape = {'boxcar','triangle','epan','gauss'};

% Shuffle the training data
new_idx = randperm(Nex);

% Randomly initialize weight for kernel width between 1 and 100
sig = (round((rand*99))+1)*ones(numel(shape),1); % (ms); dimension is due to number of kernel shapes to explore
all_sig = sig;

% Generate estimate of firing rate for each kernel shape
tk = linspace(-tspan(2),tspan(2),10000); % t for kernel, ms
t = tk;

counter = 1;

for h = 1:Nex
%     ep = ep/h;
%     disp(ep)
    for i = 1:numel(sig)
        [K{i},sup{i}] = generateKernels(tk,sig(i),shape{i});
        sup_t_idx = t >= sup{i}(1) & t <= sup{i}(2);
        sup_t = t(sup_t_idx);
        frate{i} = 0;
        bin_st = zeros(1,numel(tin));
        for j = 1:numel(strains{h})
            % get vector for current kernel
            if i == 1 % different if boxcar case
                currK = K{i}*ones(1,numel(sup_t));
            else
                currK = K{i}(sup_t_idx);
            end

            % find closest value in t to current spike time
            spikeloc = find(t > strains{h}(j)*1000,1);
            % place it appropriately in time
            idx1 = find(sup_t_idx,1);
            idx3  = find(sup_t_idx,1,'last');
            idx2   = floor((idx1 + idx3)/2);

            fullK = zeros(1,numel(t));
            fullK(idx1 + (spikeloc-idx2) : idx3 + (spikeloc-idx2)) = currK;

            % sum kernels to get estimate of firing rate (in KHz)
            frate{i} = frate{i} + fullK(1:numel(t)); % prevent overflow
            
            % Generate binary train to represent spike train according to tin
            bin_loc = find(tin > strains{h}(j)*1000,1);
            bin_st(bin_loc) = 1;
            clear bin_loc
        end
        clear sup_t_idx sup_t currK spikeloc idx1 idx2 idx3 fullK
    end
    
    % Average frate across N number of spike trains
    if mod(h-1,Nst) == 0
        clear all_frate
        all_frate = frate{1};
    else
        all_frate = vertcat(all_frate,frate{1});
    end
    
    if Nst > 1
    
        % Calculate error and push sigs in correct direction
        if mod(h,Nst) == 0 && h > 1
            frate = mean(all_frate,1);
            err_new = sum(((frate - p).^2),2);

            all_err(:,counter) = err_new; % store all errors

            if h == Nst
                err = err_new;
                sig = sig + ep*10;
            else
                traj = err - err_new; % there's probably a fancier way to check the sign
                for idx = 1:numel(traj)
                    if traj(idx) < 0
                        sig(i) = sig(i) + ep*traj(idx);
                    else
                        sig(i) = sig(i) - ep*traj(idx);
                    end
                end
            end
            counter = counter + 1;
            clear err_new frate
        end
    else
              
        % Calculate error and push sigs in correct direction
        frate = vertcat(frate{:});
        err_new = sum(((frate - p).^2),2);

        all_err(:,h) = err_new; % store all errors
        
        % Update sigma appropriately
        grad = (frate - p)*bin_st';
        
        sig = sig - ep*grad;
        all_sig = horzcat(all_sig,sig);
        
        
%         if h == 1
%             err = err_new;
%             sig = sig + ep*10;
%         else
%             traj = err - err_new; % there's probably a fancier way to check the sign
%             for idx = 1:numel(traj)
%                 if traj(idx) < 0
%                     sig(i) = sig(i) + ep*traj(idx);
%                 else
%                     sig(i) = sig(i) - ep*traj(idx);
%                 end
%             end
%         end
        clear err_new frate
    end
end
figure; plot(all_err','LineWidth',1); xlabel('Iterations'); ylabel('ISE');
if Nst > 1
    legend('boxcar')
else
    legend('boxcar','triangle','epan','gauss')
end