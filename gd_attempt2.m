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
Nex = 30; % number of training examples to generate

tin       = linspace(0,1000,10000); % time vector for underlying rate function (ms)
tro       = 400; % time at response onset (ms)
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

ep = 0.01; % learning rate

% Define indices for kernel shape
% 1 = Boxcar
% 2 = Triangle
% 3 = Epanechnikov
% 4 = Gaussian
shape = {'boxcar','triangle','epan','gauss'};

% Shuffle the training data
new_idx = randperm(Nex);

% Randomly initialize weight for kernel width between 1 and 100
weights = [rand;rand];
weights = repmat(weights,1,4)';
store_weights = weights;
% sig = (round((rand*99))+1)*ones(numel(shape),1); % (ms); dimension is due to number of kernel shapes to explore
sig =0.01*ones(numel(shape),1);
% Generate estimate of firing rate for each kernel shape
tk = linspace(-1000,1000,10000); % t for kernel, ms
t = tk;
all_grad =[0;0;0;0];
for h = 1:Nex
    for i = 1:numel(sig)
        sig_in = weights(:,1)+weights(:,2)'*sig;
        [K{i},sup{i}] = generateKernels(tk,sig_in(i),shape{i});
        sup_t_idx = t >= sup{i}(1) & t <= sup{i}(2);
        sup_t = t(sup_t_idx);
        frate{i} = 0;
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
        end
        clear sup_t_idx sup_t currK spikeloc idx1 idx2 idx3 fullK
    end
    
    % Calculate error and gradient
    frate = vertcat(frate{:});
    err_new = 2*sum(((frate - p).^2),2);
    all_err(:,h) = err_new;
    gradient = sum(((frate-p)),2);
    all_grad(:,h+1) = gradient;
    

    % update weights
    
    weights(:,1) = weights(:,1) - ep*gradient;
    weights(:,2) = weights(:,2) - ep*gradient.*sig;
    
%     weights(:,1) = weights(:,1) - ep*(all_grad(:,h+1)-all_grad(:,h));
%     weights(:,2) = weights(:,2) - ep*(all_grad(:,h+1)-all_grad(:,h)).*sig;
    store_weights(:,:,h) = weights;
    clear err_new frate
    
end

[M1, I1] = min([min(all_err(1,:)),min(all_err(2,:)),min(all_err(3,:)),min(all_err(4,:))]);
[M2, I2] = min(all_err(I1,:));

weights_opt = store_weights(:,:,I2);
weights_opt = weights_opt(I1,:);
err_opt = M1
best_shape = shape{I1}
sig_opt = weights_opt(1) + weights_opt(2)*sig(I1)


figure; plot(all_err','LineWidth',1); xlabel('Iterations'); ylabel('ISE'); legend('boxcar','triangle','epan','gauss')
figure; plot(all_grad','LineWidth',1); xlabel('Iterations'); ylabel('ISE gradient'); legend('boxcar','triangle','epan','gauss')