clearvars
close all

%% Generate underlying rate function(s)

% Pulled from Nawrot paper example
tu = linspace(0,1000,10000); % t for kernel, ms
tf = tu(end);
tro = 400; % ms, time at response onset
b = 20/1e3; % MHz
A = 20; % spikes
w = 100; % ms, where w = sqrt(tau1^2 + tau2^2)
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

figure;
subplot(2,1,1)
plot(tu_new,p_new*1e3,'k')
xlabel('Time, ms')
ylabel('\rho, Hz')

%% Generate the spike trains from the underlying rate function(s) above

tf = 1000; % ms, time to stop simulation
t0 = 0;
n = 0;
while t0<tf
    n = n + 1;
    if t0 == 0
        lam = b; % essentially the baseline of p_new
    else
        lam = p_new(find(tu_new > t0,1));
    end
    t0 = t0 - (1/lam*log(1-rand));
    if t0 < tf
        t(n) = t0;
    end
end
subplot(2,1,2)
% figure;
stem(t,ones(numel(t),1),'k','Marker','none','LineWidth',1)
xlim([0 1000])
set(gca,'YTickLabel',[]) % remove y tick labels
set(gca,'XTickLabel',[]) % remove x tick labels
set(gca,'YTick',[],'XTick',[])

%% Chuck's PS3.4a
% clearvars;
% % close all;
% 
% R = 5;
% tf = 1000; % ms
% lam = [.003 .01 .03 .1]; % firing rates
% Nl = length(lam);
% for k = 1:Nl
%     for r = 1:R % reps
%         t0 = 0;
%         n = 0;
%         while (t0<tf)
%             n = n + 1;
%             t0 = t0 - (1/lam(k)*log(1-rand));
%             if t0 < tf
%                 t{k,r}(n) = t0;
%             end
%         end
%     end
% end
% 
% figure;
% pd = {'+r','+m','+b','+k'};
% plot([0 1000], [5 5; 10 10; 15 15], 'k:');
% hold on;
% for k=1:Nl
%     for r = 1:R
%         plot(t{k,r},(k-1)*R+r-0.5,pd{k});
%     end
% end
% hold off;
% set(gca,'YTick',[0.5:1:19.5],'YTickLabel','');
% for k = 1:4
%     txt(k) = text(-50,1+5*(k-1),...
%         ['\lambda = ',num2str(lam(k)*1e3),' Hz']);
% end
% set(txt,'Rotation',90);
% xlabel('Time, ms');
% title('PS3, Prob 4a');