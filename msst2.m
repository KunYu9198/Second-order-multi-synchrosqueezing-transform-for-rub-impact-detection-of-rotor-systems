function [STFT,SST1,SST2,omega,omega2] = msst2(s,gamma,sigma,num)
% msst2 : multi-second synchrosqueezing transform
%   Uses a Gaussian window.
%
% INPUTS:
%   s: real or complex signal, must be of length 2^N
%   gamma: threshold
%   sigma: window parameter
%   ft: frequency bins
%   bt: time bins
%
% OUTPUTS:
%   STFT: the short-time Fourier transform
%   SST2: vertical second-order synchrosqueezing
%   omega2: second-order instantaneous frequency

% checking length of signal
n = length(s);
s = hilbert(s);
s = s(:);

% Optional parameters

ft = 1:n/2;    %频率轴长度
bt = 1:n;      %时间轴长度

% Padding
sleft = flipud((s(2:n/2+1)));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];

%% STFT and operators tau and omega
nb = length(bt);
neta = length(ft);

t = -0.5:1/n:0.5-1/n;t=t';

g =  1/sigma*exp(-pi/sigma^2*t.^2);
gp = -2*pi/sigma^2*t .* g; % g'

% 从这里开始改
STFT = zeros(neta,nb);  %STFT的结果

omega = zeros(neta,nb); %频率扩散范围
tau2 = zeros(neta,nb);

omega2 = zeros(neta,nb); %2阶频率扩散范围

phi22p = zeros(neta,nb);

vg = zeros(neta,7);
vgp = zeros(neta,5);
Y = zeros(neta,4,4);


for b=1:nb
    % 进行了3次STFT, window t^n*g
    for i=0:2
        tmp = (fft(x(bt(b):bt(b)+n-1).*(t.^i).*g))/n;
        vg(:,i+1) = tmp(ft);
    end
    % 进行了2次STFT, window t^n*g
    for i=0:1
        tmp = fft(x(bt(b):bt(b)+n-1).*(t.^i).*gp)/n;
        vgp(:,i+1) = tmp(ft);
    end
    % second-order operator tau
    tau2(:,b) = vg(:,2)./vg(:,1);    % 要
    
    
    % Y expressions
    for i = 1:2
        for j = 1:2
            if i>=j
                Y(:,i,j) = vg(:,1).*vg(:,i+1) - vg(:,j).*vg(:,i-j+2);
            end
        end
    end
    
    % W expressions
    W2 = 1/2/1i/pi*(vg(:,1).^2+vg(:,1).*vgp(:,2)-vg(:,2).*vgp(:,1));  % 要
    
    % operator omega
    omega(:,b) = (ft-1)'-real(vgp(:,1)/2/1i/pi./vg(:,1));
    % operator hat p: estimations of frequency modulation
    % SST2
    phi22p(:,b) = W2./Y(:,2,2);  % 要
    omega2(:,b) = omega(:,b) + real(phi22p(:,b).*tau2(:,b));
     
    % Storing STFT
    STFT(:,b) = vg(:,1).* exp(1i*pi*(ft-1)'); % compensates the tranlation 1/2 of s
   
end

tfr_sst1=STFT;
tfr_sst2=STFT;

for ite=1:num
    ite
    SST1 = zeros(neta,nb);
    SST2 = zeros(neta,nb);
    for b=1:nb
        for eta=1:neta
            if abs(STFT(eta,b))>0.001*gamma
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%SST1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                k = 1+round(omega(eta,b));
                if k>=1 && k<=neta
                    SST1(k,b)  = SST1(k,b) + tfr_sst1(eta,b);
                end
                
                %%%%%%%%%%%%%%%%%%%%SST2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                k = 1+round(omega2(eta,b));
                if k>=1 && k<=neta
                    SST2(k,b)  = SST2(k,b) + tfr_sst2(eta,b);
                end
                
            end
        end
    end
    tfr_sst1=SST1;
    tfr_sst2=SST2;
end


SST1 = sigma*SST1;
SST2 = sigma*SST2;







