clear all

N =513; % number of antennas
%Kpu = 30;% the number of PU
K=3; %10 users are served by SDMA
M=20; % the number of PU
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % wavelength
d = lambda / 2; % antenna spacing
D = (N-1)*d; %apture or size
dis_Ray = 2*D^2/lambda;
dis_Fre = sqrt(D^3/lambda)/2;
noise = 10^(-100/10); %noise power
Ppu =0.01;% near field users' power
R=K*5;
ict= 200;
bf_beam = 0; %0 for ZF, 1 for beamfocusing 
dyna_location = 0;% 1 for dynamic locations and 0 for fixed
theta_step = pi/40;

radiusPU = [50 100];%[50 dis_Ray/3];

Rvec =  [1:7];
for mx =  1: length(Rvec)
     R = K*Rvec(mx);
    for ictx = 1 : 1     
        %the location of the near-field users [5 70]
        H = [];
        thetall = [];
        %PU
            theta = [0:theta_step:theta_step*(M-1)];% theta=thetax(1:M);        
            r = ones(M,1)*50;
            %theta = pi/4*ones(M,1);
            %r = [0 : 5 : (M-1)*5] + 100;
            
            for k = 1 : M
                wch(:,k) = beamfocusing(r(k), theta(k), N, d, f)/sqrt(N); % beamforming vector
                alpham(k) = lambda/4/pi/r(k);
                H(:,k) = sqrt(N)*alpham(k)*wch(:,k);
            end 
            % either ZF or BF can be used, decided by bf_beam
            invH = inv(H'*H);
            D = diag(sqrt(1./diag(invH)));
            WPU = H*invH*D;
            if bf_beam == 1 % decide which type of beamforming, ZF or BF
                WPU = wch; %beamfocusing 
            end
         
            %SU
            %theta_SUx = [-pi/4:pi/K/2:pi/4]; theta_SU=theta_SUx(1:K);             
            theta_SU = [theta(end-2:end)  ];%[theta(1) theta(2) theta(3)];
            r_SU = ones(K,1)*200; 
        
            for k = 1 : K
                wchSU(:,k) = beamfocusing(r_SU(k), theta_SU(k), N, d, f)/sqrt(N); % beamforming vector
                alpham(k) = lambda/4/pi/r_SU(k);
                G(:,k) = sqrt(N)*alpham(k)*wchSU(:,k);
            end
            % for SU, we alwasy use ZF, as they can be FF and BF is not applicable
            invG = inv(G'*G);
            DG = diag(sqrt(1./diag(invG)));
            ZF_SU = G*invG*DG;    
            WSU=ZF_SU;
            if bf_beam == 1 %beamfocusing is used,
                WSU = wchSU;
            end
            WSU_oma = WSU;% 
        
        
        % beam scheduling
        ik=[];
        theta_temp = theta;        WPUt = WPU;
        if bf_beam == 1 %beamfocusing is used, then just focus on the angel
            for k = 1 : K
                %[z1,ik(k)] = max(abs(G(:,k)'*WPUt).^2);
                %WPUt(:,ik(k)) = 0;% exclude this to be chosen by other users
                tempa1 = theta_SU(k) - theta_temp;
                [x,ik(k)] = min(abs(tempa1));
                theta_temp(ik(k)) = 10000;% exclude this to be chosen by other users
            end
        else  %ZF CASE
            for k = 1 : K
                [z1,ik(k)] = max(abs(G(:,k)'*WPUt).^2);
                WPUt(:,ik(k)) = 0;% exclude this to be chosen by other users
            end
        end
        
        % the channels in the paper
        for k = 1 : K
            for j =  1 : K
                gkj(k,j) = abs(G(:,k)'*WPU(:,ik(j)))^2/noise;
            end
            bk(k) = Ppu*(sum(abs(G(:,k)'*WPU).^2))/noise + 1;
            for m = 1 : K
                hkm(k,m) = abs(H(:,ik(k))'*WPU(:,ik(m)))^2/noise;
            end
            dk(k) =  Ppu*sum(abs(H(:,ik(k))'*WPU).^2)/noise+1;
            for i = 1 : K
                cki(k,i) = abs(G(:,k)'*WSU(:,i))^2/noise;
                cki_oma(k,i) = abs(G(:,k)'*WSU_oma(:,i))^2/noise;
            end
             
        end 
        
        %for OMA    
        A = []; % No other constraints
        b = [];
        Aeq = [];%
        beq = [];%zeros(M+1,1); 
        lb = [];
        ub = [];
        xoma0 =  ones(K,1); %initilization 
        options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
        xoma = fmincon(@(xoma) sum(xoma(1:K)),xoma0,A,b,Aeq,beq,lb,ub,...
            @(xoma) myconsoma(xoma,K,noise,cki_oma,R),options);
        
        if min(xoma)<0
            dd=1
        end
        %for NOMA
        for scai = 1 : 5
          A = []; % No other constraints
            b = [];
            Aeq = [];%
            beq = [];%zeros(M+1,1); 
            lb = [];
            ub = [];
            x0 =  [zeros(2*K,1)];%;xoma]; %initilization 
            options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
            if scai ==1
                f0=zeros(K,1); e0= zeros(K,1);
            else
                f0 = x(1:K); e0 = x(K+1:2*K);
            end
            x = fmincon(@(x) sum(x(1:K))*M+sum(x(K+1:end))*K,x0,A,b,Aeq,beq,lb,ub,...
                @(x) mycons(x,M, K, gkj, bk,f0,e0,noise, hkm, dk,cki,R)  ,options);
            %feasibility(x,M, K, gkj, bk,f0,e0,noise, hkm, dk,cki,R)
            %sum(x(1:K))*M+sum(x(K+1:end))*K
        end
        feasibility(x,M, K, gkj, bk,f0,e0,noise, hkm, dk,cki,R)
        powerNx(ictx) = sum(x(1:K))*M+sum(x(K+1:end))*K;
        powerOMAz(ictx) = sum(xoma)*K;
    end
    powerN(mx) = mean(powerNx);
    if min(xoma)<0
        powerOMA(mx) = inf;
    else
        powerOMA(mx) = mean(powerOMAz);
    end
end


%plot(Mvec, powerOMA, Mvec, powerN)
plot(Rvec, powerOMA, Rvec, powerN)
ylim([0 10])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [c,ceq] = mycons(x,M, K, gkj, bk,f0,e0,noise, hkm, dk,cki,R)  

f = x(1:K); e = x(K+1:2*K);
for k = 1:K
    hk = [hkm(k,:)]';
    hktilde = hk; hktilde(k) = 0;
    ck = [cki(k,:)]';
    cktilde = ck; cktilde(k) = 0;
    c(k,1) = R - M*log(hk'*f+dk(k)) + M*log(hktilde'*f0+dk(k)) ...
        -M/(hktilde'*f0+dk(k))*hktilde'*(f-f0)...
        - K*log(ck'*e+1) + K*log(cktilde'*e0+1) ...
        -K/(cktilde'*e0+1)*cktilde'*(e-e0);
end


for k = 1:K
    gk = [gkj(k,:)]';
    gktilde = gk; gktilde(k) = 0;
    ck = [cki(k,:)]';
    cktilde = ck; cktilde(k) = 0;
    c(K+k,1) = R - M*log(gk'*f+bk(k)) + M*log(gktilde'*f0+bk(k)) ...
        -M/(gktilde'*f0+bk(k))*gktilde'*(f-f0)...
        - K*log(ck'*e+1) + K*log(cktilde'*e0+1) ...
        -K/(cktilde'*e0+1)*cktilde'*(e-e0);
end
 
c(2*K+1:4*K) = -x;
ceq = [];
end   

function [c,ceq] = myconsoma(x,K,noise,cki,R)  

for k = 1 : K 
    ck = [cki(k,:)]';
    cktilde = ck; cktilde(k) = 0;
    gammax = exp(R/K)-1;
    c(k,1) =  gammax + gammax*cktilde'*x-cki(k,k)*x(k); 
end
c(K+1:2*K) = -x;
ceq = [];
end

%%%%%%%% to check the feasibility of the solution
function [c,ceq] = feasibility(x,M, K, gkj, bk,f0,e0,noise, hkm, dk,cki,R)   

f = x(1:K); e = x(K+1:2*K);
for k = 1:K
    hk = [hkm(k,:)]';
    hktilde = hk; hktilde(k) = 0;
    ck = [cki(k,:)]';
    cktilde = ck; cktilde(k) = 0;
    c(k,1) = R - M*log(hk'*f+dk(k)) + M*log(hktilde'*f+dk(k)) ...
        - K*log(ck'*e+1) + K*log(cktilde'*e+1);
end


for k = 1:K
    gk = [gkj(k,:)]';
    gktilde = gk; gktilde(k) = 0;
    ck = [cki(k,:)]';
    cktilde = ck; cktilde(k) = 0;
    c(K+k,1) = R - M*log(gk'*f+bk(k)) + M*log(gktilde'*f+bk(k)) ... 
        - K*log(ck'*e+1) + K*log(cktilde'*e+1);
end
for k = 1:K
    hk = [hkm(k,:)]';
    hktilde = hk; hktilde(k) = 0;
    ck = [cki(k,:)]';
    cktilde = ck; cktilde(k) = 0;
    c(2*K+k,1) =   M*log((hk'*f+dk(k))/(hktilde'*f+dk(k)))  ...
        + K*log((ck'*e+1)/(cktilde'*e+1));
end


for k = 1:K
    gk = [gkj(k,:)]';
    gktilde = gk; gktilde(k) = 0;
    ck = [cki(k,:)]';
    cktilde = ck; cktilde(k) = 0;
    c(3*K+k,1) =   M*log((gk'*f+bk(k))/(gktilde'*f+bk(k))) ... 
        + K*log((ck'*e+1)/(cktilde'*e+1));
end

end