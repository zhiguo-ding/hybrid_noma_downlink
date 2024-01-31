clear all;
%figure
M = 10;
R=1;
Rvec = [1:0.5:5];
ct = 2000;
for mx = 1 : length(Rvec)
    R = Rvec(mx);
    bad_count=0;
    for ict = 1 : ct
        hm = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1)); 
        tempm = ones(1,M);
        hm = kron(hm,tempm);
        hm = abs(hm).^2;
        %hm = sort(hm,'descend');
        test_bad = 0;
        for m = 1 : M
            for i = 1 : m
                hbar(m,i) = min(hm(i:m,i));
                if hbar(m,i) <0.01
                    test_bad=1;
                end
            end
        end
        if test_bad==1
            bad_count = bad_count+1;
            continue; 
        end
        
        Pij = zeros(M,M);
        Pij(1,1) = (2^R-1)/hm(1);
        for m =2 : M
            bmi=[];
            for i = 1 : m
                bmi(i) = 1/(1+hbar(m,i)*sum(Pij(i:m-1,i)));
            end

            A = []; % No other constraints
            b = [];
            Aeq = [];%
            beq = [];%zeros(M+1,1); 
            lb = [];
            ub = [];
            x0 =  ones(m,1); %initilization 
            options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
            x = fmincon(@(x) sum(x),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,m,bmi,hbar,R),options);
            Pij(m,1:m) = x';

        end
        % 
        % Pija = zeros(M,M);
        % Pija(1,1) = (2^R-1)/hm(1,1);
        % for m = 2 : M
        %     for i = 1 : m                
        %         prod1 = 1;
        %         for p =1 : m
        %             prod1 = prod1 * hbar(m,p)/(hbar(m,p)*sum(Pija(p:m-1,p))+1);
        %         end
        % 
        %         Pija(m,i) = (exp(R)/prod1)^(1/m) - ...
        %             (hbar(m,i)*sum(Pija(i:m-1,i)) +1 )/hbar(m,i);
        % 
        % 
        %     end
        % end    
        % 
        Eallt(ict) = sum(sum(Pij));
        for m = 1 : M
            Poma(m) = (exp(R)-1)/hm(m);
        end
        Eallomat(ict) = sum(Poma);    
    end
    Eall(mx) =  sum(Eallt)/(ct-bad_count);
    Ealloma(mx) =  sum(Eallomat)/(ct-bad_count);
end
plot(Rvec,Ealloma,Rvec, Eall)
    

function [c,ceq] = mycons(x,m,bmi,hbar,R)

sum1 = 0;
for i =1 : m
    sum1 = sum1 + log(1+bmi(i)*hbar(m,i)*x(i));
end

c(1,1) = R-sum1;
c(2:m+1,1) = -x;
ceq = [];
end
