clear all;
figure
M = 5;
R=1;
Rvec = [1:0.5:5];
ct = 5000;
for mx = 1 : length(Rvec)
    R = Rvec(mx);
    bad_count=0;
    for ict = 1 : ct
        hm = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));   %m  
        hm = abs(hm).^2*10;
        hm = sort(hm,'descend');
        if min(hm)<0.01
            bad_count = bad_count+1;
            continue; 
        end
        
        Pij = zeros(M,M);
        Pij(1,1) = (2^R-1)/hm(1);
        for m =2 : M
            bmi=[];
            for i = 1 : m
                bmi(i) = 1/(1+hm(m)*sum(Pij(i:m-1,i)));
            end

            A = []; % No other constraints
            b = [];
            Aeq = [];%
            beq = [];%zeros(M+1,1); 
            lb = [];
            ub = [];
            x0 =  ones(m,1); %initilization 
            options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
            x = fmincon(@(x) sum(x),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,m,bmi,hm,R),options);
            Pij(m,1:m) = x';

        end
        
        Pija = zeros(M,M);
        Pija(1,1) = (2^R-1)/hm(1);
        for m = 2 : M
            for i = 1 : m
                prod1 = 1;
                for p =1 : m
                    prod1 = prod1 * hm(m)/(hm(m)*sum(Pija(p:m-1,p))+1);
                end
        
                Pija(m,i) = (exp(R)/prod1)^(1/m) - (hm(m)*sum(Pija(i:m-1,i)) +1 )/hm(m);
        
        
            end
        end    
        
        Eallt(ict) = sum(sum(Pija));
        Eallsimt(ict) = sum(sum(Pij));
        for m = 1 : M
            Poma(m) = (exp(R)-1)/hm(m);
        end
        Eallomat(ict) = sum(Poma);    
    end
    Eall(mx) =  sum(Eallt)/(ct-bad_count);
    Eallsim(mx) =  sum(Eallsimt)/(ct-bad_count);
    Ealloma(mx) =  sum(Eallomat)/(ct-bad_count);
end
plot(Rvec,Ealloma,Rvec, Eallsim,Rvec, Eall)
    

function [c,ceq] = mycons(x,m,bmi,hm,R)

sum1 = 0;
for i =1 : m
    sum1 = sum1 + log(1+bmi(i)*hm(m)*x(i));
end

c(1,1) = R-sum1;
c(2:m+1,1) = -x;
ceq = [];
end
