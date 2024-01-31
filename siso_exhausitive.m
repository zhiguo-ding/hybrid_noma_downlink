clear all;
figure
M = 3;
hm = [1:M];
hm = sort(hm,'descend');
R=2;
Pij = zeros(M,M);
Pij(1,1) = (2^R-1)/hm(1);
 

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

Eall = sum(sum(Pija)) 

%exhaustive search
Pvec = [0:0.05:2];
for i =1 : length(Pvec)
    for j = 1 : length(Pvec)
        P31 = Pvec(i); P32 = Pvec(j);
        R31 = log(1+P31*hm(3)/(hm(3)*(Pija(1,1)+Pija(2,1)) +1));
        R32 = log(1+P32*hm(3)/(hm(3)*(Pija(2,2)) +1));
        P33 = (exp(R-R31-R32)-1)/hm(3);
        Eex(i,j) = P31 + P32 + P33;
        Enoma(i,j) = sum(Pija(3,:));
    end
end

 

surf(Pvec,Pvec, Enoma )
alpha 0.5
hold

surf(Pvec,Pvec, Eex )
