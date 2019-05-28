%constants
w_c = 3.5;
w_q = 12;
w_mm = 0.2;
w_ff = 0.2;
w_n = 1.05;
v = 3.2;
ro = 0.68;
q_c = 7.6;
x = 0.8;
beta = 0.3;
e_c = 3.8;
alfa = 0.4;


Iter=100;

%LOWHouseholds
f = 1;
b = 0.5;
d = 1;

T_e = beta*v*(w_mm+w_ff)/(w_n-beta*(w_mm+w_ff));

e_m=zeros(1,Iter);
e_f=zeros(1,Iter);
N=zeros(1,Iter);
K=zeros(1,Iter);

e_m(1) = 8;
e_f(1)= 6.7;
N(1) = 1;
K(1)=0.1;


h=0.3*ones(1,Iter);
n=3*ones(1,Iter);


for i=1:Iter
    
    for j=1:10
        Y(i) = ((0.5*e_m(i)*(1-f*h(i))*N(i))^alfa)*((0.5*e_f(i)*(1-(h(i)+v*n(i)+T_e*n(i)))*N(i))^alfa)*(K(i))^(1-2*alfa);
        w_m(i) = (alfa*Y(i))/(e_m(i)*(1-f*h(i))*(N(i)/2));
        w_f(i) = (alfa*Y(i))/(e_f(i)*(1-(h(i)+v*n(i)+T_e*n(i))*(N(i)/2)));
        r(i) = (1-2*alfa)*Y(i)/K(i);
        c(i) = (e_m(i)*w_m(i)+d*e_f(i)*w_f(i))/((1+((1+r(i))*w_c)/(1+ro))+2*w_q/w_c);
        h(i) = (e_m(i)*w_m(i)+d*e_f(i)*w_f(i))/((w_c/w_q)*(1+((1+r(i))*w_c)/(1+ro))+2)*(f*e_m(i)*w_m(i)+d*e_f(i)*w_f(i));
        n(i)=((e_m(i)*w_m(i)+d*e_f(i)*w_f(i))/((1+(((1+r(i))*w_c)/(1+ro)))+2*w_q/w_c))   *    ((w_n-beta*(w_mm+w_ff))/(v*w_c*e_f(i)*w_f(i)));
    end
    
    K(i+1) = 0.5*N(i)*((w_c/(ro+1))*c(i));
    N(i+1) = n(i)*N(i)/2;
    e_m(i+1) = e_c*((e_m(i))^(1-beta))*(2*b*T_e)^beta;
    e_f(i+1) = e_c*((e_m(i))^(1-beta))*(2*(1-b)*T_e)^beta;
end





plot((1:20),Y(1:20));


