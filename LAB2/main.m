fs=1000;
ts=1/fs;
T=2;
t=-T:ts:T;


fm=1;
Tm=1/fm;
m=square(2*pi*fm*t);

k_harmonics=3;

k_vec=-k_harmonics:k_harmonics;
cnt=1;
m_approx=0;
while(cnt<=length(k_vec))
    k=k_vec(cnt);
    fun1=@(x) 1*exp(-2*j*pi*k*fm*x);
    fun2=@(x) -1*exp(-2*j*pi*k*fm*x);
    q1=integral(fun1,0,Tm/2);
    q2=integral(fun2,Tm/2,Tm);
    FC(cnt)=q1+q2;
    m_approx=m_approx+FC(cnt)*exp(2*pi*j*k*fm*t);
    cnt=cnt+1;
end
figure(1);
plot(t,m);
hold on;
plot(t,m_approx);


