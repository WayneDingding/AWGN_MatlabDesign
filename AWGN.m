%AWGN generator
clear all 
clc
load('C:\Users\dingw_000\Desktop\Verilog Design\new_urng.mat');
%input interval (xa)
B_x_e_a=8; 
B_x_f_a=6;
B_x_g_aa=7;
B_x_g_ba=7;
%%%%%%%%%%%calculate values of all original BCi%%%%%%%%%%%%%
%coefficients for e
xe=[1:2^-8:2-2^-8;1+2^-9:2^-8:2-2^-9;1+2^-8:2^-8:2]';
ye=-log(xe);
pe=zeros(size(xe));
for i=1:size(xe,1)
pe(i,1:3)=polyfit(xe(i,:),ye(i,:),2);
end
C2_e=floor(pe(:,1)*2^13)/2^13;
C1_e=floor(pe(:,2)*2^22)/2^22;
C0_e=floor(pe(:,3)*2^30)/2^30;
%coefficients for f
xf=[1:2^-6:2-2^-6,2:2^-5:4-2^-5;1+2^-6:2^-6:2,2+2^-5:2^-5:4;]';
yf=sqrt(xf);
pf=zeros(size(xf));
for i=1:size(xf,1)
pf(i,1:2)=polyfit(xf(i,:),yf(i,:),1);
end
C1_f=floor(pf(:,1)*2^12)/2^12;
C0_f=floor(pf(:,2)*2^20)/2^20;
%coefficients for sin/cos
xg=[0:2^-7:1-2^-7;2^-7:2^-7:1]';
yg=sin(pi/2*xg);
pg=zeros(size(xg));
for i=1:size(xg,1)
pg(i,1:2)=polyfit(xg(i,:),yg(i,:),1);
end
C1_g=floor(pg(:,1)*2^18)/2^18;
C0_g=floor(pg(:,2)*2^18)/2^18;
%%%%%%%%%%% Generate u0 & u1 %%%%%%%%
a=taus(urng_seed(:,1), urng_seed(:,2), urng_seed(:,3));
b=taus(urng_seed(:,4), urng_seed(:,5), urng_seed(:,6));
u0=(a*2^16+floor(b/2^16));% u0 & u1 must be over interval [0,1)
u1=bitand(b,65535);
%if input u0 is 0, make output 0. 
if u0==0
    x0=0;
    x1=0;
else 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Evaluate e=-2ln(u0) %%%%%%%%%%%
%range reduction
exp_e =LZD(u0)+1;
x_e   =u0.*2.^exp_e;
%Approximate
x_e_xa =floor(x_e/2^(48-B_x_e_a))/2^(B_x_e_a);
x_e_xb =(x_e-x_e_xa*2^48)/2^(48-B_x_e_a);%%
xexa=x_e_xa*2^(B_x_e_a)-256;%Coefficient index

C22_e =C2_e(xexa+1)/(4^B_x_e_a);
C11_e=(2*C2_e(xexa+1).*x_e_xa+C1_e(xexa+1))/(2^B_x_e_a);    
C00_e=C2_e(xexa+1).*x_e_xa.^2+C1_e(xexa+1).*x_e_xa+C0_e(xexa+1);
y_e=floor((floor((floor(C22_e.*x_e_xb*2^30)/2^30+C11_e).*x_e_xb*2^30)/2^30+C00_e)*2^27)/2^27;
%range reconstruction
ln2   =floor(log(2)*2^32)/2^32;
ee    =floor(exp_e*ln2*2^28)/2^28;
e     =floor((ee+y_e)*2*2^24)/2^24;
if e<0
    e=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Evaluate f=sqrt(e) %%%%%%%%%%%%%
e_f=e*2^(41);
%range reduction
exp_f =5-LZD(e_f);
x_ff  =e_f./2.^(exp_f);
x_f=x_ff./(2.^mod(exp_f,2));
%Approximate
x_f_xa =floor(x_f/2^(48-7-B_x_f_a))/2^(B_x_f_a);
x_f_xb =(x_f-x_f_xa*2^41)/2^(48-7-B_x_f_a);%%
xfxa=(x_f_xa*2^(B_x_f_a)-64).*mod(exp_f,2)+...
      (floor((x_f_xa*2^(B_x_f_a)-128)/2)+64 ).*(1-mod(exp_f,2));%Coefficient index
  
C11_f  =C1_f(xfxa+1)/(2^B_x_f_a);   C00_f=C1_f(xfxa+1).*x_f_xa+C0_f(xfxa+1);
y_f    =floor((C11_f.*x_f_xb+C00_f)*2^16)/2^16;
%Range Reconstruction
exp_ff=(exp_f+mod(exp_f,2))/2;
f     =floor(y_f.*2.^(exp_ff)*2^13)/2^13;
%%%%%%%%%%% Evaluate g0=sin(2*pi*u1)%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% g1=cos(2*pi*u1)%%%%%%%%%%%
quad  =floor(u1/2^14);
x_g_a =bitand(u1,2^14-1);
x_g_b =(2^14-1)-x_g_a;
%Approximate
%For part-a
x_g_a_xa =floor(x_g_a/2^(14-B_x_g_aa))/2^B_x_g_aa;%%
x_g_a_xb =(x_g_a-x_g_a_xa*2^14)/2^(14-B_x_g_aa);%%
xgaxa=x_g_a_xa*2^B_x_g_aa;%Coefficient index

C11_g_a=C1_g(xgaxa+1)./(2^B_x_g_aa);   
C00_g_a=C1_g(xgaxa+1).*x_g_a_xa+C0_g(xgaxa+1);
y_g_a =floor((C11_g_a.*x_g_a_xb+C00_g_a)*2^15)/2^15;
%For part-b
x_g_b_xa =floor(x_g_b/2^(14-B_x_g_ba))/2^B_x_g_ba;
x_g_b_xb =(x_g_b-x_g_b_xa*2^14)/2^(14-B_x_g_ba);%%
xgbxa=x_g_b_xa*2^B_x_g_ba;%Coefficient index

C11_g_b=C1_g(xgbxa+1)./(2^B_x_g_ba);   
C00_g_b=C1_g(xgbxa+1).*x_g_b_xa+C0_g(xgbxa+1);
y_g_b =floor((C11_g_b.*x_g_b_xb+C00_g_b)*2^15)/2^15;

g0=zeros(size(u1));
g1=g0;
for i=1:size(quad,1)
    switch quad(i)
        case 0 
            g0(i)= y_g_b(i);  g1(i)= y_g_a(i);
        case 1 
            g0(i)= y_g_a(i);  g1(i)= -y_g_b(i);
        case 2 
            g0(i)= -y_g_b(i); g1(i)= -y_g_a(i);
        case 3 
            g0(i)= -y_g_a(i); g1(i)= y_g_b(i);
    end
end
end
%AWGN output
x0=floor(f.*g0*2^11)/2^11; 
x1=floor(f.*g1*2^11)/2^11;

figure
plot(x0);
plot(x1);
[mu,s]=normfit(x0)

figure
hold on;
[Fx0,X_values] = ecdf(x0);
F_F = plot(X_values,Fx0,'r-');
set(F_F,'LineWidth',1.5);
G = plot(X_values,normcdf(X_values,0,1),'b-');
set(G,'LineWidth',1);
legend([F_F G],...
       'Empirical CDF','Standard Normal CDF',...
       'Location','SE');
   
figure
hold on
[fx0,x_values]=ksdensity(x0);
F_f =plot(x_values,fx0,'r-'); 
G_f=plot(x_values,normpdf(x_values,0,1),'b-');
set(F_f,'LineWidth',1.5);
set(G_f,'LineWidth',1);
legend([F_f G_f],...
       'Empirical PDF','Standard Normal PDF',...
       'Location','SE');