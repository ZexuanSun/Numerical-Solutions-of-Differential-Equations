function lbrb(Dt,Tend,choice)
%���ұ߽綼�ǽڵ��ʽ����߼�һ��ghost point
%11����
%hw , t = 0.1,0.9,2.0

%------------------------------------------
%step 1: set parameter: viacocity
%------------------------------------------
nu = 0.1;

%------------------------------------------
%step 2: define Spacial Domain and gride
%------------------------------------------
Xl = 0;
Xr = 1;
Ne = 10;%10������
Nn = Ne + 1;%11���߽�ڵ�
Dx = (Xr - Xl)/Ne;
i = 1:Nn;
x(i) = (i-1) * Dx + Xl;%����0, ... ,1

%------------------------------------------
%step 3: set the initial condition
%------------------------------------------
time = 0;
new(2:Nn+1) = x.*(1-x);%������ǹ��,
                       %��new(1),2-12������0-1��11���ڵ�
Crt = nu * Dt / Dx^2;
%new(1) = new(3) - 2*Dx* 10*sin(time);
%------------------------------------------
%step 5: solution
%------------------------------------------
It = 1;
Nt = Tend/Dt;
while It <= Nt
    old = new; 
    %����ghost point
    old(1) = old(3) - 2*Dx* 10*sin(time);
    %time = time + Dt;
    %------------------------------------------
    %step 4.1 : update  right boundary point
    %------------------------------------------
    %new(Nn+1) = 4*sin(6*time);
   
    %------------------------------------------
    %step 4.2 update interior point
    %------------------------------------------
    i = 2:Nn;
    new(i) =  old(i) + Dt .* sin(2*pi.*x(i-1)) * sin(4*pi*time) + Crt*(old(i-1)- 2*old(i) + old(i+1));
    %------------------------------------------
    %step 4.3 update ghost point
    %------------------------------------------
    %new(1) = new(3) - 2*Dx* 10*sin(time);
    time = time + Dt;
    new(Nn+1) = 4*sin(6*time);
    It = It +1; 
end
%11����
pnew = new(2:Nn+1);
%[xe,ue] = matlab_sol(x,Dt,Tend);
[xe,ue,error] = matlab_sol(x,pnew,Dt,Tend);
if choice == 1
     plot(x,pnew,'o--',xe,ue,'-') 
     legend('Num. Sol.','Pdepe. Sol.');
else
    
     plot(x,error,'o-')
     legend('Error')
end
sum(error)/size(x,2)
%plot(x,error,'o-')
%plot(x,ue-pnew,'o--')
 %plot(x,pnew,'o--',xe,ue,'-') 
 %legend('Num. Sol.','Pdepe. Sol.');    
end

