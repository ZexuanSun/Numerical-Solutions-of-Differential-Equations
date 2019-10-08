function lcrc(Dt,Tend,choice)
%左右边界都是中心格式
%10个点

%------------------------------------------
%step 1: set parameter: viacocity
%------------------------------------------
nu = 0.1;

%------------------------------------------
%step 2: define Spacial Domain and gride
%------------------------------------------
Xl = 0;
Xr = 1;
Ne = 10;%10个网格，
Nn = Ne + 1;%11个边界节点
Dx = (Xr - Xl)/Ne;
i = 1:Ne;
x = Dx/2 + Dx*(i-1);


%------------------------------------------
%step 3: set the initial condition
%------------------------------------------
time = 0;
new(2:Nn) = x.*(1-x);%最左边是 Newman边界的ghost point
%new(1) = new(2) - Dx * 10*sin(time);
Crt = nu * Dt / Dx^2;

%------------------------------------------
%step 5: solution
%------------------------------------------
It =1 ;
Nt = Tend/Dt;
while It <= Nt
    old = new;
    old(1) = old(2) - Dx * 10*sin(time);
    
    i = 2: Nn-1;
    new(i) = old(i) + Dt .* sin(2*pi.*x(i-1)) * sin(4*pi*time) + Crt*(old(i-1)- 2*old(i) + old(i+1));
    
    %new(1) = new(2) - Dx * 10*sin(time);
    new(Nn) = old(Nn) + Dt * sin(2*pi.*x(Ne)) * sin(4*pi*time) + Crt * (old(Nn-1) - ...
        3*old(Nn) + 2 * 4*sin(6*time));
     
    time = time + Dt;
    It = It + 1;
end
%10个点
pnew = new(2:Nn);
[xe,ue,error] = matlab_sol(x,pnew,Dt,Tend);
if choice == 1
     
    plot(x,pnew,'o--',xe,ue,'-') 
     legend('Num. Sol.','Pdepe. Sol.');
else
     plot(x,error,'o-')
     legend('Error')
end
sum(error)/size(x,2)
% plot(x,pnew,'o--',xe,ue,'-') 
% legend('Num. Sol.','Pdepe. Sol.');
end
