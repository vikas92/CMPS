% Program for Gauss - Seidel Load Flow Analysis
% Vikas Gupta
%https://github.com/vikas92
%Change the acceleration factor and number of iterations as desired.
%Edit busdata.m and lineimp.m to enter the values corresponding to the
%system under study

%All values to be entered in pu only.

clc;
alpha=1.4;          %Acceleration Factor
iterations=10;      %Number of iterations
data=lineimp();
fb = data(:,1);     % From bus number...
tb = data(:,2);     % To bus number...
LineCharging=data(:,4);
nbus=max([max(fb);max(tb)]);
nbranch=length(fb);
admittanceToGround=[];
for ii=1:nbus
    index=[find(data(:,1)==ii);find(data(:,2)==ii)];
    temp=0;
    for jj=1:size(index)
        temp=temp+data(index(jj),4);
    end
    admittanceToGround=[admittanceToGround;temp];
    index=[];
end
disp('Bus Code                     admittanceToGround')
disp([(1:nbus)' admittanceToGround])
b=admittanceToGround(:,1);
LineAdmittance=1./data(:,3);
disp('Line admittances')
disp('FromBus                ToBus               LineAdmittance')
y=LineAdmittance;
disp([fb tb y])
ybus = zeros(nbus,nbus);        % Initialise YBus...

% Formation of the Off Diagonal Elements...
for k=1:nbranch
    ybus(fb(k),tb(k)) = -y(k);
    ybus(tb(k),fb(k)) = ybus(fb(k),tb(k));
end

% Formation of Diagonal Elements....
for m=1:nbus
    for n=1:nbranch
        if fb(n) == m | tb(n) == m
            ybus(m,m) = ybus(m,m) + y(n);
        end
    end
    ybus(m,m)=ybus(m,m)+b(m);
end
ybus                  % Bus Admittance Matrix
zbus = inv(ybus);      % Bus Imp
busdata=busdata();
GenMW = busdata(:,3);       % PGi, Real Power injected into the buses.
GenMVAR = busdata(:,4);     % QGi, Reactive Power injected into the buses.
LoadMW = busdata(:,5);      % PLi, Real Power Drawn from the buses.
LoadMVAR = busdata(:,6);    % QLi, Reactive Power Drawn from the buses.
P = GenMW - LoadMW;         % Pi = PGi - PLi, Real Power at the buses.
Q = GenMVAR - LoadMVAR;     % Qi = QGi - QLi, Reactive Power at the buses.
KLp=[];
for m=1:nbus
    temp=(P(m)-Q(m)*i)/ybus(m,m);
    KLp=[KLp;temp];
end
disp('Bus Code                     KLp')
disp([(1:nbus)' KLp])
YLpq=zeros(nbus,nbus);
for m=1:nbus
    for n=1:nbus
        YLpq(m,n)=ybus(m,n)/ybus(m,m);
    end
end
disp('Diagonal elements in the following matrix carry no meaning.Consider only off-diagonal ones')
YLpq
E=busdata(:,2);
delE=zeros(nbus,1);
itertable=[];
delitertable=[];
for iter=1:iterations
    iter;
    for m=2:nbus
        temp=E(m);
        E(m)=(KLp(m)/conj(E(m)))-YLpq(m,:)*E+YLpq(m,m)*E(m);
        delta=E(m)-temp;
        E(m)=temp+alpha*delta;
        delE(m)=alpha*delta;
    end
    itertable=[itertable;transpose(E)];
    delitertable=[delitertable;transpose(delE)];
end
disp('Bus Voltages from the Gauss-Seidel iterative solution')
disp(itertable)
disp('Changes in bus Voltages from the Gauss-Seidel iterative solution')
disp(delitertable)
LineFlow=zeros(nbus,nbus);
for n=1:nbranch
    fi=fb(n);
    ti=tb(n);
    LineFlow(fi,ti)=conj(E(fi))*(E(fi)-E(ti))*(-ybus(fi,ti))+conj(E(fi))*E(ti)*LineCharging(n);
    fi=tb(n);
    ti=fb(n);
    LineFlow(fi,ti)=conj(E(fi))*(E(fi)-E(ti))*(-ybus(fi,ti))+conj(E(fi))*E(ti)*LineCharging(n);
end
LineFlow=conj(LineFlow);
disp('Calculated line flows in pu are:')
disp(LineFlow)
