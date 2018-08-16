# Particle-swarm-optimization-algorithm
For study of carrier dynamics
ImZ=[];
f=ImZ(:,1);
RI=[ImZ(:,2) -ImZ(:,3)];
n=size(RI,1);
c1=2;          
c2=2.4;         
wmax=0.9;
w0=0.2; 
M=20;       
D=3;                   
N=60;                  
xl=[0.1,0.1,1];     
xu=[3,3,3000];  
vu=(xu-xl)*0.1;   
vl=-vu;        

%------------------

for j=1:10
    x=rand(N,D).*(ones(N,1)*(xu-xl))+ones(N,1)*xl;
    v=(rand(N,D)-0.5*ones(N,D)).*(ones(N,1)*vu);   
for i=1:N
    x(i,x(i,:)<xl)=xl(x(i,:)<xl);
    x(i,x(i,:)>xu)=xu(x(i,:)>xu);
    
    Y=dnY1(x(i,:),f);
    p(i)=sqrt(sum(sum((Y-RI).^2))/n);
    y(i,:)=x(i,:);
end
pg=x(N,:);
pbest=p(N);
for i=1:(N-1)
       if p(i)<pbest; 
          pg=x(i,:);
          pbest=p(i);
       end
end

%------------------

for t=1:M
    for k=1:N
         Y=dnY1(x(i,:),f);
         fv(k)=sqrt(sum(sum((Y-RI).^2))/n);
    end
    fvag=sum(fv)/N;
    fmin=min(fv);
    w=wmax-w0*t/M;
    for i=1:N   
        if fv(i)<=fvag
            w=w0+(fv(i)-fmin)*(wmax-w0)/(fvag-fmin);
        else
            w=wmax;
        end
        v(i,:)=w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        x(i,:)=x(i,:)+v(i,:);
        x(i,x(i,:)<xl)=xl(x(i,:)<xl);
        x(i,x(i,:)>xu)=xu(x(i,:)>xu);      
        %----------------------------------------
        Y=dnY1(x(i,:),f);
        pt=sqrt(sum(sum((Y-RI).^2))/n);
        
        if pt<p(i)
            p(i)=pt;
            y(i,:)=x(i,:);
        end      
        %-----------------------------------------
        if p(i)<pbest
            pg=x(i,:);
            pbest=p(i);
        end
        pgrec(t,:)=pg;
        pfunrec(t)=pbest;
    end
end
Y0=dnY1(pg,f);
Solution(j,:)=pg;
pt1(j)=sqrt((sum(sum((Y0-RI).^2))/n));
end
[nn,mm]=min(pt1);
pg1=Solution(mm,:);
Y1=dnY1(pg1,f);
psobest=sqrt(sum(sum((Y1-RI).^2))/n);

figure(1)
subplot(2,1,1)
semilogx(f,RI(:,1),'g.','markersize',10)
hold on
semilogx(f,Y1(:,1),'r','linewidth',2)
xlabel('Frequence(Hz)');
ylabel('ReZ(Ω)');
legend('experimental data','PSO');
subplot(2,1,2)
semilogx(f,RI(:,2),'g.','markersize',10)
hold on
semilogx(f,Y1(:,2),'r','linewidth',2)
xlabel('Frequence(Hz)');
ylabel('-ImZ(Ω)');
legend('experimental data','PSO');
