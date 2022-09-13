%% 
clear all
clc
load('Particle_RadiusD.mat')
Contacts=readmatrix('bp_stardist_RSAJulia_Contacts5.csv');
n_c=round(Contacts(:,2));
% Load fully populated network, with all neighbors included.
%N0=loadnetwork('Nodes_101.nod','Edges_101.edg');
N0=loadnetwork('Nodes_RSAJuliaM.nod','Edges_RSAJuliaM.edg');
P=N0.P;
for m=1:size(P,1)
    nn=nnz(N0.cons(m,:));
    Idx=knnsearch(P(N0.cons(m,1:nn)',:),P(m,:),'K',nn);
    N0.cons(m,1:nn)=N0.cons(m,Idx);
end
% Target valency for testing, which is taken as a reduction
% of the valency of the original network.
%target_valency = round(N0.valency .* rand(size(N0.valency)));
target_valency =n_c;

% Simulated annealing to find target network
nr_of_attempts = 1e7;
[N, E] = reducenetwork_mod(N0, target_valency, nr_of_attempts);
d=0;
for i1=1:size(N.cons,1)
    nn=nnz(N.cons(i1,:));
    for j1=1:nn
        par=N.cons(i1,j1);
        if par>i1
            d=d+norm(P(i1,:)-P(par,:));
        end
    end
end
%% 

figure(1);
subplot(2, 2, 1);
plot(N.valency-target_valency, '+');
xlabel('Node index');
ylabel('Residual');
subplot(2, 2, 3);
plot(E/E(1));
ylim([0 1]);
xlabel('Iteration');
ylabel('Square error');
subplot(2, 2, [2 4]);
plotnetwork(N);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

%% 
stretch=size(N.cons);
for i2=1:size(N.cons,1)
    conn=nnz(N.cons(i2,:));
    for j2=1:conn
        par=N.cons(i2,j2);
        stretch(i2,j2)=norm(P(i2,:)-P(par,:))-(Radius(par)+Radius(i2));
    end
end
disp(max(stretch,[],'all'))
%% 
index=find(N.valency==1);
par_moved=[];
for i3=size(index,1):-1:1
    par1=index(i3);
    par2=N.cons(par1,1);
    if ismember(par2,par_moved)
        continue
    end
    attempt=1;
    move=[1 0.75 0.5];
    while attempt<4
        Centres=P;
        x=Centres(par2,1)-Centres(par1,1);
        y=Centres(par2,2)-Centres(par1,2);
        z=Centres(par2,3)-Centres(par1,3);
        [azimuth,elevation,r] = cart2sph(x,y,z);
        tmp1 = azimuth < 0;
        azimuth(tmp1) = azimuth(tmp1) + 2*pi;
        elevation = pi/2-elevation ;
        Orient(1,:)=[azimuth elevation];
        R=stretch(par1,1)*move(attempt);
        if R<0
            break
        end
        Centres(par1,1)=Centres(par1,1)+R*sin(Orient(1,2))*cos(Orient(1,1));
        Centres(par1,2)=Centres(par1,2)+R*sin(Orient(1,2))*sin(Orient(1,1));
        Centres(par1,3)=Centres(par1,3)+R*cos(Orient(1,2));
        mask = (Centres(:,1)-Centres(par1,1)).^2+(Centres(:,2)-Centres(par1,2)).^2 ...
            +(Centres(:,3)-Centres(par1,3)).^2<=(Radius(par1)+Radius).^2;
        mask(par1)=0;
        if ~any(mask)
            P(par1,:)=Centres(par1,:);
            par_moved=[par_moved par1];
            break
        end
        attempt=attempt+1;
    end
end
    
%% 
stretch_n=size(N.cons);
for i2=1:size(N.cons,1)
    conn=nnz(N.cons(i2,:));
    for j2=1:conn
        par=N.cons(i2,j2);
        stretch_n(i2,j2)=norm(P(i2,:)-P(par,:))-(Radius(par)+Radius(i2));
    end
end
disp(max(stretch_n,[],'all'))    
%% 
par_pairs_d=stretch_n(index,1);
[~,index2]=sort(par_pairs_d,'descend');
distant_par=index(index2);
for i4=1:200
    par1=distant_par(i4);
    par2=N.cons(par1,1);
    if par2==0
        continue
    end
    temp=N.cons(par2,:);
    temp2=N.valency(par2)-1;
    N.cons(par2,:)=0;
    N.cons(par2,1:temp2)=setdiff(unique(temp),[par1 0]);
    N.valency(par2)=N.valency(par2)-1;
    N.valency(par1)=0;
    N.cons(par1,1)=0;
end

%% 
figure
edges=-0.5:20.5;
%histogram(n_c,edges)
histogram(props.no_contacts,edges)
hold on
histogram(N.valency,edges)
%histogram(props.no_contacts,edges)
legend('desired','recons')

