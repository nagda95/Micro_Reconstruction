%% Microstructure Characterization
tic
S=load('Segmented_Volume.mat'); 
S=S.S;
% Volumetric matrix obtained from segmentation of tomography data for particles microstructure 
no_particles=max(S,[],'all');
Properties=regionprops3(S,'Centroid','EquivDiameter','Orientation',...
    'PrincipalAxisLength','Volume','SurfaceArea');
X_coord=Properties.Centroid(:,1);Y_coord=Properties.Centroid(:,2);
Z_coord=Properties.Centroid(:,3);
Particles_radius= Properties.EquivDiameter/2;
Particles_Volume=cat(1,Properties.Volume);
SurfArea=cat(1,Properties.SurfaceArea);
Const=pi^(1/3)*6^(2/3);
Sphericity=(Const*Particles_Volume.^(2/3))./SurfArea;  
connections=cell(no_particles,1);
no_contacts=zeros(no_particles,1);
parpool
parfor i1=1:no_particles
    s1 = S==i1;
    se = ones([3 3 3]); 
    % dilate the particle with structuring element and finding it's intersection with the neighbors
    neighbors = imdilate(s1, se) & ~s1; 
    neighborLabels = unique(S(neighbors));
    neighborLabels = setdiff(neighborLabels,0);
    connections{i1}=neighborLabels;
    no_contacts(i1)=length(neighborLabels);
end
toc
fprintf('Particle Morphology Characterization started...');
%% Particle Morphology Characterization
tic
load('geode4.mat')
n_sh=4; %length of expansion limited by maximum degree L
a_lm=zeros((n_sh+1)^2,1);
SH_coeff=zeros((n_sh+1)^2,no_particles);
theta=Angles(:,2);
phi=Angles(:,1);
parfor i2=1:no_particles
    Obj=S==i2;
%     n=1;
%     centre=regionprops3(Obj,'Centroid','EquivDiameter');
%     if length(centre.EquivDiameter(:))>1
%         n=find(centre.EquivDiameter==max(centre.EquivDiameter(:)));
%     end
    [I1,I2,I3] = ind2sub(size(S),find(Obj - imerode(Obj, true(3))));
    x=I2-X_coord(i2,1);
    y=I1-Y_coord(i2,1);
    z=I3-Z_coord(i2,1);
    [azimuth,elevation,r] = cart2sph(x,y,z);
    tmp1 = azimuth < 0;
    azimuth(tmp1) = azimuth(tmp1) + 2*pi;
    elevation = pi/2-elevation ;
%     if length(unique(elevation))==1
%         continue
%     end
    F = scatteredInterpolant(elevation,azimuth,r); %radius function
    R_p=F(theta,phi); 
    a_lm=leastSquaresSHT(4,R_p,Angles,'complex');
    R_SH = real(inverseSHT(a_lm, Angles, 'complex'));
    Vol=Volume_Geode(R_SH,Structure,Angles);
    if abs(Vol-Particles_Volume(i2))<0.25*Particles_Volume(i2)
        SH_coeff(:,i2)=a_lm;
    end
end

%% Ordering SH Coefficients in terms of volume 
[Vols_order,Ind] = sort(Particles_Volume, 'ascend');
New_SH_order=SH_coeff(:,Ind);


%% Creating bins for histogram
edges=[0 1000 2000 4000 6000 8000 10000 12000 14000 16000 18000 20000 max(Vols_order(:))];

Y_bin = discretize(Vols_order,edges)';


%% Final SH Coefficients
columnsWithAllZeros = all( New_SH_order== 0);

Final_SH=New_SH_order(:,~columnsWithAllZeros);
Y_bin=Y_bin(:,~columnsWithAllZeros);

%% %% Angular power spectrum
A_ps=zeros(n_sh,size(Final_SH,2));
for o=1:size(A_ps,2)
    q=2;
    for l=1:n_sh
        m=q:q+2*l;
        A_ps(l,o)=sum(abs(Final_SH(m,o)).^2)/length(m);
        q=q+2*l+1;
    end
end

%% Average Angular power spectrum
n_bin=max(Y_bin);
A_ps_avg=zeros(size(A_ps,1),n_bin);
for b=1:n_bin
    A_ps_avg(:,b)=mean(A_ps(:,find(Y_bin==b)),2);
end
toc
fprintf('Characterization complete...');

%% relation between first SH coefficient a0 (x) and volume of the particle(y)
%finding coefficients of equation y=ax^3+bx by curve fitting to get the relation
tic
ab_coeff=zeros(2,n_bin);
for b=1:n_bin
    R=zeros(size(Angles,1),1);
    all_a0s=abs(Final_SH(1,find(Y_bin==b)));
    fn_x_values=[min(all_a0s) mean(all_a0s) max(all_a0s)]';
    Vol_par=zeros(3,1);
    Sigma_rand=Sp_Coeffs_PS(A_ps_avg(:,b),n_sh);
    mu_rand=zeros(size(Sigma_rand,1),1);
    a_rand=mvnrnd(mu_rand,Sigma_rand);
    for i3=1:length(fn_x_values)
        a0_tmp=fn_x_values(i3);
        for i4=1:size(Angles,1)
            R(i4,1)=a0_tmp*0.2821+a_rand*(Yharmonic_coeffs(n_sh,Angles(i4,1),Angles(i4,2)))';
        end
        if any(R<=0)
            a_rand=mvnrnd(mu_rand,Sigma_rand);
            for i5=1:size(Angles,1)
                R(i5,1)=a0_tmp*0.2821+a_rand*(Yharmonic_coeffs(4,Angles(i5,1),Angles(i5,2)))';
            end
        end
        Vol_par(i3)=Volume_Geode(R,Structure,Angles);
    end
    fn_y_values=Vol_par;
    ft = fittype('a*x^3+b*x');
    fn_fit=fit(fn_x_values,fn_y_values,ft);
    ab_coeff(1,b)=fn_fit.a;
    ab_coeff(2,b)=fn_fit.b;
end
toc
fprintf('Reconstruction started...');

%% Random Sequential Adsorption (RSA)
tic
Xlim=size(S,1);
Ylim=size(S,2);
Zlim=size(S,3);
Width=[Xlim, Ylim,Zlim];
Radius = sort(Particles_radius, 'descend');
n=no_particles;
P= zeros(n, 3);
R2=Radius;
iLoop  = 1;                     
index2=[];
nTrials = 4*n;
nIter = 100;
while any(R2) && iLoop <= nIter
    if iLoop==1
        Radi = Radius(1);
        P(1,:) = rand(1,3) .* (Width - 2 * Radi) + Radi;
        index2=[index2;1];
        R2(1)=0;
        packed = 1;
    end
    XX=find(R2);
    
    for kLoop=1:size(XX,1)
        Radi =  Radius(XX(kLoop));
        newP_G = rand(nTrials,3) .* (Width - 2 * Radi) + Radi;
        for jLoop=1:nTrials
            newP = newP_G(jLoop,:);
            Dist2 = sum((P(index2, :) - newP) .^ 2, 2)-(Radius(index2) + Radi).^2;

            if ~any(Dist2<0)
                P(XX(kLoop),:) = newP;
                index2=[index2;XX(kLoop)];
                R2(XX(kLoop))=0;
                packed =packed+1;
                break
            end
        end
    end
    iLoop =iLoop+1;
end
toc
fprintf('Initial Structure obtained...');

%% Plot initial random microstructure
tic
figure
axes('NextPlot', 'add', ...
    'XLim', [0, Xlim], 'YLim', [0, Ylim], 'ZLim', [0, Zlim]);
xlabel('x (pixel)')
ylabel('y (pixel)')
zlabel('z (pixel)')
view(3);
[X, Y, Z] = sphere();
for k = 1:n
    FV=surf(X * Radius(k) + P(k, 1), Y * Radius(k) + P(k, 2), Z * Radius(k) + P(k, 3));
    %         ^^^                 ^^^                 ^^^
    % Bugfix based on Jared Cole's comment. Thanks!
end

%% Data for Copula
Copula_inp(:,1)=Particles_radius;
Copula_inp(:,2)=no_contacts;
x3=zeros(no_particles,1);
%finding number of neighboring particles within 10*Radius distance
for i6=1:no_particles
   mask1 = (X_coord(:,1)-X_coord(i6,1)).^2+(Y_coord(:,1)-Y_coord(i6,1)).^2 ...
   +(Z_coord(:,1)-Z_coord(i6,1)).^2<=(10*Particles_radius(i6))^2;
   mask1(i6)=0;
   x3(i6)=nnz(mask1);
end
Copula_inp(:,3)=x3;
%jittering the data
Copula_inp(:,1)=Copula_inp(:,1)+ 0.01.*(rand(length(x3(:,1)),1));
Copula_inp(:,2)=Copula_inp(:,2)+ 0.01.*(rand(length(x3(:,1)),1));
Copula_inp(:,3)=Copula_inp(:,3)+ 0.01.*(rand(length(x3(:,1)),1));
writematrix(Copula_inp,'Copula_rawData_inp.csv')

x3=zeros(no_particles,1);
for i6=1:no_particles
   mask2 = (P(:,1)-P(i6,1)).^2+(P(:,2)-P(i6,2)).^2+(P(:,3)-P(i6,3)).^2<=(10*Radius(i6))^2;
   mask2(i6)=0;
   x3(i6)=nnz(mask2);
end
Copula_pred(:,1)=Radius;
Copula_pred(:,2)=x3;
Copula_pred(:,1)=Copula_pred(:,1)+ 0.01.*(rand(length(x3(:,1)),1));
Copula_pred(:,2)=Copula_pred(:,2)+ 0.01.*(rand(length(x3(:,1)),1));
writematrix(Copula_pred,'Copula_Contacts_pred.csv')

%% Fit bounded Voronoi to Random seed points from RSA to get neighbors
L1=Xlim; % length of the box in the x-direction
L2=Ylim; % length of the box in the y-direction
L3=Zlim; % length of the box in the z-direction
c=[round(L1/2);round(L2/2);round(L3/2)]; % centre of the box
[V_L,C_L] = croppedLaguerre3D(P(:,1),P(:,2),P(:,3),Radius,L1,L2,L3,c);
[TargVolumes,TargVoisins]=PropVoronoi3D(P(:,1),P(:,2),P(:,3),V_L,C_L);
Neighbors=TargVoisins;


%% Run R-Copula Script
% Rpath ='C:\Program Files\R\R-4.1.1\bin';
% RscriptFileName='Copula_Contacts_Pred.R';
% RunRcode(RscriptFileName, Rpath);
exec=system('Rscript.exe Copula_Contacts_Prediction.R');
toc
fprintf('#Contacts predicted from Copula...');

%% Nod;Edg file for mk_network
%contact network/connection graph
tic
connections=Neighbors;
s=[];
t=[];
for i8=1:size(connections,1)
    if isempty(connections{i8})
        continue
    end
    for j=1:size(connections{i8},2)
        if connections{i8}(j)>i8
            s=[s i8];
            t=[t connections{i8}(j)];
        end
    end
end
G = graph(s,t);
deg = degree(G);
%plot(G,'XData',x,'YData',y,'ZData',z)
fid = fopen('Nodes_101.txt','w');
for j1=1:size(deg,1)
    fprintf(fid,'%-4d %-8.3f %-8.3f %-8.3f %-2d\n',j1,P(j1,1),P(j1,2),P(j1,3),deg(j1));
end
fclose(fid);
fid2 = fopen('Edges_101.txt','w');
for j2=1:size(s,2)
    fprintf(fid2,'%-4d %-4d %-4d\n',j2,s(j2),t(j2));
end
fclose(fid2);

%% Rename files
movefile 'Edges_101.txt' 'Edges_101.edg'
movefile 'Nodes_101.txt' 'Nodes_101.nod'
toc
fprintf('Optimization for connection graph started...');

%% MK_NETWORK
tic
Contacts=readmatrix('Copula_Contacts_out.csv');
n_c=round(Contacts(:,2));
% Load fully populated network, with all neighbors included.
N0=loadnetwork('Nodes_101.nod','Edges_101.edg');
Ps=N0.P;
for i9=1:size(Ps,1)
    nn=nnz(N0.cons(i9,:));
    Idx=knnsearch(Ps(N0.cons(i9,1:nn)',:),Ps(i9,:),'K',nn);
    N0.cons(i9,1:nn)=N0.cons(i9,Idx);
end

target_valency =n_c;

% Simulated annealing to find target network
nr_of_attempts = 1e7;
[N, E] = reducenetwork_mod(N0, target_valency, nr_of_attempts);
% d=0;
% for i1=1:size(N.cons,1)
%     nn=nnz(N.cons(i1,:));
%     for j1=1:nn
%         par=N.cons(i1,j1);
%         if par>i1
%             d=d+norm(Ps(i1,:)-Ps(par,:));
%         end
%     end
% end
%% Bring particles with single contact closer to each other to avoid stretching
index=find(N.valency==1);
particle_moved=[];
for i11=size(index,1):-1:1
    par1=index(i11);
    par2=N.cons(par1,1);
    stretch=norm(P(par1,:)-P(par2,:))-(Radius(par1)+Radius(par2));
    if ismember(par2,particle_moved)
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
        Orient1(1,:)=[azimuth elevation];
        R_gap=stretch*move(attempt);
        if R_gap<0
            break
        end
        Centres(par1,1)=Centres(par1,1)+R_gap*sin(Orient1(1,2))*cos(Orient1(1,1));
        Centres(par1,2)=Centres(par1,2)+R_gap*sin(Orient1(1,2))*sin(Orient1(1,1));
        Centres(par1,3)=Centres(par1,3)+R_gap*cos(Orient1(1,2));
        mask = (Centres(:,1)-Centres(par1,1)).^2+(Centres(:,2)-Centres(par1,2)).^2 ...
            +(Centres(:,3)-Centres(par1,3)).^2<=(Radius(par1)+Radius).^2;
        mask(par1)=0;
        if ~any(mask)
            P(par1,:)=Centres(par1,:);
            particle_moved=[particle_moved par1];
            break
        end
        attempt=attempt+1;
    end
end
toc
fprintf('Contact network obtained...');
    
%% Optional manipulation
stretch_n=size(N.cons);
for k3=1:size(N.cons,1)
    conn=nnz(N.cons(k3,:));
    for k4=1:conn
        par=N.cons(k3,k4);
        stretch_n(k3,k4)=norm(P(k3,:)-P(par,:))-(Radius(par)+Radius(k3));
    end
end
%disp(max(stretch_n,[],'all'))    

par_pairs_d=stretch_n(index,1);
[~,index11]=sort(par_pairs_d,'descend');
distant_par=index(index11);
n_break=200; %break some contacts to align the histogram with experimental data
for k5=1:n_break
    par1=distant_par(k5);
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
connections2=N.cons;
tic
fprintf('Realistic microstructure reconstruction started...');

%% Microstructure Reconstruction
FileName='stardist_tot169_PubliFacet.mat';                     %Name of the created MATLAB data file for the sample
MeshName='stardist_tot169_PubliFacetMesh.mat';                            %Name of the created MATLAB data file for the constrained mesh

Domain='Domain.mat';  
load(Domain)


[Vtot]=Volume_Domain(Points_Domain,Struc_Domain);
load('geode3.mat')
Structure3=Structure;
Angles3=Angles;
geode='geode4.mat';
load(geode)
Grains=cell(1,1);
Count_s=0;
Radius_par=Radius;
VolumesGrains=zeros(size(Radius_par,1),1);

S_harmonics=zeros(size(Angles,1),(n_sh+1)^2-1);
for k1=1:size(Angles,1)
    S_harmonics(k1,:)=Yharmonic_coeffs(n_sh,Angles(k1,1),Angles(k1,2));
end

for k2=1:size(Radius_par,1)
    disp(['Filling Particle ',num2str(k2),'/',num2str(size(Radius_par,1))])
    Vol_particle=(4*pi*(Radius_par(k2))^3)/3;
    bin_no = discretize(Vol_particle,edges)';
    coeffs=[ab_coeff(1,bin_no) 0 ab_coeff(2,bin_no) -1*Vol_particle]; %Coefficients of equation relation volume and a0
    a0_roots = roots(coeffs);
    a0 = a0_roots(find(a0_roots==real(a0_roots)));
    a0final=a0(a0>0);
    avg_rad(k2)=a0final;
    Sigma=Sp_Coeffs_PS(A_ps_avg(:,bin_no),n_sh);
    mu=zeros(size(Sigma,1),1);
    conns=nnz(connections2(k2,:));
    Centre=N.P(k2,:);
    if conns>0
        Y_sh=zeros(conns,(n_sh+1)^2-1);
        Orient=zeros(size(Y_sh,1),2);
        first_coeffs=zeros(size(Y_sh,1),1);
        Radius_con=zeros(conns,1);
        for j3=1:conns
            par=connections2(k2,j3);
            x=P(par,1)-P(k2,1);
            y=P(par,2)-P(k2,2);
            z=P(par,3)-P(k2,3);
            [azimuth,elevation,r] = cart2sph(x,y,z);
            tmp1 = azimuth < 0;
            azimuth(tmp1) = azimuth(tmp1) + 2*pi;
            elevation = pi/2-elevation ;
            Orient(j3,:)=[azimuth elevation];
            [Y_sh(j3,:),first_coeff]=Yharmonic_coeffs(n_sh,azimuth,elevation);
            first_coeffs(j3)=first_coeff;
            d=norm(P(k2,:)-P(par,:));
            R1=Radius_par(k2);
            R2=Radius_par(par);
            stretch=(d-(R1+R2))*(R1/(R1+R2));
            Radius_con(j3)=R1+stretch-a0final*0.2821;
            R_stretch(k2,j3)=stretch;
        end
        [U,S,V] = svd(Y_sh);
        mu=transpose(V)*mu;
        Sigma=V*Sigma*transpose(V);
        conditional_vals=zeros(size(Radius_con,1),1);
        U=transpose(U);
        for j4=1:size(Radius_con,1)
            conditional_vals(j4)=U(j4,:)*Radius_con/S(j4,j4);
        end
        [mubar,Sigmabar]=mvncondP(mu,Sigma,conditional_vals);
        %drawing SH coefficients from a multivariate gaussian distribution with conditional values.

        Cond_rand_vector=mvnrnd(mubar,Sigmabar);
        a=V*[conditional_vals;Cond_rand_vector'];
        Rgrain(:,1)=a0final*0.2821+S_harmonics*a;
    else
        Rgrain=zeros(size(S_harmonics,1),1);
        iter=1;
        while any(Rgrain<=0) && iter<6
            a=mvnrnd(mu,Sigma)';
            Rgrain(:,1)=a0final*0.2821+S_harmonics*a;
            iter=iter+1;
            if iter==6
                Count_s=Count_s+1;
                disp('Ooops')
            end
        end
    end
    %Rgrain(:,1)=a0final*0.2821+S_harmonics*a;
    Xgrain=Rgrain.*sin(Angles(:,2)).*cos(Angles(:,1))+Centre(1);
    Ygrain=Rgrain.*sin(Angles(:,2)).*sin(Angles(:,1))+Centre(2);
    Zgrain=Rgrain.*cos(Angles(:,2))+Centre(3);
    Grains{k2,1}=cat(2,Xgrain,Ygrain,Zgrain,Rgrain);
    VolumeGrain=Volume_Geode(Rgrain,Structure,Angles);
    VolumesGrains(k2,1)=VolumeGrain; 
end
save(MeshName)
disp(' ')
save(FileName)
Plot_Particles('Geode4.mat',Domain,Centres,Grains,VolumesGrains,FileName)
toc

