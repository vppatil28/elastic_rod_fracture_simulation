%%%%%
% Code simulating fracture in extensible rod
% Input: end-to-end distance, in cm, and twist Tw, in degrees, of rod
% This tool assumes rod has length 24cm and is on the point of fracture 
% Simulation show future fracture events.


% Plotting tool obtained from Matlab file exchange, used under the following license:
% 
% Copyright (c) 2016, Janus H. Wesenberg
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


input = csvread('input_data.txt', 0, 0, [0,0,0,7]);

h=input(1); %radius of the rod

rho=input(2); %density in g/cm^3
E=input(3); %young's modulus
nu=input(4); %poisson's ratio

omega0=input(5); %twist damping
lambda_f = input(6); %min fragment size
endDist = input(7); %end-to-end distance in cm, should be between 13.5 and 23.5
Twist=input(8)*pi/180; %total twist in rod

mutwist=E/(2*(1+nu)); %shear modulus
A=0.5*pi*h^2; %cross section area
I=0.25*pi*h^4; %moments of inertia
J=0.5*pi*h^4; %twisting energy

beta=mutwist*J;
alpha=E*I;
gamma=E*A;

load('initial_rod_shapes.mat','X');
X0=X(:,:,round(1176-50*endDist));


n=size(X0);
n=n(2)-2;
L0=sqrt(sum((X0(:,2:end)-X0(:,1:end-1)).^2));
Ls0=Vdom(X0);


%quasistatically relax the rod so it has the right shape given the twist
%angle
%%%%%%%%%%%%%%
BC='n';
N=10000; %number of timesteps
dt=10^(-8);
Ftot=bForce1(X0,alpha,BC)+0.1*sForce(X0,L0,gamma,BC)+tForce(X0,qtAngle(X0,Twist),beta,BC);
X(:,:,1)=X0+(dt)*Ftot; %update velocities
X(3,2,1)=0; X(3,n+1,1)=0; %clamping
dX=project(X(:,:,1),L0);
dX=reshape(dX,[3,n+2]);
X(:,2:n+1,1)=X(:,2:n+1,1)+dX(:,2:n+1);

for i=2:N
    Ftot=bForce1(X(:,:,i-1),alpha,BC)+0.1*sForce(X(:,:,i-1),L0,gamma,BC)+tForce(X(:,:,i-1),qtAngle(X(:,:,i-1),Twist),beta,BC);
    X(:,:,i)=X(:,:,i-1)+(dt)*Ftot; %update velocities
    X(3,2,i)=0; X(3,n+1,i)=0; %clamping
    dX=project(X(:,:,i),L0);
    dX=reshape(dX,[3,n+2]);
    X(:,2:n+1,i)=X(:,2:n+1,i)+dX(:,2:n+1);
    for ii=1:10
        Ftot=bForce1(X(:,:,i),alpha,BC)+(1/3)*sForce(X(:,:,i),L0,gamma,BC)+tForce(X(:,:,i),qtAngle(X(:,:,i),Twist),beta,BC);
        X(:,:,i)=X(:,:,i)+(dt)*Ftot; %update velocities
        X(3,2,i)=0; X(3,n+1,i)=0; %clamping
        dX=project(X(:,:,i),L0);
        dX=reshape(dX,[3,n+2]);
        X(:,2:n+1,i)=X(:,2:n+1,i)+dX(:,2:n+1);
    end
    if mod(i,1000)==0
        fprintf('\nTwist relaxtaion %d%% done\n', i/100)
    end
end
X(3,2,N)=0; X(3,n+1,N)=0; %clamping
%%%%%%%%%%%%%%


X0 = X(:,:,end);

N=2000; %number of timesteps
framenumber=100;
s=round(N/framenumber); %N/s is number of frames to capture
dt=10^(-6); %size of timesteps

X=zeros(3,n+2,N); X(:,:,1)=X0;
V0=zeros(3,n+2);
V=zeros(3,n+2,N);

BC='n'; %no free ends initially

Ftot1=bForce1(X0,alpha,BC)+sForce(X0,L0,gamma,BC)+tForce(X0,qtAngle(X0,Twist),beta,BC);
V(:,:,1)=V0; %verlet step
X(:,:,2)=X0+(dt)*V0+0.5*dt^2*Ftot1./(Ls0*rho*A);
X(:,2,2)=X0(:,2); X(:,end-1,2)=X0(:,end-1); %fix first/last two points
dX=project(X(:,:,2),L0);
dX=reshape(dX,[3,n+2]);
X(:,3:n,2)=X(:,3:n,2)+dX(:,3:n);


%%%%initialize frame
d1=zeros(3,n+1,N); %initial frame will be bishop frame
d3=X0(:,2:end)-X0(:,1:end-1);
No=cross(d3(:,1:end-1),d3(:,2:end));
No=No./sqrt(dot(No,No));
d3=d3./sqrt(dot(d3,d3));
ph=phi(X0);
ph=ph(2:end-1);
d1(:,1,1)=[0,0,1]';
for i=2:n+1
    d1(:,i,1)=cos(ph(i-1)).*d1(:,i-1,1)+((1-cos(ph(i-1))).*dot(No(:,i-1),d1(:,i-1,1))).*No(:,i-1) + (sin(ph(i-1)).*cross(No(:,i-1),d1(:,i-1,1)));
end

theta=zeros(1,n+1,N);
m0=qtAngle(X0,Twist);
theta(:,:,1)=[0,cumsum(m0)]; %initial angle distribution is 

%twist Verlet step 1
Vt=zeros(1,n+1); %initially at rest
theta(:,:,2)=theta(:,:,1)+dt*Vt+0.5*dt^2*angleForce(X0,m0,beta,'n')./(L0*rho*pi*h^4*0.5);
d1(:,:,2)=fupdate(X0,X(:,:,2),d1(:,:,1));
mAngle=tAngle(X(:,:,2),d1(:,:,2),theta(:,:,2));
%%%%

critStress=max(stress(X(:,:,2), mAngle, E,I,mutwist,J)); %set critical stress
TwD=zeros(1,N); %max twist density at i'th timestep
fullK=zeros(N,n+2); %full curvature stress
fullTwD=zeros(N,n+2); %full twist stress
fullStress=zeros(N,n+2); %full stress profile
K=zeros(1,N); %max curvature at the i'th timestep
K(1)=max(kappa(X0));
fullmAngle=zeros(N,n); fullmAngle(1,:)=m0;
ctrlPoints=zeros(N,n+2);
Ls1=Vdom(X(:,:,2)); %initial voronoi domains
b=zeros(1,N)+n+2; %number of active points, decreases for each breaking event
for i=2:N      
    K(i)=max(kappa(X(:,1:b(i),i))); %max curvature
    %check for fracture
    if max(stress(X(:,1:b(i),i), mAngle, E,I,mutwist,J)) >= critStress
        [m,ii]=max(stress(X(:,1:b(i),i), mAngle, E,I,mutwist,J)); %update b
        %fractures can't be too close
        if sum(Ls1(ii+1:b(i)-1))>lambda_f && sum(Ls1(1:ii-1))>lambda_f
            b(i:N)=ii; %set new cutoffs
            BC='r'; %after the break we have a free end
        end
    end
    %New quantities to calculate forces
    Ls1=Vdom(X(:,1:b(i),i)); %new voronoi domains
    mAngle=tAngle(X(:,1:b(i),i),d1(:,1:b(i)-1,i),theta(:,1:b(i)-1,i)); %new twist density (post fracture)
   
    %store values for curvature, twist density
    fullK(i,1:b(i))=(kappa(X(:,1:b(i),i))); %store the curvature profile
    fullTwD(i,2:b(i)-1)=mutwist^2*J*(mAngle./Ls1(2:b(i)-1)).^2; %twist density profile
    fullStress(i,1:b(i))=stress(X(:,1:b(i),i), mAngle, E,I,mutwist,J);
    fullmAngle(i,1:b(i)-2) = mAngle;
    ctrlPoints(i,1:b(i))=cumsum(Vdom(X(:,1:b(i),i))); %ctrlPoints(i,1:b(i))=ctrlPoints(i,1:b(i))-ctrlPoints(i,1);

    %forces on material points
    Ftot=bForce1(X(:,1:b(i),i),alpha,BC)+0.1*sForce(X(:,1:b(i),i),L0,gamma,BC)+tForce(X(:,1:b(i),i),mAngle,beta,BC);
    X(:,1:b(i),i+1)=2*X(:,1:b(i),i)-X(:,1:b(i),i-1)+(dt^2)*Ftot./(Ls1(1:b(i))*rho*A);
    X(:,2,i+1)=X(:,2,i); %fix first 2 points
    %%%manifold projection
    dX=project(X(:,1:b(i),i+1),L0(1:b(i)-1));
    dX=reshape(dX,[3,b(i)]);
    X(:,3:b(i),i+1)=X(:,3:b(i),i+1)+dX(:,3:b(i));
    %%%
    X(:,b(i)+1:end,i+1)=X(:,b(i)+1:end,i); %update second half by symmetry
    
    %forces on material frame
    

    theta(:,1:b(i)-1,i+1)=theta(:,1:b(i)-1,i) + (1+omega0*dt)^(-1) * ((1-omega0*dt) * (theta(:,1:b(i)-1,i)-theta(:,1:b(i)-1,i-1))+dt^2*angleForce(X(:,1:b(i),i),mAngle,beta,BC)./(L0(1:b(i)-1)*rho*pi*h^4*0.5));  
    d1(:,1:b(i)-1,i+1)=fupdate(X(:,1:b(i),i),X(:,1:b(i),i+1),d1(:,1:b(i)-1,i));
    mAngle=tAngle(X(:,1:b(i),i+1),d1(:,1:b(i)-1,i+1),theta(:,1:b(i)-1,i+1)); %new twist density (pre fracture)
    
    %%%
    if mod(i, 200)==0
        fprintf('\nIteration %d%%\n', i/20)
    end
end

save('rod_fracture.mat','X');


%redefine colourscheme


%set colour scheme, axis
smax=2*10^7; vc=[0,smax];
shape=[X0(1,(n+3)/2)-11.25, X0(1,(n+3)/2)+11.25, -0.5, 10, -5, 5];

 

%%%%%create AVI object
nFrames = N;
vidObj = VideoWriter('rod_fracture.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);


%%%%create movie
for j=[zeros(1,11)+1,s+1:s:nFrames] %initial frame is repeated 10 times so initial condition is visible in final video
    c=[0,unique(b(1:j))]; %c contains the fracture points
    S=[10^(0)*stress(X(:,1:c(2),j), fullmAngle(j,1:c(2)-2), E, I, mutwist, J), zeros(1,n+2-c(2))-1];
    S(1)=S(2); S(end)=S(end-1);
    pieceNo = [1:max(length(c)-2,1), 1:length(c)-2];
    for mm = 1:length(pieceNo)
        jj=pieceNo(mm);
        
        %increase gaps between pieces to aid visualization
        X(:,26,j)=0.2*X(:,25,j) + 0.8*X(:,26,j);
        if mm<=max(length(c)-2,1)
            pX=X(:,max(c(jj)+1,1):c(jj+1),j);
        else
            v=0.5*(X(:,end,j)-X(:,1,j));
            pX=[-1,0,0;0,1,0;0,0,-1]*(X(:,max(c(jj)+1,1):c(jj+1),j)-v)+v;
        end
        if j==1
            pX=X0;
            S=10^(0)*stress(X0, qtAngle(X0,Twist), E, I, mutwist, J);
            S(1)=S(2); S(end)=S(end-1);
        end
        
        
        %third party plotting tool, see disclaimer in function description    
        [x,y,z,C]=tubeplot2(pX,0.2,8,h/2,0.1*sqrt(1/(pi*h^2))*S(max(c(jj)+1,1):c(jj+1)));
        
        A=surf(x,y,z,C);
        set(A, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
        shading interp;
        set(A,'meshstyle','row');
        axis(shape); view(2); caxis(vc); 
        cb2=colorbar;
        cb2.Position=[0.8550, 0.275, 0.02, 0.485];
        cb2.Ticks=[0,1,2]*10^7;
        cb2.Label.String='Eff. Stress (N/m)';
        set(cb2, 'FontSize',16,'FontName','Helvetica');
        
        
        set(gca, 'XTick',[]); set(gca, 'YTick',[]);
        daspect([1,1,1]);
        grid off; box on;
        material shiny;
        set(gca,'Color',[0.6 0.6 0.6]);
        set(gca,'position',[0.060 0.1100 0.7750 0.8150]);
        drawnow;
        hold on;
    end
    camlight left; camlight right;
    fill([15.3,18.3,18.3,15.3],[9.5,9.5,9.1,9.1],[1,1,1],'EdgeColor','none'); %scale bar 3cm
    writeVideo(vidObj, getframe(gcf));
    hold off;       
    fprintf('\nPlotting frame %d\n', j)
end
close(vidObj);
hold off;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%elastic rod functions


% function L = rodLength(X) %output length of rod described by 3 x n matrix
%     n=size(X);
%     n=n(2)-2;
%     A=X(:,2:n+2)-X(:,1:n+1);
%     L=trace((A'*A).^0.5);
% end


function kb = darboux(X) %calculate darboux frame at each point along the rod
    d3=(X(:,2:end)-X(:,1:end-1));
    t=d3./sqrt(dot(d3,d3));
    kb=2*cross(t(:,1:end-1),t(:,2:end))./(1+dot(t(:,1:end-1),t(:,2:end)));
end




function Fs = sForce(X0,L0,gamma,BC) %stretching forces at each point of deformed rod
    n=size(X0);
    n=n(2)-2;
    Ys=zeros(3,n+2);
    if BC=='l' || BC=='b'
        t0=X0(:,2)-X0(:,1);
        Ys(:,1)=(gamma/(norm(t0)*L0(1)))*(norm(t0)-L0(1))*t0;
    end
    if BC=='r'|| BC=='b'
        tf=X0(:,n+2)-X0(:,n+1);
        Ys(:,n+2)=-(gamma/(norm(tf)*L0(n+1)))*(norm(tf)-L0(n+1))*tf;
    end
    d3=X0(:,2:n+2)-X0(:,1:n+1);
    Ln=sqrt(sum(d3.^2));
    Fu=(gamma./(Ln(2:n+1).*L0(2:n+1))).*(Ln(2:n+1)-L0(2:n+1)).*d3(:,2:n+1);
    Fd=-(gamma./(Ln(1:n).*L0(1:n))).*(Ln(1:n)-L0(1:n)).*d3(:,1:n);    
    Ys(:,2:n+1)=Fu+Fd;
    Fs=Ys;
end





function ph = phi(X) %turning angle from one rod segment to the next 
    n=size(X,2)-2;
    d3=X(:,2:n+2)-X(:,1:n+1);
    Y=atan2(sqrt(sum(cross(d3(:,1:n),d3(:,2:n+1)).^2,1)),dot(d3(:,1:n),d3(:,2:n+1)));
    ph=[0,Y,0];
end



function Ls = Vdom(X) %vector containing size of voronoi domain at each X
    n=size(X,2)-2;
    Y=(X-[X(:,1),X(:,1:n+1)]);
    Z=([X(:,2:n+2),X(:,n+2)]-X);
    Ls=0.5*(sqrt(sum(Y.^2))+sqrt(sum(Z.^2)));
end


function k = kappa(X) %curvature at each point of rod
    ph=phi(X);
    k=ph./Vdom(X);
end

      
function dp = dphi(X) %dp(i,j,:) = dphi(i)/dx(j)
    n=size(X);
    n=n(2)-2;
    Yp=zeros(n+2,n+2,3);
    d3=X(:,2:n+2)-X(:,1:n+1);
    cosph = (dot(d3(:,1:n), d3(:,2:n+1))./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))));
    dcos = -1./sqrt(1-cosph.^2);
    ld = -d3(:,2:n+1)./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))) + (d3(:,1:n)./(sum(d3(:,1:n).^2))).*cosph;
    ud = d3(:,1:n)./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))) - (d3(:,2:n+1)./(sum(d3(:,2:n+1).^2))).*cosph;
    d=[[0,0,0]',-dcos.*(ld+ud),[0,0,0]'];
    ld=[dcos.*ld,[0,0,0]']; ud=[[0,0,0]',dcos.*ud];
    Yp(:,:,1)=(diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1));
    Yp(:,:,2)=(diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1));
    Yp(:,:,3)=(diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1));
    Yp=real(Yp);
    Yp(isnan(Yp))=0; Yp(isinf(Yp))=0;
    dp=Yp;
end


function derivs = dL(X) %derivs(i,j,:) = dL(i)/dx(j) where L(i) is twice the voronoi domain of i'th vertex
    n=size(X);
    n=n(2)-2;
    Yp=zeros(n+2,n+2,3);
    ld=(X(:,1:n+1)-X(:,2:n+2))./(sqrt(sum((X(:,1:n+1)-X(:,2:n+2)).^2)));
    ud=-ld;
    d=[[0,0,0]',ud]-[ud,[0,0,0]'];
    Yp(:,:,1)=diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1);
    Yp(:,:,2)=diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1); 
    Yp(:,:,3)=diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1);
    derivs=Yp;
end



function Fb=bForce1(X,alpha,BC) %bending forces at each point of deformed rod
    n=size(X);
    n=n(2)-2;
    ph=phi(X); dph=dphi(X);
    L=2*Vdom(X); dl=dL(X);
    Yb=zeros(3,n+2);
    Yb(1,:)= -((2*alpha*ph)./L)*dph(:,:,1) + ((alpha*ph.^2)./L.^2)*dl(:,:,1);
    Yb(2,:)= -((2*alpha*ph)./L)*dph(:,:,2) + ((alpha*ph.^2)./L.^2)*dl(:,:,2);
    Yb(3,:)= -((2*alpha*ph)./L)*dph(:,:,3) + ((alpha*ph.^2)./L.^2)*dl(:,:,3);
    if BC=='r' || BC=='n' 
        %fix left end
        Yb(:,1)=[0,0,0]';
    end
    if BC=='l' || BC=='n'
        %fix right end
        Yb(:,n+2)=[0,0,0]';
    end
    Fb=Yb;
end





function m = qtAngle(X, Twist) %turning angle in quasistatic case, %Twd=twist density
    n=size(X);
    n=n(2)-2;
    Ls=Vdom(X);
    Y=(Twist/sum(Ls(2:n+1))).*Ls;
    m=Y(2:n+1);
end
    
function m = tAngle(X, d1, theta) %turning angle in non-quasistatic case
%d1(i) is reference frame vector, theta(i) is angle on link i 
    n=size(X);
    n=n(2)-2;
    d3=X(:,2:end)-X(:,1:end-1);
    N=cross(d3(:,1:end-1),d3(:,2:end));
    N=N./sqrt(dot(N,N));
    d3=d3./sqrt(dot(d3,d3));
    ph=phi(X);
    ph=ph(2:end-1);
    pTrans=cos(ph).*d1(:,1:end-1)+((1-cos(ph)).*dot(N,d1(:,1:end-1))).*N + (sin(ph).*cross(N,d1(:,1:end-1)));  
    frameTwist=atan2(dot(d1(:,2:end),cross(d3(:,2:end),pTrans)),dot(d1(:,2:end),pTrans));
    m=theta(2:end)-theta(1:end-1)+frameTwist;
end



function Ft = tForce(X,m,beta,BC) %force on x_i due to twisting, m is qtAngle(X,Twist) if quasistatic or tAngle if dynamic
    n=size(X);
    n=n(2)-2;
    Ls=Vdom(X);
    kb=darboux(X);
    Es=sqrt(dot(X(:,2:end)-X(:,1:end-1),X(:,2:end)-X(:,1:end-1))); %Es(i) is distance between x_i and x_i+1
    %forces on x_i due to changing the frame due to changing m at i-1, i+1,
    %i respectively
    Fd=0.5*[[0,0,0]', [0,0,0]', (m(1:end)./Ls(2:end-1)).*(kb(:,1:end)./Es(2:end))]; 
    Fu=0.5*[(m(1:end)./Ls(2:end-1)).*(-kb(:,1:end)./Es(1:end-1)), [0,0,0]',[0,0,0]'];
    Fe=0.5*[[0,0,0]', (m(1:end)./Ls(2:end-1)).*(kb(:,1:end)./Es(1:end-1) - kb(:,1:end)./Es(2:end))  ,[0,0,0]'];
    twist1=-beta*(Fd+Fu+Fe);
    %force on x_i due to changing the lengths of links in the twist energy
    %expression
    L=2*Ls; dl=dL(X); 
    twist2=[((beta*[0,m,0].^2)./L.^2)*dl(:,:,1);((beta*[0,m,0].^2)./L.^2)*dl(:,:,2);((beta*[0,m,0].^2)./L.^2)*dl(:,:,3)];
    Yt=twist1+twist2;
    if BC=='r' || BC=='n' 
        %fix left end
        Yt(:,1)=[0,0,0]';
    end
    if BC=='l' || BC=='n'
        %fix right end
        Yt(:,n+2)=[0,0,0]';
    end
    Ft=Yt;
end


function Fa=angleForce(X,m,beta,BC) %scalar force on the twist angles of each rod segment
    n=size(X);
    n=n(2)-2;
    Ls=Vdom(X);
    Ya=-beta*([0,m./Ls(2:end-1)]-[m./Ls(2:end-1),0]);
    if BC=='r' || BC=='n' 
        %fix left end
        Ya(1)=0;
    end
    if BC=='l' || BC=='n'
        %fix right end
        Ya(n+1)=0;
    end
    Fa=Ya;
end

function d1New = fupdate(X0,X1,d1) %calculate new bishop frame
    d30=X0(:,2:end)-X0(:,1:end-1); d31=X1(:,2:end)-X1(:,1:end-1);
    ph=atan2(sqrt(sum(cross(d30,d31).^2,1)),dot(d30,d31));
    N=cross(d30,d31);
    N=N./sqrt(dot(N,N));
    Y=cos(ph).*d1+((1-cos(ph)).*dot(N,d1)).*N + (sin(ph).*cross(N,d1));  
    Y(isnan(Y))=d1(isnan(Y));
    d1New=Y;
end



function s = stress(X, m, E, I, mutwist, J) %stress squared at each interior point of rod
    n=size(X);
    n=n(2)-2;
    Yb=kappa(X); %curvature
    Yt=[0,m,0]./Vdom(X); %twisting stress
    s = sqrt((E^2)*I*Yb.^2 + mutwist^2*J*Yt.^2);
end





function C = constraints(X,L0) %constraints on rod length under inextensibility
    Y=(1./L0).*dot(X(:,2:end)-X(:,1:end-1),X(:,2:end)-X(:,1:end-1))-L0;
    C=Y';
end

function dC = dCons(X,L0) %derivative of length constraint matrix
    n=size(X);
    n=n(2)-2;
    i=repmat(1:n+1,3,1);
    s=(2./L0).*(X(:,1:end-1)-X(:,2:end));
    dC=sparse([i(:);i(:)],[(1:3*(n+1))',(4:3*(n+2))'],[s(:);-s(:)], n+1, 3*(n+2));
end


function dX = project(X,L0) %manifold projection
    n=size(X);
    n=n(2)-2;
    C=constraints(X,L0);
    dC=dCons(X,L0);
    A=dC*dC';
    lambda=-A\C;
    dX=dC'*lambda;
end



function [x,y,z,C]=tubeplot2(curve,r,n,ct,S)
% Usage: [x,y,z]=tubeplot(curve,r,n,ct,S)
% 
% Modified version of tubeplot from matlab file exchange. (see license at
% top)
% Tubeplot constructs a tube, or warped cylinder, along
% any 3D curve, much like the build in cylinder function.
% If no output are requested, the tube is plotted.
% Otherwise, you can plot by using surf(x,y,z);
%
%
% Arguments:
% curve: [3,N] vector of curve data
% r      the radius of the tube
% n      number of points to use on circumference. Defaults to 8
% ct     threshold for collapsing points. Defaults to r/2 
% S      scalar field plotted along curve as colour
%
% The algorithms fails if you have bends beyond 90 degrees.
% Janus H. Wesenberg, july 2004

  if nargin<3 || isempty(n), n=8;
     if nargin<2, error('Give at least curve and radius');
    end;
  end;
  if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
  end;
  if nargin<4 || isempty(ct)
    ct=0.5*r;
  end

  
  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct;
      npoints=npoints+1;
      curve(:,npoints)=curve(:,k);
    end
  end
  %Always include endpoint
  if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
  end

  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);

  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;

  xyz=repmat([0],[3,n+1,npoints+2]);
  Col=repmat([0],[3,n+1,npoints+2]); 

  %precalculate cos and sing factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
  
  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
    Col(:,:,k+1)=S(k);
  end;
  %finally, cap the ends:
  xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
  xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  Ct=squeeze(Col(1,:,:));
  C=Ct;
  %... and plot:
  if nargout<3, surf(x,y,z,C); end
end
















