for sub=1:21
    sub
    clear y ye x Yr  err
    y=squeeze(data(:,:,sub));
    p=2;
    UC=0.001;
    mode=0;
%     [x2,e,Kalman2,Q2]=mvaar(squeeze(data(:,:,sub)),bic_smooth,ff_actual,0,Kalman2);
y=single(y);

[M,LEN] = size(y');		%number of channels, total signal length
L = M*M*p;


ye = zeros(size(y));	%prediction of y


        x=zeros(L,LEN);


%upd = eye(L,'single')/single(L*UC);		%diagonal matrix containing UC

for n = 2:LEN,
    
        if(n<=p)
                Yr=[y(n-1:-1:1,:)' zeros(M,p-n+1)];	%vector of past observations
                Yr=single(Yr(:)');
        else
                Yr=y(n-1:-1:n-p,:)';						%vector of past observations
                Yr=single(Yr(:)');
        end
        
        %Update of measurement matrix
        Kalman.H=kron(eye(M,'single'),Yr);
        
        
        %calculate prediction error
        ye(n,:)=(Kalman.H*Kalman.x)';
        err=y(n,:)-ye(n,:);
        
        if ~any(isnan(err(:))),
                %update of Q2 using the prediction error of the previous step
                Kalman.Q2=(1-UC)*Kalman.Q2+UC*(err'*err);
                
                
                KpH=Kalman.Kp*Kalman.H';
                % HKp=Kalman.H*Kalman.Kp;
                
                %Kalman gain
                Kalman.G=(KpH)/(Kalman.H*(KpH)+Kalman.Q2);
               
                %calculation of the a-posteriori state error covariance matrix
                %K=Kalman.Kp-Kalman.G*KpH'; Althouh PK is supposed to be symmetric, this operation makes the filter unstable
                %K=Kalman.Kp-Kalman.G*HKp; 
                
                %mode==0 no update of Q1
                %update of Q1 using the predicted state error cov matrix
                if(mode==1)      
                        Kalman.Q1=diag(diag(K)).*UC;
                elseif(mode==2)
                        Kalman.Q1=upd*trace(K);
                elseif(mode==3)
                        Kalman.Q1=diag(sum((Block*diag(diag(K)))'))/(p*M)*UC;
                elseif(mode==4)
                        avg=trace(K(index,index))/(p*M)*UC;
                        Kalman.Q1=Block1*UC+Block0*avg;
                end
                
                %a-priori state error covariance matrix for the next time step
                Kalman.Kp=Kalman.Kp-Kalman.G*(Kalman.H*Kalman.Kp);
                for t=1:L, Kalman.Kp(t,t)=Kalman.Kp(t,t)+single(UC);end
                
                %current estimation of state x
                Kalman.x=Kalman.x+Kalman.G*(err)';
        end; % isnan>(err)   
         x(:,n) = Kalman.x;
       
end
e = y - ye;
x = x';
temp1=reshape(x,size(data,1),size(data,2),size(data,2),bic_smooth);

    conn_individual(:,:,:,sub)=squeeze(sum(temp1,4));
    clear x temp1 Kalman;
    load('Kalman.mat');
end
