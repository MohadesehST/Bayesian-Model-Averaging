clc
clear
close all
S=xlsread('e:\w_july');    %input data
D=S(:,1:7);                        %Ensemble of forecastes
Y=S(:,8);                            %Ensemble of observation

% In the original BMA algorithm, it is assumed that posterior distribution is Gaussian. This assumption,
% however, is not applicable to all types of forecast variables, such as stream flow data and evaporation data.
% Therefore, the Box-Cox transformation must be used to map both model prediction and observed data from
% their original spaces to Gaussian spaces .

y=boxcox(Y);
[transdat,lambda]=boxcox(Y);
[n,k]=size(D);       

 % Bias correction forecasts ensemble and observation
  for i=1:7              
   x=[ones(length(D(:,i)),1) D(:,i)];  
   b(i,:)=x\Y;                  
   bcD(:,i)=boxcox(D(:,i));       
   d(:,i)=b(i,1)+b(i,2)*D(:,i);    
   yhat(:,i)=b(i,1)+b(i,2)*bcD(:,i);
   mu(:,i)=mean(yhat(:,i));
   sd(:,i)=std(yhat(:,i));
   yhat22(:,i)=normpdf(y,mu(i),sd(i));     % Normal pdf of models
  end
  
%%Optimization method: Expectation-Maximization algorithm
  
  w=ones(1,k)/k;                  % Initial weights      
 for t=1:18            
     for k=1:7
         s(k,t)=(y(t)-yhat(t,k)).^2; 
         s1(t,k)=s(k,t)';
     end
    s2=sum(s1);
 end
  M=sum(s2);
  sigma2=(1/k)*(1/n)*M*ones(1,k);    %Initial Variances
  sigma=sqrt(sigma2);       
  
 tt=1;
 max_tt=100;
 while(tt<=max_tt)
 
      if tt==1    %%first step
      for k=1:7
          z(:,k)=w(k)*normpdf(y,mu(k),sigma(k));
      end
      z_t=bsxfun(@rdivide,z,sum(z,2)); 
      w_t=sum(z_t)/18;
      l(tt)=sum(log(sum(z,2)));            
      sigma2_t= sum(z_t.*bsxfun(@minus,yhat,y).^2)./sum(z_t);
      sd1=sqrt(sigma2_t);
      end
    if tt>1
      for k=1:7
           zt(:,k)=w_t(k)*normpdf(y,mu(k),sd1(k));
      end
      z_tt=bsxfun(@rdivide,zt,sum(zt,2)); 
      w_tt=sum(z_tt)/18;
      l(tt)=sum(log(sum(zt,2)));
      sigma2_tt= sum(z_tt.*bsxfun(@minus,yhat,y).^2)./sum(z_tt);
      w_t=w_tt;
   weight_total(tt,1:k)=w_t;
   sd1 = sqrt(sigma2_tt);
    err(tt-1)=abs(l(tt)-l(tt-1));
    if err(tt-1)<0.1
        break
    end
    end
      tt= tt+1;
 end
[u,v]=find(err<0.1);   %%tt_test

for k=1:7
    yhat2(:,k)=normpdf(y,mu(k),sd1(k));
     pdfkol(k,:)=w_tt(k)*yhat2(:,k);
end
t=1:1:18;
for k=1:7
    pdfkol1(k,:)=w_tt(k)*d(:,k);      
end
totalsim=sum(pdfkol1);
