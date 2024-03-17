% random walks

        clear all;
        close all
        npart = 1000;
        random.low=-5.;
		random.high=5.;
		random.delta=random.high-random.low;
        
        
        for i=1:1000
        r = random.low + random.delta*rand(npart,1);
        
        w=zeros(npart,1);
        
        w(1)=r(1);
        
        for x=2:npart
        w(x)=w(x-1)+r(x);
        end
        
        y(i)=w(npart);
        plot(w)
        hold on;
        
        end
        
        figure (2)
        histfit(y,20)
        
        
        
% take a random walk to find the slope of this regression       

x=[.1:.1:100];
x=x(:);   % this transposes x, as would x' or transpose(x)

% simple linear regression model

% y = b(1) + b(2) * x;

b1=100;
b2=50;

y=b1+ b2 * x;

% add a random component that is about 10% of the respective values of y

yrand=randn(1000,1);  % normal distribution, mean of zero, n by 1 matrix

yerr= yrand .* y *.1 ;  % apply error that is proportion to 0.1 cv

yy=y+ yerr;

figure(1)
clf
plot(x,yy,'.')

xlabel('x')
ylabel('y')
title('y = 50 x + noise')

h = lsline;
set(h(1),'color','k','LineWidth',1.5)


% prior for slope..assume between 0 and 100
      
        % prior for slope, assume 0 to 100
        % sample from a uniform distribution
        
        rng;
        rslope.low=0;
		rslope.high=100;
		rslope.delta=rslope.high-rslope.low;
        rslope.rand = rslope.low + rslope.delta*rand(1,1);
        
        
        % prior for intercept, -50 to 50
        % sample from a uniform distribution
        
        rint.low=-50;
		rint.high=50;
		rint.delta=rint.high-rint.low;
        rint.rand = rint.low + rint.delta*rand(1,1);
        
        
        n=length(x);
        
        
        
        % take random walk and find least cost of beta
        beta=[rint.rand(1), rslope.rand(1)];
        
      cnt1=0;
      cnt2=0;
        
        for i=1:10000
           
            
        ypoly = beta(1) + beta(2)*x;   
        
        se=std(yy)/sqrt(n);
        
        y1 = normpdf(mean(ypoly), mean(yy),se); % likelihood function
           
        
        
        % pick new beta and perturb it with persistance from tha past and a
        % random forcing
        rint.rand = rint.low + rint.delta*rand(1,1);
        
        %betatst= beta+ 0.5* r * beta;
        
        betatst(1)= 0.5* (beta(1)+ rint.rand * beta(1));
        
        ypolytst = betatst(1) + beta(2)*x;  
        y2 = normpdf(mean(ypolytst), mean(yy),se); 
        
        while y2==0
        y2 = normpdf(mean(ypolytst), mean(yy),se); % likelihood function
         rint.rand = rint.low + rint.delta*rand(1,1);
        beta(1)= 0.5* (beta(1)+ rint.rand);
         ypolytst = betatst(1) + beta(2)*x;  
          
        end
        
        
        
        
        if y2 > y1  % if beta for likelihood y2 is greater than y1, then accept beta for y2
        beta(1)=betatst(1);
        cnt1=cnt1+1;
        bb1(cnt1)=beta(1);
        end
        
        rslope.rand = rslope.low + rslope.delta*rand(1,1);
        betatst(2)= 0.5* (beta(2)+ rslope.rand);
        % perturb beta +/- 1% of its value
        
        ypolytst = beta(1) + betatst(2)*x;  
        y3 = normpdf(mean(ypolytst),mean(yy),se); % likelihood function
        
        if y3 > y2  % if beta for likelihood y2 is greater than y1, then accept beta for y2
        beta(2)=betatst(2);
        cnt2=cnt2+1;
         bb2(cnt2)=beta(2);
        end
        
       
                
        if y3 == 0
             beta(2)= 0.5* (beta(2)+ rslope.rand);
        end
        
        
       
            
            
            
        end
                    
        
        
        
             
        
        figure(5)
        clf
      semilogx(bb1)
      hold on
      semilogx(bb2)