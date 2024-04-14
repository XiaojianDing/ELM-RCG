% clear
% clc
 eps1=1.0000e-08;
%eps1=1;


load iono;

% rand_sequence=randperm(size(iono,1));
% temp_dataset=iono;
% iono=temp_dataset(rand_sequence, :);


P=iono(1:100,2:(size(iono,2)))';
T=iono(1:100,1)';

pos5=find(T==0);
T(pos5)=-1;



P1=iono(101:size(iono,1),2:(size(iono,2)))';
T1=iono(101:size(iono,1),1)';
pos6=find(T1==0);
T1(pos6)=-1;


xapp=P';
yapp=T';

xtest=P1';
ytest=T1';
 
 


n=size(P,2);
n1=size(P1,2);
dim=size(P,1);

 % C=a1;
 
 C=2^(2);

 
 
  kerneloption=rand(1,dim)*dim;


 
%  C=a1(aa1);
% n2=a2(aa2);
%  
 

 
 %--------------------------random compact kernel---------------------
 
 epsilon = .000001;
kernel='rcg';
verbose = 0;


 tic; 
 
%ps  =  zeros(n,n);		
ps = randkernel(xapp,kernel,kerneloption);

H =ps.*(yapp*yapp'); 
 
 % t=toc;
 
 % tic
 
 %---------------------------------solve the QP problem------------------%
 f=ones(n,1);


%[alpha_pos , lambda , position] = monqp(K1,f,zeros(n,1),0,C,eps1,1,P',K,[]);      
  [alpha_pos, position,iter] =monqp_ELM(H,f,C,eps1,0);    
t=toc;


xsup = xapp(position,:);
ysup = yapp(position);
w = (alpha_pos.*ysup);

Hsup=H(:,position);
ps1=H'*Hsup;
y=ps1*w;

       % actual output

tic;
%------------------------Compute test data output----------------%

ps2=randkernel(xtest,kernel,kerneloption,xsup);
y1=ps2*w;


t1=toc;

     err=0;
   for i=1:n
    if(T(i)==1)
        if(y(i)<=0)
            err=err+1;
        end
    else
        if(y(i)>=0)
            err=err+1;
        end
    end
   end
   
   accu=(n-err)/n;
   
   %----------------------Error of test data--------------------------------%
    err1=0;
   for i=1:n1
    if(T1(i)==1)
        if(y1(i)<=0)
            err1=err1+1;
        end
    else
        if(y1(i)>=0)
            err1=err1+1;
        end
    end
   end
   
   accu1=(n1-err1)/n1
   

   