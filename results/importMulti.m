function [N_vals,T_vals,data] = importMulti(filename)

numOfTemps=readmatrix(filename,"FileType","text","Delimiter","\t","Range",'B1:B1');
numOfIterations=readmatrix(filename,"FileType","text","Delimiter","\t","Range",'B2:B2');
%it=numOfTemps*numOfIterations;

import=readmatrix(filename,"FileType","text","Delimiter","\t","Range",[5,1,1E8,4]); %assume file is never longer than 1E8 rows

N_vals=import(1,1);
T_vals(1,1)=import(1,2);
data=NaN(5,numOfTemps,3,numOfIterations); %this is maybe too big, but is reduced at the end

n=1;
t=1;
dataIndex=1;
for i=1:size(import,1)
   if import(i,1)~=N_vals(n)
      n=n+1; %later we can search if this N have been added earlier 
      t=0; %assume last N and next N donst end and start with same temp
            
      N_vals(n)=import( i,1); % energy
   end
   
   if (t==0)||(import(i,2)~=T_vals(t))
      t=t+1; %later we can search if this t has been found earlier
      
      dataIndex=1;
      T_vals(t,n)=import( i,2); % temps
   end
   data(n,t,1,dataIndex)=import( i,4); % energy
   data(n,t,2,dataIndex)=abs(import( i,3)); % absolute magnetization
   data(n,t,3,dataIndex)=import( i,2); % temp again, only for convenience
   
   dataIndex=dataIndex+1;
    
end

data(isnan(data))=[];% remove NaN's 

%% Check integrity
%if(max(N)~=min(N))
%    error("wrong separation of data! (correct dataLineStart and dataLineEnd)")
%end
end