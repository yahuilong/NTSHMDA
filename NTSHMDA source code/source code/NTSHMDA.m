function NTSHMDA(gamma,phi,delte,beta1,beta2)
  %NTSHMDA(0.7,0.9,0.3,0.8,0.2)
  %random walk with restart on the heterogeneous network to calculate association scores for each microbe-disease pair.
  %gamma: restart probability
  %phi: jumping probability
  %delte: weight factor
  %beta1, beta2: weight factor 
  %MS: adjacency matrix of the microbe similarity network
  %DS: adjacency matrix of the disease similarity network
  %interaction: adjacency matrix of the microbe-disease association network
  %A: known disease microbe interaction
  
  A=textread('knowndiseasemicrobeinteraction.txt');
  [pp,qq]=size(A);
   
  load interaction;
  [nd nm]=size(interaction);
  %calculate microbe similarity and disease similarity;
  [DS,MS]=gaussiansimilarity(interaction,nd,nm);
      
  %reobtain adjacency matrix
  Am=interaction*MS;
  Ad=DS*interaction;
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  based on Am   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
      %establish heterogeneous network based on reobtained adjacent matrix Am
      HN1=ones(nd+nm);
      HN1(1:nd,1:nd)=DS;
      HN1(1:nd,(nd+1):(nd+nm))=Am;
      HN1((nd+1):(nd+nm),1:nd)=Am';
      HN1((nd+1):(nd+nm),(nd+1):(nd+nm))=MS;
  
      %the transition probability
      Am=Am.*interaction;
      W_md1=row_norm(Am');        %the transition probability from microbe nodes to disease nodes
      W_md1=W_md1';
      W_dm1=row_norm(Am);         %the transition probability from disease nodes to microbe nodes
      W_dm1=W_dm1';
     
  
      W_dd=row_norm(DS);         %the transition probability from disease nodes to disease nodes
      W_mm=row_norm(MS);         %the transition probalibitly from microbe nodes to microbe nodes
 
      %transition probability matrix
      W1=HN1;
      
      %calculate disease-disease transition probability matrix
      for i=1:nd                  
         if sum(interaction(i,:))==0
              W1(i,1:nd)=W_dd(i,:);
         else
              W1(i,1:nd)=(1-phi)*W_dd(i,:);
         end
      end
      
      W1(1:nd,(nd+1):(nd+nm))=phi*W_md1;
      W1((nd+1):(nd+nm),1:nd)=phi*W_dm1;
      
      %calculate microbe-microbe transition probability matrix
      for j=1:nm
         if sum(interaction(:,j))==0
              W1(nd+j,nd+1:nd+nm)=W_mm(j,:);
         else
              W1(nd+j,nd+1:nd+nm)=(1-phi)*W_mm(j,:);
         end
      end
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  based on Ad  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     %establish heterogeneous network based on reobtained adjacent matrix Ad
      HN2=ones(nd+nm);
      HN2(1:nd,1:nd)=DS;
      HN2(1:nd,(nd+1):(nd+nm))=Ad;
      HN2((nd+1):(nd+nm),1:nd)=Ad';
      HN2((nd+1):(nd+nm),(nd+1):(nd+nm))=MS;
  
      %the transition probability
      Ad=Ad.*interaction;
      W_md2=row_norm(Ad');        %the transition probability from microbe nodes to disease nodes
      W_md2=W_md2';
      W_dm2=row_norm(Ad);         %the transition probability from disease nodes to microbe nodes
      W_dm2=W_dm2';
    
  
      %W_dd=row_norm(DS);         %the transition probability from disease nodes to disease nodes
      %W_mm=row_norm(MS);         %the transition probalibitly from microbe nodes to microbe nodes
 
      %transition probability matrix
      W2=HN2;
      
      %calculate disease-disease transition probability matrix
      for i=1:nd                  
         if sum(interaction(i,:))==0
              W2(i,1:nd)=W_dd(i,:);
         else
              W2(i,1:nd)=(1-phi)*W_dd(i,:);
         end
      end
      
      W2(1:nd,(nd+1):(nd+nm))=phi*W_md2;
      W2((nd+1):(nd+nm),1:nd)=phi*W_dm2;
      
      %calculate microbe-microbe transition probability matrix
      for j=1:nm
         if sum(interaction(:,j))==0
              W2(nd+j,nd+1:nd+nm)=W_mm(j,:);
         else
              W2(nd+j,nd+1:nd+nm)=(1-phi)*W_mm(j,:);
         end
      end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      
      %initial walk probability matrix
      M0=eye(nd+nm); %p0 is a matrix, representing the initial probability matrix
      M0(1:nd,1:nd)=(1-delte)*(1/nd).*M0(1:nd,1:nd);
      M0(nd+1:nd+nm,nd+1:nd+nm)=delte*(1/nm).*M0(nd+1:nd+nm,nd+1:nd+nm);
      Pt1=M0;
      Pt2=M0;
      
      %random walk on the heterogeneous network
      for i=1:10  
         P1 = (1-gamma) * W1 * Pt1 + gamma * M0;
         P2 = (1-gamma) * W2 * Pt2 + gamma * M0;
         Pt1=P1;
         Pt2=P2;
      end
     
      P=(beta1*P1+beta2*P2)/(beta1+beta2);
      %save ('Result.mat','P');
       
      %obtain predicted disease-microbe association scores;
      prediction_score=P(1:nd,nd+1:nd+nm);

end