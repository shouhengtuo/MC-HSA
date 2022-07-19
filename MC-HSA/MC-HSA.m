clear;
clc
tic

%%%%%%%%%%%% Predefine the parameters of MCHS algorithm %%%%%%%%%%%%%%%%%% 
D = 5;  %%%The dimension of solution
HMS1 = 50;%%%The size of HM1,HM2,HM3 in three elementary membranes
HMS2 = 50;%%%The size of HM4 in membrane "4"
%HMS2 = 100;
HMS3 = 100;%%%The size of HM5 in membrane "5"
%HMS3 = 100;
MaxFEs = 100000;%%%The maximum number of iterations in the searching stage
%MaxFEs = 150000;
lambda = 0.1;%%% penalty factor
theta=0.95;%%%A threshold in the verifying stage

%%%%%%%%%%%% Predefine the parameters of Harmony search algorithm %%%%%%%%%%%%%%%%%%
HMCR = 0.98;
PAR = 0.35;
BW = 10;

%%%%%%% Upload the dataset %%%%%
% x = load('filteredGeneExpression.mat');
% data = x.diffData';
% label = x.labels;
% x = load('filteredMethylation.mat');
% data = [x.downCgData;x.upCgData]';
% label = x.LABELS';
% [~,n] = size(data);
% xL = 1;
% xU = n;


%%%%%%%%%%%%% Data balance %%%%%%%%%%
%  data = [data(1:16,:);data];
%  label = [label(1:16,:);label]; 
%  sId = randperm(48);
%  data = data(sId,:);
%  label = label(sId,:);
 load('filteredleukemia.mat');
data=p2_leukemia(1:72,1:1970);
label=p2_leukemia(1:72,1971);
  
   data1=p2_leukemia;
  [label1,ind]=sort(data1(:,1971));
  data2=data1(ind,:);
  r=[];
  for i=1:22
     ran=randi([48,72]);
     while ismember(ran,r)
        ran=randi([48,72]);
     end
     r=[r;ran];
     data3(i,:)=data2(ran,:);
  end
data3=[data2;data3];
data=data3(1:94,1:1970);
label=data3(1:94,1971);


%%%%%%% Randomly initialize the solutions of three elementary membranes
for k = 1:3
    for i = 1:HMS1
        X0(i,:) = rand(1,D).*(xU-xL) + xL;
        HM(k).X(i,:) = sort(round(X0(i,:)));
        SX = HM(k).X(i,:); 
        X = length(SX);
        SX = unique(SX);
        XL = length(SX);
        penalty = XL/X;
        HM(k).Fit(i,:)=multiObjFitness(data(:,SX),label);
        
        %%%%Calculate the fitness of solutions in each membrane
        HM(k).Fit_penalty(i,:) = HM(k).Fit(i,:) - lambda * penalty;
    end
end

Fes=0;
gen=0;
EliteNum2=0;
EliteNum=0;
OutEliteNum=0;
HM(4).X=[];
HM(4).Fit=[];
HM(4).Fit_penalty=[];
HM(5).X=[];
HM(5).Fit=[];
HM(5).Fit_penalty=[];
HM(6).X=[];
HM(6).Fit=[];
HM(6).Fit_penalty=[];

%%%%%  MCHS in the searching stage
while Fes < MaxFEs
     gen=gen+1;
     %%%% Information exchange in the searching stage
    if mod(gen,HMS1)==1
        if Fes < 0.8 * MaxFEs %%%%Information exchange between membrane "4"and three elementary membranes
             for k=1:3
                [HM(k).X,idx,~]=unique(HM(k).X,'rows');
                HM(k).Fit = HM(k).Fit(idx,:);
                HM(k).Fit_penalty = HM(k).Fit_penalty(idx,:);
                [n,~]=size(HM(k).X);
                while n<HMS1
                     n=n+1;
                     X0(n,:)=rand(1,D).*(xU-xL) + xL;
                     HM(k).X(n,:) = sort(round(X0(n,:)));
                     SX = HM(k).X(n,:); 
                     X = length(SX);
                     SX = unique(SX);
                     XL = length(SX);
                     penalty = XL/X;
                     HM(k).Fit(n,:)=multiObjFitness(data(:,SX),label);
                     HM(k).Fit_penalty(n,:) = HM(k).Fit(n,:) - lambda * penalty;
                end
                 [~,Kbestid]=max(HM(k).Fit_penalty(:,1));
                 [~,Sbestid]=max(HM(k).Fit_penalty(:,2));
                 [~,Gbestid]=max(HM(k).Fit_penalty(:,3));
                 [~,Tbestid]=max(HM(k).Fit_penalty(:,4));
                 if Kbestid==Sbestid && Sbestid==Gbestid && Gbestid==Tbestid
                     fbest_penalty=HM(k).Fit_penalty(Kbestid,:);
                     bestId=Kbestid;
                 else
                     [fbest_penalty,bestId] = max(sum(HM(k).Fit_penalty,2));
                 end
                xbest = HM(k).X(bestId,:);
                Fbest_penalty = HM(k).Fit_penalty(bestId,:);
                Fbest = HM(k).Fit(bestId,:);
                 % The best solutions of three elementary membranes flow into membrane "4"
                 if isempty(HM(4).X)||~ismember(xbest,HM(4).X,'rows') 
                    if EliteNum < HMS2    
                       HM(4).X = [HM(4).X; xbest];
                       HM(4).Fit = [HM(4).Fit; Fbest];
                       HM(4).Fit_penalty = [HM(4).Fit_penalty; Fbest_penalty];
                       EliteNum = EliteNum + 1;
                             X0(bestId,:) = rand(1,D).*(xU-xL) + xL;
                %              X(bestId,:) = sort(round(X0));
                             HM(k).X(bestId,:) = sort(round(X0(bestId,:)));
                              SX = HM(k).X(bestId,:); 
                                X = length(SX);
                                SX = unique(SX);
                                XL = length(SX);
                                penalty = XL/X;
                                HM(k).Fit(bestId,:)=multiObjFitness(data(:,SX),label);
                                HM(k).Fit_penalty(bestId,:) = HM(k).Fit(bestId,:) - lambda * penalty;
                             Fes = Fes + 1;
                    else
                        [~,Kworstid]=min(HM(4).Fit_penalty(:,1));
                        [~,Sworstid]=min(HM(4).Fit_penalty(:,2));
                        [~,Gworstid]=min(HM(4).Fit_penalty(:,3));
                        [~,Tworstid]=min(HM(4).Fit_penalty(:,4));
                        if Kworstid==Sworstid && Sworstid==Gworstid && Gworstid==Tworstid
                           eworstFit_penalty=HM(4).Fit_penalty(Sworstid,:);
                           worstId=Sworstid;
                        else
                           [eworstFit_penalty,worstId] = min(sum(HM(4).Fit_penalty,2));
                        end
                        if (eworstFit_penalty(1) <= Fbest_penalty(1)) && (eworstFit_penalty(2) <= Fbest_penalty(2)) && (eworstFit_penalty(3) <= Fbest_penalty(3)) && (eworstFit_penalty(4) <= Fbest_penalty(4))
                            if ~ismember(xbest,HM(4).X,'rows')
                                HM(k).X(bestId,:) = HM(4).X(worstId,:);
                                HM(k).Fit(bestId,:)=HM(4).Fit(worstId,:);
                                HM(k).Fit_penalty(bestId,:)=HM(4).Fit_penalty(worstId,:);
                                
                                HM(4).X(worstId,:) = xbest;
                                HM(4).Fit(worstId,:) = Fbest;
                                HM(4).Fit_penalty(worstId,:) = Fbest_penalty;
                            elseif ~ismember(HM(4).X(worstId,:),HM(k).X,'rows')
                                HM(k).X(bestId,:) = HM(4).X(worstId,:);
                                HM(k).Fit(bestId,:)=HM(4).Fit(worstId,:);
                                HM(k).Fit_penalty(bestId,:)=HM(4).Fit_penalty(worstId,:);
                            end
                        end            
                    end
                 end
             end
        end
         [~,eKbestid]=max(HM(4).Fit_penalty(:,1));
         [~,eSbestid]=max(HM(4).Fit_penalty(:,2));
         [~,eGbestid]=max(HM(4).Fit_penalty(:,3));
         [~,eTbestid]=max(HM(4).Fit_penalty(:,4));
        if eKbestid==eSbestid && eSbestid==eGbestid && eGbestid==eTbestid
            efbest_penalty=HM(4).Fit_penalty(eSbestid,:);
            ebestId=eSbestid;
        else
            [efbest_penalty,ebestId] = max(sum(HM(4).Fit_penalty,2));
        end
           
        eXbest = HM(4).X(ebestId,:);%将最好的x找出来
        eFbest_penalty = HM(4).Fit_penalty(ebestId,:);
        eFbest = HM(4).Fit(ebestId,:);
%%% Information exchange between membrane "4" and membrane "5" 
        xbest2 = unique(eXbest);  
        if EliteNum == HMS2
                    if  isempty(HM(5).X) || EliteNum2 < HMS3  && ~ismember(eXbest,HM(5).X,'rows')
                              HM(5).X = [HM(5).X; eXbest];
                               HM(5).Fit = [HM(5).Fit; eFbest];
                               HM(5).Fit_penalty = [HM(5).Fit_penalty; eFbest_penalty];
                               EliteNum2 = EliteNum2 + 1;

                               HM(4).X(ebestId,:)=[];
                               HM(4).Fit(ebestId,:) = [];
                               HM(4).Fit_penalty(ebestId,:) = [];
                               EliteNum = EliteNum - 1;
                    else
                             [~,eKworstid]=min(HM(5).Fit_penalty(:,1));
                             [~,eSworstid]=min(HM(5).Fit_penalty(:,2));
                             [~,eGworstid]=min(HM(5).Fit_penalty(:,3));
                             [~,eTworstid]=min(HM(5).Fit_penalty(:,4));
                             if eKworstid==eSworstid && eSworstid==eGworstid && eGworstid==eTworstid
                                 eworstFit_penalty=HM(5).Fit_penalty(eSworstid,:);
                                 worstId=eKworstid;
                             else
                                [eworstFit_penalty,worstId] = min(sum(HM(5).Fit_penalty,2)); 
                             end
                                if (eworstFit_penalty(1) <= eFbest_penalty(1)) && (eworstFit_penalty(2) <= eFbest_penalty(2)) && (eworstFit_penalty(3) <= eFbest_penalty(3)) && (eworstFit_penalty(4) <= eFbest_penalty(4))
                                    if ~ismember(eXbest,HM(5).X,'rows')
                                        HM(4).X(ebestId,:) = HM(5).X(worstId,:);
                                        HM(4).Fit(ebestId,:) = HM(5).Fit(worstId,:);
                                        HM(4).Fit_penalty(ebestId,:) = HM(5).Fit_penalty(worstId,:);
                                        HM(5).X(worstId,:) = eXbest;
                                        HM(5).Fit(worstId,:) = eFbest; 
                                        HM(5).Fit_penalty(worstId,:) = eFbest_penalty;
                                    elseif ~ismember(HM(5).X(worstId,:),EliteX,'rows')
                                        HM(4).X(ebestId,:) = HM(5).X(worstId,:);
                                        HM(4).Fit(ebestId,:) = HM(5).Fit(worstId,:);
                                        HM(4).Fit_penalty(ebestId,:) = HM(5).Fit_penalty(worstId,:);
                                    end
                                end            
                         end
        end

                
                 
    %%%%%%% Optimize the solution in three elementary membranes 
     if Fes < MaxFEs*0.8  
           for k = 1:3
                for i = 1:D
                    if rand < HMCR
                        a = ceil(rand*D);
                        Xnew(i) = HM(k).X(ceil(rand*HMS1),a);
                        if rand <PAR
                            Xnew(i) = round(Xnew(i) + BW*(rand-0.5));
                            Xnew(i) = max(min(Xnew(i),xU),xL);
                        end
                    else
                        Xnew(i) = round(xL + rand*(xU - xL));
                    end            
                end
                Xnew = sort(Xnew);
                 SX = Xnew; 
                X = length(SX);
                SX = unique(SX);
                XL = length(SX);
                penalty = XL/X;
                fnew=multiObjFitness(data(:,SX),label);
                Fit_penaltyNew = fnew - lambda * penalty;
                
                  Fes = Fes + 1;
                  [~,Kworstid]=min(HM(k).Fit_penalty(:,1));
                 [~,Sworstid]=min(HM(k).Fit_penalty(:,2));
                 [~,Gworstid]=min(HM(k).Fit_penalty(:,3));
                 [~,Tworstid]=min(HM(k).Fit_penalty(:,4));
                 if Kworstid==Sworstid && Sworstid==Gworstid && Gworstid==Tworstid
                     Fworst_penalty=HM(k).Fit_penalty(Kworstid,:);
                     idworst=Kworstid;
                 else
                     [Fworst_penalty, idworst] = min(sum(HM(k).Fit_penalty,2));
                 end 
                

                    if ((Fit_penaltyNew(1)) >= HM(k).Fit_penalty(idworst,1)) && ((Fit_penaltyNew(2)) >= HM(k).Fit_penalty(idworst,2)) && ((Fit_penaltyNew(3)) >= HM(k).Fit_penalty(idworst,3)) && ((Fit_penaltyNew(4)) >= HM(k).Fit_penalty(idworst,4))%%%选取准确率较高的
                        HM(k).Fit(idworst,:) = fnew;
                        HM(k).Fit_penalty(idworst,:) = Fit_penaltyNew;
                        HM(k).X(idworst,:) = Xnew;                        
                    end
           end
     end
   
     %%%%% Optimize the solutions in HM4 and HM5
    if gen > 2*HMS1
        
         for i = 1:D
            if rand < HMCR
                 a = ceil(rand*D);
                Xnew(i) = HM(4).X(ceil(rand*HMS2),a);
                if rand <PAR
                    Xnew(i) = round(Xnew(i) + BW*(rand-0.5));
                    Xnew(i) = max(min(Xnew(i),xU),xL);
                end
            else
                Xnew(i) = round(xL + rand*(xU - xL));
            end            
        end
        Xnew = sort(Xnew);
        while  ismember(Xnew,HM(4).X,'rows')
            a = ceil(rand*D);
            Xnew(a) = round(xL + rand*(xU - xL));            
        end
        Xnew = sort(Xnew);
        
         SX = Xnew; 
        X = length(SX);
        SX = unique(SX);
        XL = length(SX);
        penalty = XL/X;
        efnew =multiObjFitness(data(:,SX),label);
        Fit_penaltyNew = efnew - lambda * penalty;
              Fes = Fes + 1;
             [~,Kworstid]=min(HM(4).Fit_penalty(:,1));
             [~,Sworstid]=min(HM(4).Fit_penalty(:,2));
             [~,Gworstid]=min(HM(4).Fit_penalty(:,3));
             [~,Tworstid]=min(HM(4).Fit_penalty(:,4));
                 if Kworstid==Sworstid && Sworstid==Gworstid && Gworstid==Tworstid
                     eFworst_penalty=HM(4).Fit_penalty(Kworstid,:);
                     eidworst=Kworstid;
                 else
                     [eFworst_penalty, eidworst] = min(sum(HM(4).Fit_penalty,2));
                 end   
          

             if (Fit_penaltyNew(1) >= HM(4).Fit_penalty(eidworst,1)) && (Fit_penaltyNew(2) >= HM(4).Fit_penalty(eidworst,2)) &&(Fit_penaltyNew(3) >= HM(4).Fit_penalty(eidworst,3)) && (Fit_penaltyNew(4) >= HM(4).Fit_penalty(eidworst,4))  %%%选取准确率较高的
                 HM(4).Fit(eidworst,:) = efnew;
                 HM(4).X(eidworst,:) = Xnew;  
                 HM(4).Fit_penalty(eidworst,:) = Fit_penaltyNew;
             end
        
         if EliteNum2>0
                 for i = 1:D
                    if rand < HMCR
                         a = ceil(rand*D);
                        Xnew(i) = HM(5).X(ceil(rand*EliteNum2),a);
                        if rand <PAR
                            Xnew(i) = round(Xnew(i) + BW*(rand-0.5));
                            if Xnew(i) > xU 
                               Xnew(i) = 2*xU - Xnew(i);
                            elseif Xnew(i) < xL
                               Xnew(i) = 2*xL - Xnew(i); 
                            end
                        end
                    else
                        Xnew(i) = round(xL + rand*(xU - xL));
                    end            
                end
                Xnew = sort(Xnew);            
                while  ismember(Xnew,HM(5).X,'rows')
                    a = ceil(rand*D);
                    Xnew(a) = round(xL + rand*(xU - xL));            
                end
                Xnew = sort(Xnew);
                SX = Xnew; 
                X = length(SX);
                SX = unique(SX);
                XL = length(SX);
                penalty = XL/X;
                efnew =multiObjFitness(data(:,SX),label);
                Fit_penaltyNew = efnew - lambda * penalty;  
                minScore = min(Fit_penaltyNew);
                Fes = Fes + 1;
                 [~,Kworstid]=min(HM(5).Fit_penalty(:,1));
                 [~,Sworstid]=min(HM(5).Fit_penalty(:,2));
                 [~,Gworstid]=min(HM(5).Fit_penalty(:,3));
                 [~,Tworstid]=min(HM(5).Fit_penalty(:,4));
                 if Kworstid==Sworstid && Sworstid==Gworstid && Gworstid==Tworstid
                     eFworst_penalty=HM(5).Fit_penalty(Kworstid,:);
                     eidworst=Kworstid;
                 else
                     [eFworst_penalty, eidworst] = min(sum(HM(5).Fit_penalty,2));
                 end 
  
                if (Fit_penaltyNew(1) >= HM(5).Fit_penalty(eidworst,1)) && (Fit_penaltyNew(2) >= HM(5).Fit_penalty(eidworst,2)) && (Fit_penaltyNew(3) >= HM(5).Fit_penalty(eidworst,3)) && (Fit_penaltyNew(4) >= HM(5).Fit_penalty(eidworst,4)) %%%选取准确率较高的
                      HM(5).Fit(eidworst,:) = efnew;
                      HM(5).X(eidworst,:) = Xnew; 
                      HM(5).Fit_penalty(eidworst,:) =  Fit_penaltyNew;
                end

         end
    end
   
     
    end
end
toc

%%%% MCHS in the verifying stage 
HMS4=EliteNum+EliteNum2;
HM(5).X=[HM(4).X;HM(5).X];
HM(5).Fit=[HM(4).Fit;HM(5).Fit];
HM(5).Fit_penalty=[HM(4).Fit_penalty;HM(5).Fit_penalty];
for i=1:HMS4
    if min(HM(4).Fit(i,:))>0.95  && (isempty(HM(6).X) || ~ismember(HM(5).X(i,:),HM(6).X,'rows'))
         HM(6).X = [HM(6).X; HM(5).X(i,:)];
         HM(6).Fit = [HM(6).Fit; HM(5).Fit(i,:)];
         HM(6).Fit_penalty = [HM(6).Fit_penalty; HM(5).Fit_penalty(i,:)];
         OutEliteNum = OutEliteNum + 1;
     end
end
