function accuracy = fitnessSVM(X,Y)  
 
  SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','RBF',...
     'KernelScale','auto');
%    t = templateSVM('Standardize',1,'KernelFunction','RBF');
%    SVMModel = fitcecoc(X,Y,'Learners',t); 
  CVSVMModel = crossval(SVMModel);
  classLoss = kfoldLoss(CVSVMModel);
 
  accuracy = 1 - classLoss;
%   Y1 = predict(SVMModel,X);
%    accuracy = 1 - sum(abs(Y1-Y))/length(Y);
end