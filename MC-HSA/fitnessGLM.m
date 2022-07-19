function accuracy = fitnessGLM(X,Y) 
    warning off
     %% Ωª≤Ê—È÷§
        L = length(Y);
        trainSize = fix(L*0.8);
        for i = 1:10
           IDX =  randperm(L);
           trainIdx = IDX(1:trainSize);
           testIdx = IDX(trainSize+1 : end);
           Xtrain = X(trainIdx,:);
           Ytrain = Y(trainIdx);
           Xtest = X(testIdx,:);
           Ytest = Y(testIdx);
             Model = fitglm(Xtrain,Ytrain);
            [yp,~] = predict(Model,Xtest);
            y = round(yp);
            y = max(min(y,1),0);
            [s,~]=size(Ytest);
            c = sum( y~= Ytest)/s;
            acc(i)=1-c;
        end
        accuracy = mean(acc);
     
end