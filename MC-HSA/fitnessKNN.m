function accuracy = fitnessKNN(X,Y) 

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
                knnModel=fitcknn(Xtrain,Ytrain,'NumNeighbors',12,'Standardize',1);
            [y1,~] = predict(knnModel,Xtest);
            [s,~]=size(Ytest);
            c = sum( y1~= Ytest)/s;
            acc(i)=1-c;
        end
        accuracy = mean(acc);
     
end