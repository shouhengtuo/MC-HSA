function Fit = multiObjFitness(X,Y)

   Fit(1,1) = fitnessKNN(X,Y) ;
     
   Fit(1,2) = fitnessSVM(X,Y) ;
  
   Fit(1,3) = fitnessGLM(double(X),double(Y)) ; 
   %Fit(1,3) = fitnessMLP(X,Y) ; 
   Fit(1,4) = fitnessTreeBagger(X,Y) ;
  
end