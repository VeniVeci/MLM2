% Result_demostration using boxplot 

RecRateSsErr_kpca=load('RecSSErr_kpca_nptrain100.txt');
RecRateSsErr_svd=load('RecSSErr_svd_nptrain100.txt');

new_dim_start=2;    
new_dim_step =2;
new_dim_end  =10;

figure(1)
boxplot(RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));

figure(2)
boxplot(RecRateSsErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of LPCA-GPR'));