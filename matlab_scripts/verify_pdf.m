dat=csvread('../Result/slope.csv');
IQR=quantile(dat,0.75)-quantile(dat,0.25);
thresh=quantile(dat,0.75)+3*IQR;
Idx=find(dat<thresh);
dat=dat(Idx);
min(dat)
max(dat)