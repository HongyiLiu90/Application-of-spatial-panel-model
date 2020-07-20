

% Hongyi Liu
% 03/14/2014
% James Lesage was so kind to make his code available to estimate SLM, SDM, and SEM with direct, indirect and total effects.
% This matlab code is based on James Lesage (2009) and his toolbox.

% We construct a novel economic-geographical spatial weighted matrix.


load p2.txt % input data
load w2.txt % input the economic-geographical spatial weighted matrix
T=15; % number of time periods
N=30; % number of regions
% row-normalize W
W=normw(w2); 
y=p2(:,[3]); 
x=p2(:,[8:17]); 
for t=1:T
t1=(t-1)*N+1;t2=t*N;
wx(t1:t2,:)=W*x(t1:t2,:);
end
xconstant=ones(N*T,1);
[nobs K]=size(x);
% ----------------------------------------------------------------------------------------
% No fixed effects + spatially lagged dependent variable
info.lflag=0; % required for exact results
info.model=0;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
% New routines to calculate effects estimates
results=sar_panel_FE(y,[xconstant x],W,T,info); 
vnames=char('g','intercept','logLgdp','logsk','logsh','logn','logpopden','lnwage','bus','road','toil','llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=0;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sar(results,vnames,W);
% ----------------------------------------------------------------------------------------
% No fixed effects + spatially lagged dependent variable + spatially
% independent variables
info.lflag=0; % required for exact results
info.model=0;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
% New routines to calculate effects estimates
results=sar_panel_FE(y,[xconstant x wx],W,T,info); 
vnames=char('g','intercept','logLgdp','logsk','logsh','logn','logpopden',...
  'lnwage','bus','road','toil','llogcla','W*logLgdp','W*logsk','W*logsh','W*logn',...
  'W*logpopden','W*lnwage','W*bus','W*road','W*toil','W*llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=1;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sdm(results,vnames,W);
% ----------------------------------------------------------------------------------------
% Spatial fixed effects + spatially lagged dependent variable
info.lflag=0; % required for exact results
info.model=1;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
% New routines to calculate effects estimates
results=sar_panel_FE(y,x,W,T,info); 
vnames=char('g','logLgdp','logsk','logsh','logn','logpopden','lnwage','bus','road','toil','llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=0;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sar(results,vnames,W);
% ----------------------------------------------------------------------------------------
% Spatial fixed effects + spatially lagged dependent variable + spatially
% independent variables
info.lflag=0; % required for exact results
info.model=1;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
% New routines to calculate effects estimates
results=sar_panel_FE(y,[x wx],W,T,info); 
vnames=char('g','logLgdp','logsk','logsh','logn','logpopden',...
  'lnwage','bus','road','toil','llogcla','W*logLgdp','W*logsk','W*logsh','W*logn',...
  'W*logpopden','W*lnwage','W*bus','W*road','W*toil','W*llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=1;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sdm(results,vnames,W);
% ----------------------------------------------------------------------------------------
% Time period fixed effects + spatially lagged dependent variable
info.lflag=0; % required for exact results
info.model=2;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
% New routines to calculate effects estimates
results=sar_panel_FE(y,x,W,T,info); 
vnames=char('g','logLgdp','logsk','logsh','logn','logpopden','lnwage','bus','road','toil','llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=0;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sar(results,vnames,W);
% ----------------------------------------------------------------------------------------
% Time period fixed effects + spatially lagged dependent variable + spatially
% independent variables
info.lflag=0; % required for exact results
info.model=2;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
% New routines to calculate effects estimates
results=sar_panel_FE(y,[x wx],W,T,info); 
vnames=char('g','logLgdp','logsk','logsh','logn','logpopden',...
  'lnwage','bus','road','toil','llogcla','W*logLgdp','W*logsk','W*logsh','W*logn',...
  'W*logpopden','W*lnwage','W*bus','W*road','W*toil','W*llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=1;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sdm(results,vnames,W);
% ----------------------------------------------------------------------------------------
% Spatial and time period fixed effects + spatially lagged dependent variable
info.lflag=0; % required for exact results
info.model=3;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
% New routines to calculate effects estimates
results=sar_panel_FE(y,x,W,T,info); 
vnames=char('g','logLgdp','logsk','logsh','logn','logpopden','lnwage','bus','road','toil','llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=0;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sar(results,vnames,W);
% ----------------------------------------------------------------------------------------
% Spatial and time period fixed effects + spatially lagged dependent variable + spatially
% independent variables
% No bias correction
info.bc=0;
info.lflag=0; % required for exact results
info.model=3;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
% New routines to calculate effects estimates
results=sar_panel_FE(y,[x wx],W,T,info); 
vnames=char('g','logLgdp','logsk','logsh','logn','logpopden',...
  'lnwage','bus','road','toil','llogcla','W*logLgdp','W*logsk','W*logsh','W*logn',...
  'W*logpopden','W*lnwage','W*bus','W*road','W*toil','W*llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=1;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sdm(results,vnames,W);
% Wald test for spatial Durbin model against spatial lag model
btemp=results.parm;
varcov=results.cov;
Rafg=zeros(K,2*K+2);
for k=1:K
Rafg(k,K+k)=1; % R(1,3)=0 and R(2,4)=0;
end
Wald_spatial_lag=(Rafg*btemp)'*inv(Rafg*varcov*Rafg')*Rafg*btemp
prob_spatial_lag=1-chis_cdf (Wald_spatial_lag, K) % probability greater than 0.05 points to insignificance
% LR test spatial Durbin model against spatial lag model (requires
% estimation results of the spatial lag model to be available)
resultssar=sar_panel_FE(y,x,W,T,info); 
LR_spatial_lag=-2*(resultssar.lik-results.lik)
prob_spatial_lag=1-chis_cdf (LR_spatial_lag,K) % probability greater than 0.05 points to insignificance
% Wald test spatial Durbin model against spatial error model
R=zeros(K,1);
for k=1:K
R(k)=btemp(2*K+1)*btemp(k)+btemp(K+k); % k changed in 1, 7/12/2010
% R(1)=btemp(5)*btemp(1)+btemp(3);
% R(2)=btemp(5)*btemp(2)+btemp(4);
end
Rafg=zeros(K,2*K+2);
for k=1:K
Rafg(k,k) =btemp(2*K+1); % k changed in 1, 7/12/2010
Rafg(k,K+k) =1;
Rafg(k,2*K+1)=btemp(k);
% Rafg(1,1)=btemp(5);Rafg(1,3)=1;Rafg(1,5)=btemp(1);
% Rafg(2,2)=btemp(5);Rafg(2,4)=1;Rafg(2,5)=btemp(2);
end 
Wald_spatial_error=R'*inv(Rafg*varcov*Rafg')*R
prob_spatial_error=1-chis_cdf (Wald_spatial_error,K) % probability greater than 0.05 points to insignificance
% LR test spatial Durbin model against spatial error model (requires
% estimation results of the spatial error model to be available
resultssem=sem_panel_FE(y,x,W,T,info); 
LR_spatial_error=-2*(resultssem.lik-results.lik)
prob_spatial_error=1-chis_cdf (LR_spatial_error,K) % probability greater than 0.05 points to insignificance
% ----------------------------------------------------------------------------------------
% Spatial and time period fixed effects + spatially lagged dependent variable + spatially
% independent variables
info.lflag=0; % required for exact results
info.model=3;
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
info.bc=1;
% New routines to calculate effects estimates
results=sar_panel_FE(y,[x wx],W,T,info); 
vnames=char('g','logLgdp','logsk','logsh','logn','logpopden',...
  'lnwage','bus','road','toil','llogcla','W*logLgdp','W*logsk','W*logsh','W*logn',...
  'W*logpopden','W*lnwage','W*bus','W*road','W*toil','W*llogcla');
% Print out coefficient estimates
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=1;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sdm(results,vnames,W);
% Wald test for spatial lag model
btemp=results.parm;
varcov=results.cov;
Rafg=zeros(K,2*K+2);
for k=1:K
Rafg(k,K+k)=1; % R(1,3)=0 and R(2,4)=0;
end
Wald_spatial_lag=(Rafg*btemp)'*inv(Rafg*varcov*Rafg')*Rafg*btemp
prob_spatial_lag= 1-chis_cdf (Wald_spatial_lag, K) % probability greater than 0.05 points to insignificance
% LR test spatial Durbin model against spatial lag model (requires
% estimation results of the spatial lag model to be available)
resultssar=sar_panel_FE(y,x,W,T,info); 
LR_spatial_lag=-2*(resultssar.lik-results.lik)
prob_spatial_lag=1-chis_cdf (LR_spatial_lag,K) % probability greater than 0.05 points to insignificance
% Wald test for spatial error model
R=zeros(K,1);
for k=1:K
R(k)=btemp(2*K+1)*btemp(k)+btemp(K+k); % k changed in 1, 7/12/2010
% R(1)=btemp(5)*btemp(1)+btemp(3);
% R(2)=btemp(5)*btemp(2)+btemp(4);
end
Rafg=zeros(K,2*K+2);
for k=1:K
Rafg(k,k) =btemp(2*K+1); % k changed in 1, 7/12/2010
Rafg(k,K+k) =1;
Rafg(k,2*K+1)=btemp(k);
% Rafg(1,1)=btemp(5);Rafg(1,3)=1;Rafg(1,5)=btemp(1);
% Rafg(2,2)=btemp(5);Rafg(2,4)=1;Rafg(2,5)=btemp(2);
end 
Wald_spatial_error=R'*inv(Rafg*varcov*Rafg')*R
prob_spatial_error= 1-chis_cdf (Wald_spatial_error,K) % probability greater than 0.05 points to insignificance
% LR test spatial Durbin model against spatial error model (requires
% estimation results of the spatial error model to be available
resultssem=sem_panel_FE(y,x,W,T,info); 
LR_spatial_error=-2*(resultssem.lik-results.lik)
prob_spatial_error=1-chis_cdf (LR_spatial_error,K) % probability greater than 0.05 points to insignificance
% needed for Hausman test later on
logliklag=results.lik; 
blagfe=results.parm(1:end-1);
covblagfe=results.cov(1:end-1,1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% random effects estimator by ML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Spatial random effects and time period fixed effects + spatially lagged dependent variable + spatially
% independent variables
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x wx],N,T,2); % 2=time dummies
info.model=1;
results=sar_panel_RE(ywith,xwith,W,T,info); 
prt_sp(results,vnames,1);
% Print out effects estimates
spat_model=1;
direct_indirect_effects_estimates(results,W,spat_model);
panel_effects_sdm(results,vnames,W);
% Wald test for spatial lag model
btemp=results.parm(1:2*K+2);
varcov=results.cov(1:2*K+2,1:2*K+2);
Rafg=zeros(K,2*K+2);
for k=1:K
Rafg(k,K+k)=1; % R(1,3)=0 and R(2,4)=0;
end
Wald_spatial_lag=(Rafg*btemp)'*inv(Rafg*varcov*Rafg')*Rafg*btemp
prob_spatial_lag= 1-chis_cdf (Wald_spatial_lag, K) % probability greater than 0.05 points to insignificance
% Wald test for spatial error model
R=zeros(K,1);
for k=1:K
R(k)=btemp(2*K+1)*btemp(k)+btemp(K+k); % k changed in 1, 7/12/2010
% R(1)=btemp(5)*btemp(1)+btemp(3);
% R(2)=btemp(5)*btemp(2)+btemp(4);
end
Rafg=zeros(K,2*K+2);
for k=1:K
Rafg(k,k) =btemp(2*K+1); % k changed in 1, 7/12/2010
Rafg(k,K+k) =1;
Rafg(k,2*K+1)=btemp(k);
% Rafg(1,1)=btemp(5);Rafg(1,3)=1;Rafg(1,5)=btemp(1);
% Rafg(2,2)=btemp(5);Rafg(2,4)=1;Rafg(2,5)=btemp(2);
end 
Wald_spatial_error=R'*inv(Rafg*varcov*Rafg')*R
prob_spatial_error= 1-chis_cdf (Wald_spatial_error,K) % probability greater than 0.05 points to insignificance
% needed for Hausman test later on
logliklagre=results.lik;
blagre=results.parm(1:end-2);
covblagre=results.cov(1:end-2,1:end-2);
% ----------------------------------------------------------------------------------------
% Hausman test FE versus RE 
hausman=(blagfe-blagre)'*inv(covblagre-covblagfe)*(blagfe-blagre);
dof=length(blagfe);
probability=1-chis_prb(abs(hausman),dof); 
% Note: probability < 0.025 implies rejection of random effects model in favor of fixed effects model
% Use 0.025, since it is a one-sided test
fprintf(1,'Hausman test-statistic, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',abs(hausman),dof,probability);

