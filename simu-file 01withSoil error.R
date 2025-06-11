################################################################################
######   Simulate True and Error Data
################################################################################
library(mvtnorm)
library(invgamma)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(xtable)
set.seed(200)


# Set seed for reproducibility
set.seed(123)

# Define parameters
n <- 100  # Total number of samples

# Generate data for soil characteristics and land size
land_size <- runif(n, 1, 10)  # Assuming a uniform distribution for simplicity
soil_characteristics <- runif(n, land_size/15, 1)  # make it depend on land size 


# simulate model matrix
#sigma_u_L =0.1#0.05 0.10 0.25 0.50
sigma_u_S =0.25#0.05 0.10 0.25 0.50 
x<-cbind(1, land_size, soil_characteristics)

# true beta coefficients
true_beta_coef<-c(100, 50,5)

# true sigma for yield
true_sigma <- 50
I<-diag(1,n,n) # n by n identity covariance matrix 

# simulate yield in tons  
y<-t(rmvnorm(1, x%*%true_beta_coef, true_sigma*I))
# simulate many outcomes for later asymptotic evaluations
y_list<-replicate(1000, t(rmvnorm(1, x%*%true_beta_coef, true_sigma*I)),simplify = FALSE)

#error covariates
x<-cbind(1, land_size, soil_characteristics+rnorm(n,0,sigma_u_S))


############################
#   gibbs samlping 
############################
blocked_gibbs<-function(y, x, iter, burnin, thin){
  # initialize gibbs sampler scheme
  xxt_inv<-solve(t(x)%*%x) # predictor
  sigma<-numeric(iter) # def sigma vector
  beta_coef<-matrix(nrow=iter, ncol = 3) # vector for betas
  sigma[1]<-5 # initial sigma value to start the chain
  
  # sigma hyperparameters
  a<- 1/2
  b<-1.0E+2
  
  # gibbs sampling iterations
  for(i in 2:iter ){
    beta_coef[i,]<-rmvnorm(n = 1, 
                           mean = ((xxt_inv%*%t(x))%*%y), 
                           sigma = sigma[i-1]*xxt_inv )
    
    sigma[i]<-rinvgamma(n = 1, 
                        shape = (n/2 + a), 
                        rate = .5*( t((y - x%*%t(t(beta_coef[i,])) ))%*%(y - x%*%t(t(beta_coef[i,])) ) ) + b)
    
  }
  
  # apply burnin and trimming  
  keep_draws<-seq(burnin,iter,thin)
  sigma<-sigma[keep_draws]
  beta_coef<-beta_coef[keep_draws,]
  
  # format and output
  joint_post<-data.frame(beta_coef=beta_coef,sigma=sigma)
  colnames(joint_post)[1:(ncol(x))]<-paste0('beta',0:(ncol(x)-1) )
  
  joint_post_long<-gather(joint_post,keep_draws) %>%
    rename(param=keep_draws, draw=value) %>%
    mutate(iter=rep(keep_draws,ncol(joint_post)))
  
  return(joint_post_long)
}

# run Gibbs sampler with specified parameters
post_dist<-blocked_gibbs(y = y, x = x, iter = 3000, burnin = 1000, thin = 1)

###### Posterior distributions  summaries  ######
# calculate posterior summary statistics 
post_sum_stats<-post_dist %>%
  group_by(param) %>%
  summarise(median=median(draw),
            lwr=quantile(draw,.025),
            upr=quantile(draw,.975)) %>%
  mutate(true_vals=c(true_beta_coef, true_sigma))

# merge on summary statistics
post_dist <- post_dist %>%
  left_join(post_sum_stats, by='param')

# plot MCMC Chains
plot1<-ggplot(post_dist,aes(x=iter,y=draw)) +
  geom_line() +
  geom_hline(aes(yintercept=true_vals, col='red'), show.legend=FALSE)+
  facet_grid(param ~ .,scale='free_y',switch = 'y') +
  theme_bw() + 
  xlab('Gibbs Sample Iteration') + ylab('MCMC Chains') + 
  ggtitle('Gibbs Sampler MCMC Chains by Parameter')

# plot Posterior Distributions
plot2<-ggplot(post_dist,aes(x=draw)) +
  geom_histogram(aes(x=draw),bins=50) +
  geom_vline(aes(xintercept = true_vals,col='red'), show.legend = FALSE) +
  facet_grid(. ~ param, scale='free_x',switch = 'y') +
  theme_bw() + 
  xlab('Posterior Distributions') + ylab('Count') + 
  ggtitle('Posterior Distributions of Parameters (true values in red)')

# Combine the two plots in two rows
combined_plots <- grid.arrange(plot1, plot2, ncol = 1)


# Display the combined plots
print(combined_plots)

# Autocorrelation plot
acf(post_dist$draw)

######  Assessing bias and coverage  ######
bayes_res<-lapply(y_list, blocked_gibbs, x=x, iter=5000, burnin=2500, thin=1)

calc_sumstats<-function(post_dist){
  post_sum_stats<-post_dist %>%
    group_by(param) %>%
    summarise(post_median=median(draw),
              lwr=quantile(draw,.025),
              upr=quantile(draw,.975)) %>%
    mutate(true_vals=c(true_beta_coef, true_sigma))
  return(post_sum_stats)
}


all_sum_stats<-lapply(bayes_res, calc_sumstats)
all_sum_stats_stack<-bind_rows(all_sum_stats) %>%
  arrange(param) %>%
  rename(est=post_median)

eval_sum <- all_sum_stats_stack %>%
  mutate(covered=ifelse(true_vals<upr & true_vals>lwr,1,0)) %>%
  group_by(param) %>%
  summarise(est_var=var(est),
            est_mean=mean(est),
            bias=mean(est-true_vals),
            true_val=mean(true_vals),
            coverage=mean(covered)) %>%
  mutate(perc_bias=(bias/true_val)*100)

# format table
eval_sum<-eval_sum[,c(1,5,3,2,4,7,6)]
names(eval_sum)<-c('Parameter','True Value','Estimator Mean',
                   'Estimator Variance','Bias','Percent Bias (of truth)',
                   'Coverage of 95% CI')
View(eval_sum)

print.xtable(xtable(eval_sum,
                    caption =paste0("Estimator Evaluation sample size = ",n, 
                                    ", sigma_u Land Size = ", sigma_u_L)), 
                                    
             caption.placement = 'top', 
             include.rownames = F)
              
  
   