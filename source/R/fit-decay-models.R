
# clear the environment
# load functions and data
# declare path to save model files
rm(list = ls())

path <- "source/final_code/luminex_recommendation/generated_rds/decay_model_results/"

source("source/final_code/luminex_recommendation/R/packages.R")
source("source/final_code/shared/utils.R")
source("source/final_code/luminex_recommendation/R/load-data.R")


markers <- unique(analysisData$new_long_data$full_test)

########## Univariate models (no covariates) ------------

#fit all exponential and biphasic models
# 1) select exponential model
# 2) fit model
# 3) save model fit
# 4) add summary statistic to a dataframe
# 5) compute LOO for fit and store somehow


for(m in markers){ 
        
        cat(m, "\n")
        
        
        compare_fit<- compare_decay_delay_form(delays=5,marker=m,iterations=5000,
                                               data=filter(analysisData$new_wide_data,
                                                           status=="Case"
                                               ),
                                               cens_df = filter(analysisData$new_wide_censor,
                                                                status=="Case"
                                               )
        )

        # #save fit objects in folder
        write_rds(compare_fit,
                  paste0(path,"fits-compare/compare-fit-",m,".rds"))
        
}


highRhat <- summary_df %>% 
        filter(str_detect(parameter,"mu_|halflife")) %>% 
        group_by(marker,model) %>%
        filter(Rhat==max(Rhat)) %>%
        ungroup() %>%
        filter(Rhat>=1.05)%>% 
        distinct(marker,model,Rhat)

for(m in highRhat$marker){ 
        
        cat(m, "\n")
        
        
        compare_fit<- compare_decay_delay_form(delays=5,marker=m,iterations=5000,
                                               data=filter(analysisData$new_wide_data,
                                                           status=="Case"
                                               ),
                                               cens_df = filter(analysisData$new_wide_censor,
                                                                status=="Case"
                                               )
        )
        
        # #save fit objects in folder
        write_rds(compare_fit,
                  paste0(path,"fits-compare/compare-fit-",m,".rds"))
        
        
}


#creat summary df
summary_df <-data.frame()
loo_df <- data.frame()

for(m in markers){ 
        
        cat(m, "\n")
        
        compare_fit<-  read_rds(paste0(path,"fits-compare/compare-fit-",m,".rds"))
        
        #save summaries in a dataframedataframe
        summary_df <- bind_rows(
                summary_df,
                compare_fit$fits$exponential_5$summary %>% mutate(model="exponential"),
                compare_fit$fits$biphasic_5$summary%>% mutate(model="biphasic"))

        # save loo comparisons in dataframe
        loo_df <- bind_rows(
                loo_df,
                compare_fit$loo %>% as.data.frame() %>%
                        mutate(marker=m,
                               model=str_extract(rownames(.),"exponential|biphasic"))
        )

        #save summary files for simple univariate models
        write_rds(summary_df,paste0(path,"fits-compare/summary_df.rds"))
        write_rds(loo_df, paste0(path,"fits-compare/loo_df.rds"))
        
}


############### Covariate Models------------------

cov <-c("Ogawa", "Under 10","Female","O Blood")
# markers <- markers[str_detect(markers,"NetMFI")]

exp_model_cov <- select_univariate_model("exponential",covariates = cov[1])


for(m in markers){
  for (c in cov){
          
  cat("marker:",m," covariate: ",c,"\n")
          
  demoCov_fit <-analysisData$new_wide_data %>% filter(status == "Case") %>%
                        stanfit_decay_univariate(
                                cens_df=analysisData$new_wide_censor %>%
                                        filter(status == "Case"),
                                marker=m,
                                delay=5,
                                covariates=c,
                                iterations=5000,
                                model= exp_model_cov,
                                chains=4,
                                warmup = 1000
                                 )

   write_rds(demoCov_fit, paste0(path, "covariates/cov-fit-",m,"-",c,".rds"))

}
}

summary_cov_df <- data.frame()
for(m in markers){
        for (c in cov){
                
                cat("marker:",m," covariate: ",c,"\n")
                
                demoCov_fit <- read_rds(paste0(path, "covariates/cov-fit-",m,"-",c,".rds"))
                
                
                summary_cov_df <-bind_rows(summary_cov_df,
                                           demoCov_fit$summary %>%
                                             mutate(covariate=c)
                                           )

                write_rds(summary_cov_df,
                          paste0(path, "covariates/summary_cov_df.rds"))
                
                
        }
}

#rerun those with an Rhat over 1.05

filter(summary_cov_df, str_detect(parameter,"mu_|halflife|beta")) %>%
        distinct(marker,covariate)  %>% nrow()

highRhat_cov <- filter(summary_cov_df, str_detect(parameter,"mu_|halflife|beta")) %>% 
        group_by(marker,covariate) %>%
        filter(Rhat==max(Rhat)) %>%
        ungroup() %>%
        filter(Rhat>=1.01) %>% distinct(marker,covariate,Rhat)
        
 for(r in 1:nrow(highRhat_cov)){

                 m <-highRhat_cov$marker[r]
                 c <- highRhat_cov$covariate[r]

                 cat("marker:",m," covariate: ",c,"\n")
                 
                 demoCov_fit <-analysisData$new_wide_data %>% filter(status == "Case") %>%
                         stanfit_decay_univariate(
                                 cens_df=analysisData$new_wide_censor %>%
                                         filter(status == "Case"),
                                 marker=m,
                                 delay=5,
                                 covariates=c,
                                 iterations=5000,
                                 model= exp_model_cov,
                                 chains=4,
                                 warmup = 1000
                         )
                 
                 write_rds(demoCov_fit, paste0(path, "covariates/cov-fit-",m,"-",c,".rds"))
                 

 }
 
 
        # 
# 
# summary_cov_df <- filter(summary_cov_df,covariate!="Female")
# 
demoCov_fit <- read_rds(paste0(path, "covariates/cov-fit-","NetMFI_IgG_VCC","-","Female",".rds"))
# 
demoCov_fit$fit %>% shinystan::launch_shinystan()


cp <- read_rds(
          paste0(path,"fits-compare/compare-fit-","NetMFI_IgA_Flu",".rds")) 
 
cp$fits$exponential_5$fit %>% shinystan::launch_shinystan()
