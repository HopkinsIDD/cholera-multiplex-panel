##### Section 0: setup 
rm(list = ls())

path <- "source/final_code/luminex_recommendation/generated_rds/random_forest_results/"

source("source/final_code/luminex_recommendation/R/packages.R")
source("source/final_code/shared/utils.R")
source("source/final_code/luminex_recommendation/R/load-data.R")

#trees
tree_param <- 1000

#set up the four windows of highest interest
times <- c(45,120,200,300)


###### Conduct LOOCV to get serostatus----------


form_NetMFI <- analysisData$new_long_data %>%
        filter(str_detect(full_test,"NetMFI")) %>%
        filter(!(antigen %in% c("CTHT","LTh","LTB","Flu","O139BSA"))) %>%
        distinct(full_test) %>% pull(full_test) %>%
        c("age","sex","blood")

form_RAU <- analysisData$new_long_data %>%
        filter(str_detect(full_test,"RAU")) %>%
        filter(!(antigen %in% c("CTHT","LTh","LTB","Flu","O139BSA"))) %>%
        distinct(full_test) %>% pull(full_test) %>%
        c("age","sex","blood")

traditional_predictors <- casecon_analysis %>% 
        filter(test_type == "Vibriocidal"|
                       (test_type=="ELISA" &
                                isotype %in% c("IgG","IgA"))
        ) %>%
        distinct(full_test) %>%
        pull(full_test) %>%
        c("age","sex","blood")



var_list <- list(
                 "Traditional (with non-Immuno)"= traditional_predictors,
                 "Traditional"= traditional_predictors[!str_detect(traditional_predictors,"age|sex|blood")],
                 
                 "All NetMFI (with non-Immuno)" = form_NetMFI,
                 "All NetMFI IgG (with non-Immuno)"= form_NetMFI[str_detect(form_NetMFI,"IgG|age|sex|blood")],
                 
                 "All NetMFI" = form_NetMFI[!str_detect(form_NetMFI,"age|sex|blood")],
                 "All NetMFI IgG"= form_NetMFI[str_detect(form_NetMFI,"IgG") & !str_detect(form_NetMFI,"age|sex|blood")]
)

var_list[["Reduced Panel NetMFI (with non-Immuno)"]]<- var_list$`All NetMFI IgG (with non-Immuno)`[str_detect(var_list$`All NetMFI IgG (with non-Immuno)`,"CtxB|OSP|age|sex|blood")]
var_list[["Reduced Panel NetMFI"]]<- var_list$`All NetMFI IgG`[str_detect(var_list$`All NetMFI IgG`,"CtxB|OSP")]

var_list[["All RAU (with non-Immuno)"]] <- form_RAU
var_list[["All RAU IgG (with non-Immuno)"]] <-  form_RAU[str_detect(form_RAU,"IgG|age|sex|blood")]

var_list[["All RAU"]] <-  form_RAU[!str_detect(form_RAU,"age|sex|blood")]
var_list[["All RAU IgG"]]  <- form_RAU[str_detect(form_RAU,"IgG") & !str_detect(form_RAU,"age|sex|blood")]

var_list[["Reduced Panel RAU (with non-Immuno)"]]<- var_list$`All RAU IgG (with non-Immuno)`[str_detect(var_list$`All RAU IgG (with non-Immuno)`,"CtxB|OSP|age|sex|blood")]
var_list[["Reduced Panel RAU"]]<- var_list$`All RAU IgG`[str_detect(var_list$`All RAU IgG`,"CtxB|OSP")]


loocv_list_weighted <- list()
raw_df <- data.frame()


for(w in c(45,120,200,300)){
        # for(v in 1:length(var_list)){
                for(v in 9:14){
                        
                
                tmp_loocv <-loocv_ranger(wide_analysis_nomissing,
                                     end_window=w,
                                     variables=var_list[[v]],
                                     weighted=TRUE,
                                     num.trees=tree_param
                        )

                tmp_loocv$cv_df <- tmp_loocv$cv_df %>%
                        mutate(variables=names(var_list)[[v]]) %>%
                        mutate(end_window=w) %>%
                        getWeight(end_window=w)

                loocv_list_weighted[[as.character(w)]][[names(var_list)[v]]] <- tmp_loocv
                raw_df <- bind_rows(raw_df,tmp_loocv$cv_df)
                
                }
        write_rds(loocv_list_weighted,
                  paste0(path,"tvsens-spec/loocv_list_weighted.rds")
        )
        
        write_rds(raw_df,
                  paste0(path,"tvsens-spec/raw_df.rds")
        )
        
}


# raw_df <- raw_df %>% mutate(variables=factor(variables,levels=names(var_list)))



##### Estimate time-varying sensitivity  and specificity-------------


tvsens_fit <- list()
tvsens_df <- data.frame()
spec_fit <- list()
spec_df <- data.frame()

options(mc.cores = 1)


for(w in c(45,120,200,300)){
        for(v in unique(raw_df$variables)[9:14]){
                for(s in c("youden_seropos","spec90_seropos",
                           "spec95_seropos","spec99_seropos")){
                
                cat("Window: ",w, " Variables: ",v," Seropositivity: ",s,"\n")  
                
                # limit to particular window and variables
                reg_data <- raw_df %>%
                        filter(end_window==w) %>%
                        filter(variables==v)
                
                # estimate time-varying sensitivity
                tv_sens_tmp <- fit_tvsens(data = reg_data,
                                       end_window=w,
                                       seropos_type= s
                )
                
                tvsens_fit[[as.character(w)]][[v]][[s]] <- tv_sens_tmp
                tvsens_df <-bind_rows(tvsens_df,
                                      tv_sens_tmp$summary_df %>%
                                              mutate(variables=v)%>%
                                              mutate(seropos_type=s)
                                              )
                
                # estimate specificity
                spec_tmp <- fit_spec(data = reg_data,
                                          end_window=w,
                                          seropos_type= s
                )
                
                spec_fit[[as.character(w)]][[v]][[s]] <- spec_tmp
                spec_df <-bind_rows(spec_df,
                                      spec_tmp$summary_df %>% 
                                            mutate(variables=v) %>%
                                            mutate(seropos_type=s)
                                        )
                
                
                write_rds(tvsens_fit,
                          paste0(path,"tvsens-spec/tvsens_fit.rds"))
                write_rds(tvsens_df,
                          paste0(path,"tvsens-spec/tvsens_df.rds"))
                write_rds(spec_fit,
                          paste0(path,"tvsens-spec/spec_fit.rds"))
                write_rds(spec_df,
                          paste0(path,"tvsens-spec/spec_df.rds"))
                
                
                
                
}}}




# ###### Estimate average sensitivity------



tvsens_fit <- read_rds(paste0(path,"tvsens-spec/tvsens_fit.rds"))

avg_sens_df <- data.frame()

for(i in names(tvsens_fit)){
        for(j in names(tvsens_fit[[i]])){
                for(k in names(tvsens_fit[[i]][[j]])){
                        
                        times <- seq(10,as.numeric(i),7)
                        
                        tmp <- tvsens_fit[[i]][[j]][[k]][["out_df"]] %>%
                                        filter(time %in% times) %>%
                                        group_by(.draw) %>%
                                        #make approximation of integration under curve
                                        summarize(sens=sum(theta)/(max(times)-10+7)*7) %>%
                                        ungroup() %>%
                                        median_qi() %>%
                                        mutate(end_window=as.numeric(i)) %>%
                                        mutate(variables=j) %>%
                                        mutate(seropos_type=k)
                        
                        avg_sens_df <- bind_rows(avg_sens_df,tmp)
                                        
                                
                                
        
                }}}


write_rds(avg_sens_df,
          paste0(path,"tvsens-spec/avgsens_df.rds"))

###### Estimate cutoffs----------
 
# #small function to get positivity rate with hierarchichal model
# getPR <- function(this_df){
#         
#         data1 <- this_df %>% filter(inf_120==1)
#         data0 <- this_df %>% filter(inf_120==0)
#         
#         PR1 <- mean(data1$positive)
#         PR0 <- mean(data0$positive)
#         
#         #sensitivity model
#         if(!mean(data1$positive) %in% c(1,0)){ #dont fit the model if all pos/neg
#                 fit1  <- lme4::glmer(data=data1, positive ~ (1|id),
#                                      family = binomial()
#                 ) %>% summary()
#                 
#                 
#                 
#                 PR1 = fit1$coefficients[,1] %>%locfit::expit()
#                 
#         }
#         
#         #specificity model
#         if(!mean(data0$positive) %in% c(1,0)){#dont fit the model if all pos/neg
#                 fit0  <- lme4::glmer(data=data0, positive ~ (1|id),
#                                      family = binomial()
#                 ) %>% summary()
#                 
#                 PR0= fit0$coefficients[,1] %>%locfit::expit()
#         }
#         
#         data.frame(inf_120=c(1,0),
#                    PR=c(PR1,
#                         PR0
#                    )
#         )
#         
#         
# }
# 
# 
# 
# fit_heatmap_df <- data.frame()
# 
# for(i in seq(-5,-2,0.1)){
#         for(j in seq(-5,-2,0.1)){
#                 
#                 
#                 tmp <- wide_analysis %>%addInfectionWindow(120) %>% 
#                         mutate(positive=ifelse((Luminex_IgG_CtxB>= i)  & (Luminex_IgG_OgawaOSPBSA>= j ),
#                                                1,0)) %>%
#                         getPR() %>%
#                         mutate(Luminex_IgG_OgawaOSPBSA=j,
#                                Luminex_IgG_CtxB=i
#                         )
#                 
#                 fit_heatmap_df <- bind_rows(fit_heatmap_df,tmp)
#         }
# }
# 
# 
# 
# fit_ssy_df <- fit_heatmap_df %>% mutate(param=case_when(inf_120 =="0" ~ "spec",
#                                                         inf_120 =="1" ~ "sens",
# )) %>% 
#         mutate(value=ifelse(param=="spec",1-PR,PR)) %>%
#         select(-PR,-inf_120)%>%
#         spread(param,value) %>%
#         mutate(youden=sens+spec-1)
# 
# 
# 
# write_rds(fit_ssy_df,
#           paste0(path,"tvsens-spec/cutoff_fit_ssy_df.rds"))



