
##### Section 0: setup 
rm(list = ls())

path <- "source/final_code/luminex_recommendation/generated_rds/random_forest_results/"

source("source/final_code/luminex_recommendation/R/packages.R")
source("source/final_code/shared/utils.R")
source("source/final_code/luminex_recommendation/R/load-data.R")


tree_param <- 1000
folds <- 10
#set up all the days when scanning across windows
days <- c(45,seq(50,600,by=10))
#set up the four windows of highest interest
times <- c(45,120,200,300)


#List the different combination of predictions

#all net MFI and non-immunovariables
pred_NetMFI_all <-  analysisData$new_long_data %>%
        filter(str_detect(full_test,"NetMFI")) %>%
        distinct(full_test) %>% pull(full_test) %>% 
        c("age","sex","blood")

#remove control variables
pred_NetMFI_limited <-  pred_NetMFI_all[!str_detect(pred_NetMFI_all,"CTHT|LTh|LTB|Flu|O139BSA")]
#remove non-immunovariables
pred_NetMFI_immuno <-  pred_NetMFI_limited[!str_detect(pred_NetMFI_limited,"age|sex|blood")]

#all net MFI and non-immunovariables
pred_RAU_all<-  analysisData$new_long_data %>%
        filter(str_detect(full_test,"RAU")) %>%
        distinct(full_test) %>% pull(full_test) %>% 
        c("age","sex","blood")

#remove control variables
pred_RAU_limited <-  pred_RAU_all[!str_detect(pred_RAU_all,"CTHT|LTh|LTB|Flu|O139BSA")]
#remove non-immunovariables
pred_RAU_immuno <- pred_RAU_limited[!str_detect(pred_RAU_limited,"age|sex|blood")]




##### Section 1: Estimate cvAUC across windows for MBA---------

#set up prediction list
pred_list_1 <- list( 
        `Net MFI All` = pred_NetMFI_all,
        `Net MFI Limited` = pred_NetMFI_limited,
        `Net MFI Immuno Only` = pred_NetMFI_immuno,
        
        `RAU All` = pred_RAU_all,
        `RAU Limited` = pred_RAU_limited,
        `RAU Immuno Only` = pred_RAU_immuno
)


#set up for loop
scan_cvAUC_df_1 <- data.frame()
scan_cvPred_df_1 <- data.frame()

for (W in c(TRUE,FALSE)){
for (pred in names(pred_list_1)){
        for (d in days){

                cat("Predictors: ",pred,"Weighted: ", W,"Window:",d,"\n")

                cv_tmp <- kfoldcv_ranger(analysisData$new_wide_data, k=folds,
                                         weighted=W,
                                         end_window=d,
                                         variables= pred_list_1[[pred]],
                                         num.trees=tree_param)

                scan_cvAUC_df_1 <- bind_rows(scan_cvAUC_df_1,
                                     cv_tmp$cvAUC %>%
                                             mutate(predictors=pred)
                                             
                                     )
                
                
                cv_tmp$cv_df["outcome"] <- cv_tmp$cv_df[paste0("inf_",d)]
                cv_tmp$cv_df["predictors"] <- pred
                
                
                #for the ROC curve
                scan_cvPred_df_1 <- bind_rows(scan_cvPred_df_1,
                                       cv_tmp$cv_df
                )



        }
}
}

#save these files
write_rds(scan_cvAUC_df_1,
          paste0(path,"cvAUC_importance/1_scan_cvAUC_df.rds")
          )

write_rds(scan_cvPred_df_1,
          paste0(path,"cvAUC_importance/1_scan_cvPred_df.rds")
)





##### Section 2: Estimate cvAUC across windows for MBA using SuperLearner---------

#set up prediction list
pred_list_2 <- list( 
        `Net MFI All` = pred_NetMFI_all,
        `RAU All` = pred_RAU_all
)

#set up for loop
superL_cvAUC_df_2 <- data.frame()

#load superlearner
library(SuperLearner)
huh <- c()

for(t in times){
        for (pred in names(pred_list_2)){
                

        cat("Window: ",t, " Variables: ",pred)

        CV_ranger <- analysisData$new_wide_data %>% kfoldcv_ranger(end_window = t,
                                                      variables = pred_list_2[[pred]],
                                                      weighted = FALSE,
                                                      k=folds,
                                                      num.trees=tree_param
        )

        CV_SuperLearner <- analysisData$new_wide_data %>%
                kfoldcv_SuperLearner(end_window = t,
                                     variables = pred_list_2[[pred]],
                                     k=folds)


        superL_cvAUC_df_2 <- bind_rows(
                superL_cvAUC_df_2,
                CV_ranger$cvAUC %>%
                        mutate(Model="Random Forest",
                               measure=pred
                               ),
                CV_SuperLearner$cvAUC%>%
                        mutate(Model="Ensemble",
                               mesure=pred
                               )


        )


}
}

write_rds(superL_cvAUC_df_2,
          paste0(path,"cvAUC_importance/2_superL_cvAUC_df.rds")
)
##### Section 3: Estimate importance across windows---------

importance_df <- data.frame()
for (W in c(TRUE,FALSE)){
for (pred in c("Net MFI Limited","RAU Limited")){
for(t in days){
        cat("Weighted:",W,"Predictors:",pred,"Window:",t,"\n")

        #add the outcome variable to data and get weightes
        fit_data <- analysisData$new_wide_data %>%
                addInfectionWindow(end_window=t) %>%
                        getWeight(end_window = t)
        
        made_weights <- rep(nrow(fit_data),1)
        
        if(W ==TRUE){
                made_weights <- fit_data %>% pull(weight)
        }

        outcome <- paste0("inf_",t)
        

        tmp_fit<- ranger::ranger(data=fit_data,
                                          formula=make_formula(outcome,pred_list_1[[pred]]),
                                          num.trees= 5000,
                                          case.weights = made_weights,
                                          replace = TRUE,
                                          importance = "permutation")

        tmp_df <- data.frame(importance=ranger::importance(tmp_fit),
                             end_window=t) %>%
                mutate(marker=rownames(.)) %>%
                mutate(weighted=W,
                       predictors=pred
                       )

        importance_df <- bind_rows(importance_df,tmp_df)

}
}}


write_rds(importance_df,
          paste0(path,"cvAUC_importance/3_importance_df.rds")
          )


##### Section 4: Estimate cvAUC across different assay types ---------

traditional_predictors <- casecon_analysis %>% 
        filter(test_type == "Vibriocidal"|
                       (test_type=="ELISA" &
                                isotype %in% c("IgG","IgA"))
        ) %>%
        distinct(full_test) %>%
        pull(full_test) 

vibonly_predictors <- casecon_analysis %>%
        filter(test_type=="Vibriocidal") %>%
        distinct(full_test) %>%
        pull(full_test) 

ELISAonly_predictors <- casecon_analysis %>% filter(test_type=="ELISA" &
                                                            isotype %in% c("IgG","IgA")) %>%
        distinct(full_test) %>%
        pull(full_test) 

# get the luminex predictors
key_antigens <- c("OgawaOSPBSA","InabaOSPBSA","CtxB","Sialidase","VCC","TcpA")


Luminex_predictors <- casecon_analysis %>% filter(test_type=="Luminex") %>%
        filter(antigen %in% key_antigens) %>%
        distinct(full_test) %>%
        pull(full_test) 


Luminex_IgG_predictors <- casecon_analysis %>% filter(test_type=="Luminex") %>%
        filter(isotype=="IgG") %>%
        filter(antigen %in% key_antigens) %>%
        distinct(full_test) %>%
        pull(full_test) 

Luminex_IgA_predictors <- casecon_analysis %>% filter(test_type=="Luminex") %>%
        filter(isotype=="IgA") %>%
        filter(antigen %in% key_antigens) %>%
        distinct(full_test) %>%
        pull(full_test) 

Luminex_IgM_predictors <- casecon_analysis %>% filter(test_type=="Luminex") %>%
        filter(isotype=="IgM") %>%
        filter(antigen %in% key_antigens) %>%
        distinct(full_test) %>%
        pull(full_test) 


# First show which isotype is most important
pred_list_4 <- list(
        `Vibriocidal &\nELISA Markers`=traditional_predictors,
        `Vibriocidal\nMarkers`=vibonly_predictors,
        `ELISA\nMarkers`=ELISAonly_predictors,
        
        `All NetMFI\nMarkers`=Luminex_predictors[str_detect(Luminex_predictors,"NetMFI")],
        `IgG NetMFI\nMarkers`=Luminex_IgG_predictors[str_detect(Luminex_IgG_predictors,"NetMFI")],
        `IgA NetMFI\nMarkers`=Luminex_IgA_predictors[str_detect(Luminex_IgA_predictors,"NetMFI")],
        `IgM NetMFI\nMarkers`=Luminex_IgM_predictors[str_detect(Luminex_IgM_predictors,"NetMFI")],
        
        `All RAU\nMarkers`=Luminex_predictors[str_detect(Luminex_predictors,"RAU")],
        `IgG RAU\nMarkers`=Luminex_IgG_predictors[str_detect(Luminex_IgG_predictors,"RAU")],
        `IgA RAU\nMarkers`=Luminex_IgA_predictors[str_detect(Luminex_IgA_predictors,"RAU")],
        `IgM RAU\nMarkers`=Luminex_IgM_predictors[str_detect(Luminex_IgM_predictors,"RAU")]
)



cv_marker_df <- data.frame()
marker_cvAUC_list <- list()

for (i in 1:length(pred_list_4)){
        for(j in times){
                
                cat(names(pred_list_4)[i]," window: ",j,"\n")
                
                tmp_obj <- kfoldcv_ranger(wide_analysis_nomissing,
                                       k=folds,
                                       weighted=TRUE,
                                       end_window=j,
                                       variables=pred_list_4[[i]],
                                       num.trees=tree_param)
                
                cv_marker_df <- bind_rows(cv_marker_df,
                                          tmp_obj$cvAUC %>%
                                                  mutate(
                                                          model=names(pred_list_4)[i])
                )
                
                
                
                marker_cvAUC_list[[names(pred_list_4)[i]]][[as.character(j)]] <-  tmp_obj 
                
                
                
        }
}


write_rds(pred_list_4,
          paste0(path,"cvAUC_importance/4_marker_set.rds")
)

write_rds(cv_marker_df,
          paste0(path,"cvAUC_importance/4_compare_cvAUC_markers_df.rds")
)

write_rds(marker_cvAUC_list,
          paste0(path,"cvAUC_importance/4_compare_cvAUC_markers_list.rds")
          )

##### Section 5: Estimate cvAUC across combinations in all isotypes---------

###order of importance
# NetMFI_IgG_predictors <-paste0("NetMFI_IgG_", c("CtxB",
#                                           "OgawaOSPBSA",
#                                           "InabaOSPBSA",
#                                           "TcpA","VCC","Sialidase"))
# 
# #add the outcome variable to data and get weightes
# fit_data <- wide_analysis_nomissing %>%
#         addInfectionWindow(end_window=200) %>%
#         getWeight(end_window = 200)
# made_weights <- fit_data %>% pull(weight)
# outcome <- paste0("inf_",200)
# 
# netmfi_igg_fit<- ranger::ranger(data=fit_data,
#                formula=make_formula(outcome,NetMFI_IgG_predictors),
#                num.trees= 5000,
#                case.weights = made_weights,
#                replace = TRUE,
#                importance = "permutation")
# 
# netmfi_igg_fit$variable.importance %>% sort()


cvAUC_combo_5 <- data.frame()
for(measure in c("NetMFI","RAU")){
        for(iso in c("IgG","IgM","IgA")){
                ###order of importance
                predictors <-paste(measure,iso, c("CtxB",
                                                           "OgawaOSPBSA",
                                                           "InabaOSPBSA",
                                                           "TcpA","Sialidase","VCC"),
                                   sep="_"
                                   )
                for (i in 1:length(predictors)){
                        for(j in times){
                                
                                cat(predictors[i]," window: ",j,"\n")
                                
                                
                                tmp_predictors <- predictors[i:0]
                                
                                
                                tmp <- kfoldcv_ranger(wide_analysis_nomissing,
                                                      k=folds,
                                                      weighted=TRUE,
                                                      end_window=j,
                                                      variables=tmp_predictors,
                                                      num.trees=tree_param)
                                
                                
                                cvAUC_combo_5 <-bind_rows(cvAUC_combo_5,
                                                              tmp$cvAUC %>%
                                                                      mutate(added=predictors[i],
                                                                             type=measure,
                                                                             isotype=iso
                                                                      )
                                )
                                
}}}}


write_rds(cvAUC_combo_5,
          paste0(path,"cvAUC_importance/5_cvAUC_combo.rds")
          )



##### Section 6: Estimate cvAUC across different assay types [without CTXB] ---------

#Remove all ctxB markers from pred list 3 and rerun
pred_list_6 <- lapply(pred_list_4, 
                      function(x) {
                              x[!str_detect(x,"CtxB")]
                      }
                      )


marker_cvAUC_list_noctb <- list()

for (i in 1:length(pred_list_6)){
        for(j in times){
                
                cat(names(pred_list_6)[i]," window: ",j,"\n")
                
                marker_cvAUC_list_noctb[[names(pred_list_6)[i]]][[as.character(j)]] <-
                        kfoldcv_ranger(wide_analysis_nomissing,
                                       k=folds,
                                       weighted=TRUE,
                                       end_window=j,
                                       variables=pred_list_6[[i]],
                                       num.trees=tree_param)
                
        }
}

write_rds(marker_cvAUC_list_noctb,
          paste0(path,"cvAUC_importance/6_cvAUC_list_noctb.rds")
)

##### Section 7: Estimate cvAUC across combinations in all isotypes [without CTXB]---------

cvAUC_combo_noctxb_7 <- data.frame()
for(measure in c("NetMFI","RAU")){
        for(iso in c("IgG","IgM","IgA")){
                ###order of importance
                predictors <-paste(measure,iso, c(
                                                  "OgawaOSPBSA",
                                                  "InabaOSPBSA",
                                                  "TcpA","Sialidase","VCC"),
                                   sep="_"
                )
                for (i in 1:length(predictors)){
                        for(j in times){
                                
                                cat(predictors[i]," window: ",j,"\n")
                                
                                
                                tmp_predictors <- predictors[i:0]
                                
                                
                                tmp <- kfoldcv_ranger(wide_analysis_nomissing,
                                                      k=folds,
                                                      weighted=TRUE,
                                                      end_window=j,
                                                      variables=tmp_predictors,
                                                      num.trees=tree_param)
                                
                                
                                cvAUC_combo_noctxb_7 <-bind_rows(cvAUC_combo_noctxb_7,
                                                          tmp$cvAUC %>%
                                                                  mutate(added=predictors[i],
                                                                         type=measure,
                                                                         isotype=iso
                                                                  )
                                )
                                
                        }}}}


write_rds(cvAUC_combo_noctxb_7,
          paste0(path,"cvAUC_importance/7_cvAUC_combo_noctxb.rds")
)



