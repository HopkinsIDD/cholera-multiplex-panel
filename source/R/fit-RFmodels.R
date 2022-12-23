
source("source/R/packages.R")
source("source/R/utils.R")


analysisData <- read_rds("data/raw_data/analysisData.rds") 
wide_analysis <- analysisData$new_wide_data

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

traditional_predictors <- analysisData$new_long_data %>% 
        filter(test_type == "Vibriocidal"|
                       (test_type=="ELISA" &
                                isotype %in% c("IgG","IgA"))
        ) %>%
        distinct(full_test) %>%
        pull(full_test) %>%
        c("age","sex","blood")



var_list <- list(
        "Traditional with nonImmuno"= traditional_predictors,
        "Traditional"= traditional_predictors[!str_detect(traditional_predictors,"age|sex|blood")],
        
        "All NetMFI with nonImmuno" = form_NetMFI,
        "All NetMFI IgG with nonImmuno"= form_NetMFI[str_detect(form_NetMFI,"IgG|age|sex|blood")],
        
        "All NetMFI" = form_NetMFI[!str_detect(form_NetMFI,"age|sex|blood")],
        "All NetMFI IgG"= form_NetMFI[str_detect(form_NetMFI,"IgG") & !str_detect(form_NetMFI,"age|sex|blood")]
)

var_list[["Reduced Panel NetMFI with nonImmuno"]]<- var_list$`All NetMFI IgG (with non-Immuno)`[str_detect(var_list$`All NetMFI IgG (with non-Immuno)`,"CtxB|OSP|age|sex|blood")]
var_list[["Reduced Panel NetMFI"]]<- var_list$`All NetMFI IgG`[str_detect(var_list$`All NetMFI IgG`,"CtxB|OSP")]

var_list[["All RAU with nonImmuno"]] <- form_RAU
var_list[["All RAU IgG with nonImmuno"]] <-  form_RAU[str_detect(form_RAU,"IgG|age|sex|blood")]

var_list[["All RAU"]] <-  form_RAU[!str_detect(form_RAU,"age|sex|blood")]
var_list[["All RAU IgG"]]  <- form_RAU[str_detect(form_RAU,"IgG") & !str_detect(form_RAU,"age|sex|blood")]

var_list[["Reduced Panel RAU with non-Immuno"]]<- var_list$`All RAU IgG (with non-Immuno)`[str_detect(var_list$`All RAU IgG (with non-Immuno)`,"CtxB|OSP|age|sex|blood")]
var_list[["Reduced Panel RAU"]]<- var_list$`All RAU IgG`[str_detect(var_list$`All RAU IgG`,"CtxB|OSP")]


fit_list <- list()

windows <- c(45,120,200,300) 
var_names <- names(var_list)
weighted_options <- c(TRUE,FALSE)

counter <- 0 
n_models <- length(windows) * length(var_names) * length(weighted_options)


for(v in var_names){
        for(w in windows){
                for(weighted in weighted_options){

                #add outcome variable to dataframe
                outcome <- paste0("inf_",w)
                fit_df <-  wide_analysis %>%
                         addInfectionWindow(end_window=w) %>%
                        #limit to individuals with all datapoints
                         select(all_of(c(outcome,
                                         "day","day_actual","status",
                                         var_list[[v]]
                                         ))) %>%
                                        drop_na()
                
                cat(nrow(fit_df),"\n")
                

                #make weights if you want a weighted fit
                inside_w = NULL
                if(weighted){
                        inside_w <- fit_df %>%
                                getWeight(end_window = w)  %>%
                                pull(weight)
                        
                }
                
                #fit rf model
                tmp_fit <- ranger::ranger(formula=make_formula(paste0("inf_",w),
                                                                  var_list[[v]]),
                                             data=fit_df,
                                             num.trees = tree_param,
                                             case.weights = inside_w,
                                             replace=TRUE
                )
                
                
                # list labels
                window_label <- paste0(w,"-days")
                predictor_label <-v
                if(!weighted) weight_label <- "Unweighted"
                if(weighted) weight_label <-  "Weighted"
                
                
                # Create a descriptive list object
                tmp_list <- list(
                        `Infection Window` = window_label,
                        `Predictor Group` = predictor_label,
                        `Predictors` = var_list[[v]],
                        `Weighted` = weight_label,
                        `Fit Data` =  fit_df, 
                        `Model Fit` = tmp_fit
                        
                )
                
                #add fit to the list
                counter <- counter+1
                
                
                #print combination for tracking
                cat(paste(counter,"of",n_models,":",
                          window_label,
                          predictor_label,
                          weight_label),"\n")
                
                
                #add files to the list
                fit_list[[window_label]][[weight_label]] <- tmp_list
                
                
                }}
        
        #define file name
        file_name <- paste0(str_replace_all(predictor_label," ","_"),".rds")
        #save list of fits with file name
        write_rds(fit_list,paste0("data/generated_data/rf_fits/",file_name))
        
}



