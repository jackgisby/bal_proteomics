#'Function which performs a standard linear model per protein
#' @export
#' @import dplyr stringr tidyr purrr

olink_lm <- function(df,                        
                       variable,                  
                       outcome="NPX",              
                       covariates = NULL,         
                       return.covariates=F,
                       return.models=FALSE,
                       verbose=T,
                       reorder=NULL
) {  
    
    if(missing(df) | missing(variable)){
        stop('The df and variable arguments need to be specified.')
    }
    
    withCallingHandlers({  # muffles particular errors
        
        #Allow for :/* notation in covariates
        variable <- gsub("\\*",":",variable)
        if(!is.null(covariates)) covariates <- gsub("\\*",":",covariates)
        
        add.main.effects <- NULL
        if(any(grepl(":",covariates))){
            tmp <- unlist(strsplit(covariates,":"))
            add.main.effects <- c(add.main.effects,setdiff(tmp,covariates))
            covariates <- union(covariates,add.main.effects)
        }
        if(any(grepl(":",variable))){
            tmp <- unlist(strsplit(variable,":"))
            add.main.effects <- c(add.main.effects,setdiff(tmp,variable))
            variable <- union(variable,unlist(strsplit(variable,":")))
            variable <- variable[!grepl(":",variable)]
        }
        
        #If variable is in both variable and covariate, keep it in variable or will get removed from final table
        covariates <- setdiff(covariates,variable)
        add.main.effects <- setdiff(add.main.effects, variable)
        
        #Variables to be checked
        variable_testers <- intersect(c(variable,covariates), names(df))
        
        ##Remove rows where variables or covariate is NA (cant include in analysis anyway)
        removed.sampleids <- NULL
        for(i in variable_testers){
            removed.sampleids <- unique(c(removed.sampleids,df$SampleID[is.na(df[[i]])]))
            df <- df[!is.na(df[[i]]),]
        }
        
        #Not testing assays that have all NA:s
        all_nas <- df  %>%
            group_by(GeneID) %>%
            summarise(n = n(), n_na = sum(is.na(!!rlang::ensym(outcome)))) %>%
            ungroup() %>%
            filter(n == n_na) %>%
            pull(GeneID)
        
        
        if(length(all_nas) > 0) {
            
            warning(paste0('The assays ',
                           paste(all_nas, collapse = ', '),
                           ' have only NA:s. They will not be tested.'),
                    call. = F)
            
        }
        
        ##Convert character vars to factor
        converted.vars <- NULL
        num.vars <- NULL
        for(i in variable_testers){
            if (i == variable & !is.null(reorder)) {
                df[[i]] <- ordered(df[[i]], levels=reorder)
                converted.vars <- c(converted.vars,i)
            } else if(is.character(df[[i]])){
                df[[i]] <- factor(df[[i]])
                converted.vars <- c(converted.vars,i)
            } else if(is.numeric(df[[i]])){
                num.vars <- c(num.vars,i)
            }
        }
        
        
        #Not testing assays that have all NA:s in one level
        #Every sample needs to have a unique level of the factor
        
        nas_in_var <- character(0)
        
        if(!is.null(covariates)){
            factors_in_df <- names(df)[sapply(df, is.factor)] 
            single_fixed_effects <- c(variable,
                                      intersect(covariates,
                                                factors_in_df))
        }else{
            single_fixed_effects <- variable
        }
        
        if(!is.null(covariates)){
            formula_string <- paste0(outcome, "~", 
                                     paste(variable,collapse="*"),
                                     "+", 
                                     paste(covariates, sep = '', collapse = '+'))
        }else{
            
            formula_string <- paste0(outcome, "~", paste(variable,collapse="*"))
        }
        
        #Get factors
        fact.vars <- sapply(variable_testers, function(x) is.factor(df[[x]]))
        fact.vars <- names(fact.vars)[fact.vars]
        
        
        #Print verbose message
        if(verbose){
            if(!is.null(add.main.effects) & length(add.main.effects) > 0){
                message("Missing main effects added to the model formula: ",
                        paste(add.main.effects,collapse=", "))
            }
            if(!is.null(removed.sampleids) & length(removed.sampleids) >0){
                message("Samples removed due to missing variable or covariate levels: ",
                        paste(removed.sampleids,collapse=", "))
            }
            if(!is.null(converted.vars)){
                message(paste0("Variables and covariates converted from character to factors: ",
                               paste(converted.vars,collapse = ", ")))
            }
            if(!is.null(num.vars)){
                message(paste0("Variables and covariates treated as numeric: ",
                               paste(num.vars,collapse = ", ")))
            }
            message("Linear mixed effects model fit to each assay: ",formula_string)
        }
        
        if(!is.null(covariates) & any(grepl(":", covariates))){
            covariate_filter_string <- covariates[str_detect(covariates, ':')]
            covariate_filter_string <- sub("(.*)\\:(.*)$", "\\2:\\1", covariate_filter_string)
            covariate_filter_string <- c(covariates, covariate_filter_string)
            
        }else{
            covariate_filter_string <- covariates
        }
        
        if (return.models) {
            #make LMM
            lmer_model<- df %>%
                filter(!(GeneID %in% all_nas)) %>%
                filter(!(GeneID %in% nas_in_var)) %>%
                group_by(Assay, GeneID, Uniprot, Panel) %>%
                group_map(~lm(as.formula(formula_string), data=.x))
            
            mod_names <- df %>%
                filter(!(GeneID %in% all_nas)) %>%
                filter(!(GeneID %in% nas_in_var)) %>%
                group_by(Assay, GeneID, Uniprot, Panel) %>%
                group_data() %>% 
                dplyr::select("GeneID")
            
            names(lmer_model) <- as.vector(mod_names$GeneID)
            
            return(lmer_model)
            
        } else {
            ##make LMM
            lmer_model<-df %>%
                filter(!(GeneID %in% all_nas)) %>%
                filter(!(GeneID %in% nas_in_var)) %>%
                group_by(Assay, GeneID, Uniprot, Panel) %>%
                group_modify(~tidy(anova(lm(as.formula(formula_string), data=.x)))) %>%
                ungroup() %>%
                mutate(covariates = term %in% covariate_filter_string) %>% 
                group_by(covariates) %>% 
                mutate(Adjusted_pval=p.adjust(p.value,method="fdr")) %>%
                mutate(Threshold  = ifelse(Adjusted_pval<0.05,"Significant","Non-significant")) %>%
                mutate(Adjusted_pval = ifelse(covariates,NA,Adjusted_pval),
                       Threshold = ifelse(covariates,NA,Threshold)) %>% 
                ungroup() %>% 
                dplyr::select(-covariates) %>% 
                arrange(p.value)
            
            if(return.covariates){
                return(lmer_model)
            } else{
                return(lmer_model %>% filter(!term%in%covariate_filter_string))
            }
        }
        
    }, warning = function(w) {
        if (grepl(x = w, pattern = glob2rx("*not recognized or transformed: NumDF, DenDF*")) |
            grepl(x = w, pattern = glob2rx("*contains implicit NA, consider using*"))){
            invokeRestart("muffleWarning") 
        }
    })
    
}

#' simple standard differential expression analysis plot
#' @export
#' @import ggplot2 dplyr

simple_lm_de <- function(
    long, 
    case_control_to_remove=NA, 
    variable="case.control", 
    covariates=NULL,
    logp_label=6,
    fc_label=c(-0.5, 1),
    plot_xlab="Log2 Fold Change (P-N)",
    plot_title="plasma case/control",
    return_models=FALSE
) {
    # remove a particular case/control group prior to analysis
    if (!is.na(case_control_to_remove)) {
        long <- dplyr::filter(long, case.control != case_control_to_remove)
    }
    
    # get model outputs in i) dataframe format ii) the models themselves
    cctrl_cov_pvals <- olink_lm(long, variable, covariates=covariates)
    cctrl__cov_models <- olink_lm(long, variable, covariates=covariates, return.models = TRUE)
    
    if (return_models) {
        return(cctrl__cov_models)
    }
    
    cctrl_cov_pvals <- cctrl_cov_pvals[!is.na(cctrl_cov_pvals$p.value),]
    
    # calculate fold change
    fc <- sapply(cctrl__cov_models, function(lm_model) {return(lm_model$coefficients[2])})

    fc_df <- data.frame(GeneID=names(cctrl__cov_models), fc=fc)
    cctrl_cov_pvals <- left_join(cctrl_cov_pvals, fc_df)
    cctrl_cov_pvals$logp <- -log10(cctrl_cov_pvals$Adjusted_pval)
    
    # volcano plot
    print(ggplot(cctrl_cov_pvals, aes(x=fc, y=logp, colour=factor(ifelse(Adjusted_pval < 0.05, "<0.05", "NS"), levels=c("NS", "<0.05")))) +
        geom_point(alpha=0.5) +
        ylab("-log10(Adjusted Pvalue)") +
        xlab(plot_xlab) +
        labs(colour = "Adjusted Pvalue") +
        ggtitle(plot_title) +
        geom_text_repel(data=subset(cctrl_cov_pvals,
                                    logp > logp_label | 
                                        ((fc > fc_label[2] | fc < fc_label[1]) & Adjusted_pval < 0.05)),
                        aes(fc, logp, label = GeneID), size = 3, color="steelblue"))
    
    return(cctrl_cov_pvals)
}
