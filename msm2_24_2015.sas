
/*   
     MSM (24 Feb 2015)
     Miguel Hernan, Roger Logan, Eric Tchetgen, James Robins
     Copyright, 2005, President and Fellows of Harvard College
     
     MSM is a set of SAS macros to fit a marginal structural Cox model.
     The parameters of the structural Cox model are estimated by fitting 
     weighted Cox mode using inverse probability weights.
     The weighted Cox model is approximated by a weighted pooled logistic 
     regression model.

     Reference (introductory level):
     Hernán MA, Brumback B, Robins JM. Marginal structural models to estimate 
     the causal effect of zidovudine on the survival of HIV-positive men. 
     Epidemiology 2000;11:561-570.


    Note: when including censoring the following assumptions are made

    1) when event = 1 or 0 then all censoring variables are defined to be 0.
    2) when event = . then at least one of the censoring variables should be 1, others can be 0 or missing 
       depending of ordering of censoring variables.
    3  The only possible exception to (2) is at the end of followup where we do not need to model the censoring variable


*/




%macro msm(

     /* INPUT DATA */
    data=,              /* data set name */
    id=,                /* subject unique identifier (variable name) */
    time=,              /* time of follow-up (variable name)
                           first observation for each subject must be time=0 */

     /* structural model */
    outcMSM=,           /* time-varying outcome variable 1: event, 0: no event 
                                 missing if cens=1 or cens2=1 or cens3=1 */
    AMSM=,              /* time-varying exposure variable(s). This has to be a dichotomas variable  (same as A)
                           for current treatment if using doubly robust model */
    AMSM_type  = 0 ,  /* default type of function of treatment will be AMSM = A 
                         1 = A_m and sum(A_j , j = 1 to m-1)/(m-1)
                         2 = categorical variable, which is a function of time treated. For curves submacro 
                             user must list the times where the function changes. 
                         3 = user defined function of treatment. For curves macro this
                             will need to be coded by user */ 
    AMSM_knots = ,      /* when amsm_type = 2, this lists the times where amsm changes levels */
    AMSM_class = ,      /* categorical variables in AMSM (other than two levels) reference level will be first level */
    covMSMbv = ,        /* baseline variables */
    covMSMbh= ,         /* baseline hazard, functions of time (e.g., splines for time) */  
    user_function_of_time = 0 , /* this allows the user to specify their own function of time for variables listed in covMSMbh 
                                   for the outcome model and used in the calculation of the survival probabilities,   Variables 
                                   should be listed in covMSMbh. When user_function_of_time = 1, the user will need to supply a macro 
                                   called %user_function_of_time that creates all variables listed in covMSMbh which can be functions 
                                   of the time variable and any baseline covariate listed in covMSMbv.  */       
    contMSMbh = 1,      /* all variables listed in baseline hazaard variables are continuous */
  
    time_knots = ,      /* if non-missing, these are used for creating covMSMbh time categorical/ spline variables for 
                           final model and curves submacro */ 
    classMSMbv = ,      /* categorical variables in covMSMbv */
    classMSMbh = ,      /* categorical variables in covMSMbh */
    inter_MSM= ,        /* interactions with variables listed in AMSM */
    inter_time= ,       /* interactions with variables listed in covMSMbh */

     /* models for treatment weights */


    A=,                 /* time-varying outcome variable 1: treated, 0: untreated
                                 if missing, weights set to 1  */
    covAd=,             /* baseline and time-varying variables in treatment model 
                                 for denominator of weights */
    classAd=,           /* list of categorical variables in covAd */
    covAn=,             /* baseline variables in treatment model for numerator of 
                                 weights */
    classAn=,           /* list of categorical variables in covAn */
    eligible=,           /* eligibility variable for treatment 1: yes, 0: no If 0 then pA_d and pA_n will
                           be set to 1*/
    A_censoring = ,      /* functions of treatment that should be included in censoring models. For the non-natural course models
                            either A (when A_censoring is missing) or A_censoring will be included in both the numerator and denominator
                            models. When using the natural course only A_censoring will be included in the denominator models. No treatment 
                            variables should be included in the numerator models. */
     /* models for censoring weights */
    Cens=,              /* time-varying outcome variable 1: yes, 0: no */
    covCd=,             /* baseline and time-varying variables in model for denominator
                                 of weights, should not include any treatment variables ( or functions of A ) */
    classCd=,           /* categorical variables in covCd */
    covCn=,             /* baseline variables in model for denominator of weights , 
                           should not include any treatment variables (or functions of A ) */
    classCn=,           /* categorical variables in covCn */
    eligible_cens = ,   /* eligibility variable for Cens 1: yes, 0: no If 0 then pC_d and pC_n will
                           be set to 1*/
     /* models for censoring2 weights */
    Cens2=,              /* time-varying outcome variable 1: yes, 0: no if missing then P(Cens2=0|Cov2) = 1 */
    covC2d=,             /* baseline and time-varying variables in model for denominator
                                   of weights , should not include any treatment variables (or functions of A ) */
    classC2d=,           /* categorical variables in covC2d */
    covC2n=,             /* baseline variables in model for denominator of weights , 
                            should not include any treatment variables (or functions of A )*/
    classC2n=,           /* categorical variables in covC2n */
    eligible_cens2 =  , /* eligibility variable for Cens2 1: yes, 0: no If 0 then pC2_d and pC2_n will
                           be set to 1*/
     /* models for censoring3 weights */
    Cens3=,              /* time-varying outcome variable 1: yes, 0: no */
    covC3d=,             /* baseline and time-varying variables in model for denominator 
                                   of weights , should not include any treatment variables (or functions of A ) */
    classC3d=,           /* categorical variables in covC3d */
    covC3n=,             /* baseline variables in model for denominator of weights ,
                             , should not include any treatment variables (or functions of A )*/
    classC3n=,           /* categorical variables in covC3n */
    eligible_cens3= , /* eligibility variable for Cens3 1: yes, 0: no If 0 then pC3_d and pC3_n will
                           be set to 1*/
    use_stabilized_wts = 1 , /* include the numerator model in the calculation of the weights */
    user_defined_weights = 0, /* use user defined weights in analysis, skips calculate_weights macro. For the outcome model, the weights should
                                 be contained the the variable "sw" and if the natural course is being used the corresponding weights should be 
                                 be contained in the variable "natural_wt". The natural course weight should be for censoring only. Any truncation 
                                 of the weights should be done prior to calling the macro.*/

    /* interactions in hazard model */ 
   
     /* DATA ANALYSIS SPECIFICATIONS */
    class_ref = first ,    /* method for coding levels for binary dummy variables for class variables.
                           0 default coding of {0 1}
                           1 is similar to param = ref option in class statement { 1, -1} coding
                          */
  
    use_genmod = 0 ,    /* use proc genmod for final model */
   
    msm_var = 1,        /* calculate analytic estimate of variance  */
    truncate_weights = 0,  /*  0 default: use all weights 
                           1 truncates weights below and above 1 and 99 percentile
                           2 truncates  weights below and above 2 and 98 percentile, etc*/ 
    user_defined_limits = , /* user defined lower and upper limits for the weights. This must be a list with two numbers
                               separated with a space: user_defined_limits = lower_limit upper_limit. this will only be used 
                               when truncate_weights = 1 and user_defined_limits is not missing. If only one entry is listed or the 
                               list is impty then the method will  used the percentile values given in the truncate_weights option.
                             */
    bootstrap= 0,       /* bootstrap estimate of variance */
    
    nboot= 200,         /* number of bootstrap samples */
    bootstart = 0 ,
    bootend = ,
    bootlib =   work , /* libname for saving  bootstrap results */
    bootname = bootstrap , /* name of dataset for holding bootstrap results */
    bseed= 12345,       /* random number seed for boostrap */


    /* COMPUTATIONAL SPECIFICATIONS */
    just_models = 0,    /* estimate weights only */
    override_warnings= 0 , /* override the warnings from the weight models */
    logistic_options = , /* options for proc logistic, in model statement */
    inv_type = inv , /* When msm_var = 1, inv_type is the proc iml function that is used in calculating inverse of required matrices, possible choices are 
                        inv = usual inverse or ginv = generalized inverse. Due to potential ill conditioned matrices one 
                        function may work when the other does not. Information concerning the determinant and 
                        estimated condtion number are supplied in the log file. */
  
  
    /* OUTPUT OPTIONS */
    no_listing = 0 ,    /* suppress all listing output , useful when calling macro from a loop   */ 
    print_opt= 1,       /* print out results for the optimization method when est_method = 1 */
    optn_level= 1,      /* output level for Nelder-Mead Simplex method (est_method = 1)
                           0 = none 5 = all */
    print_boot= 0,      /* print out the indivudal model results for each bootstrap sample */
 
    save_results= 1,    /* save bootstrap results into sas data sets */
    save_weights= 0,    /* save the weights */ 


    libsurv= work   ,      /* User defined libname for saving results, work directory will be default */
    results_name= _results,  /* name of file for saving results */
    weights_name = _weights,
    
    rnseed= ,
    debug= 0 ,           /* 1 = keep intermediary data sets, 0 = delete intermediary data sets */

            /* new variables for incorporating Sarah's/Jessica's msm code 
               TEMPORARY HOLDING SPACE
               */
   
    use_natural_course = 0, /* for natural course analysis that remves the treatment weights and AMSM from
                              weighted models */
    survival_curves = 0 ,  /* create survival curves for always treated versus never treated */
     
    competing_risk = ,     /* additional outcome simulated in survival curves (0-1 outcome variable) will use
                               same model statement as in the outcome model */
    riskdiff1 = ci_treat, 
    riskdiff0 = ci_nontreat,  /* select from ci_nontreat ci_nc for calcluating ci_treated - &riskdiff0 */
    print_ci = 0 ,
    time_start = 0 ,
    time_end = 60,
    time_step = 1 
    
     );

     /*
     The program is structured into 9 main modules:

     %MSM, the main driver, inputs data, analysis and computing 
     specifications

   
     %SETDATA prepares the data to be analyzed by the macro %MSMCOX

     %CALCULATE_WEIGHTS, calculate models for treatment and censoring that are
     used in the calculation of the stabalized weights used the the structural model

     %MSMCOX estimates the parameters of a marginal structural Cox model

     %INITIALIZE_MSM_VAR, initializes macro variables that will be used in the calculation of the analytic
     variance under the non-doubly robust marginal structural Cox model
  
     %CALCULATE_MSM_VAR, calculates the analytic variance of the parameters of  a non-doubly robust
     marginal structural Cox model         

     Utilitiy macros:

     %NUMARGS, auxillary macro for finding number of entries in a macro list

     %CREATE_CLASS_LIST, create a list of class variables which appear in a list of varialbes 

     %REMOVE_CLASS,  remove entries in one list  from a second list

     %CONVERT_CLASS, converts class variables to numeric variables so they can be used
     in the models for E[delta * f(H)| covDelta] in %SNMAFT

     %REMOVE_EXTRANEOUS, removes entries in a class list if they are not contained in
     the corresponding covariate list, e.g. classDelta and covDelta

     

     */

      

     %local bsample bootstrap nboot ;
     %local  individuals chilevel cov_betas n_betas num_events  events_used  ;
     %local time_var comptime_var scores modelAn   modelCd modelCn modelCd2 modelCn2 modelCd3 modelCn3;
     %local class_inter_MSM class_inter_time class_hazard for_emppest  ;
     %local interAd interAn interCd interCn interC2d interC2n interC3d interC3n ;
     %local treatment_weights class_upper class_lower beta_bsample0 just_models_holder varnametmp ;
     %local cov_betas_nc cov_betas_cr ;
     %local inter_def ndef_time pred inter_def_amsm ndef_amsm _interactions _interactions_bh covbh_inter_def amsm_inter_def ;
     %local idname lower_limit upper_limit ;
   
    
    
     %let idname = &id ;      
    
     /* consistency checks and input error management (under development) */
          
     %if &user_defined_weights = 1 %then %let msm_var = 0 ;
 
     %let treatment_weights = 1; /* was used to allow for only censoring weights , now always create treatment weights */
     %let just_models_holder = &just_models ;

     %if &msm_var = 1 %then %let use_genmod = 0 ;

     %if &survival_curves = 1 %then %do;
         %let use_genmod = 0 ;
         %let msm_var = 0 ;
     %end;
     %else %do;
         %let competing_risk = ;
         %let use_natural_course = 0 ;
         %let simulate_censoring = 0;
         %let user_function_of_time = 0;
     %end;

     %if %bquote(&A) =  %then %do;
          %let covAd = ;
          %let classAd = ;
          %let covAn = ;
          %let classAn = ; 
          %let treatment_weights = 0;
                      
     %end;

     %if %bquote(&Cens)= %then %do;
          %let covCd = ;
          %let classCd = ;
          %let covCn = ;
          %let classCn = ;          
     %end;

     %if %bquote(&cens2) = %then %do;
          %let covC2d = ;
          %let classC2d = ;
          %let covC2n = ;
          %let classC2n = ;
     %end;    

     %if %bquote(&cens3) = %then %do;
          %let covC3d = ;
          %let classC3d = ;
          %let covC3n = ;
          %let classC3n = ;
     %end;    

     %local use_spines knots nknots ntime ;
     %let use_splines = 0 ;
     %if %bquote(&time_knots)^= and &contMSMbh = 1 %then %do;

        %let ntime = %numargs(&covMSMbh) ;
        %let use_splines = 1 ;
        %let nknots = %numargs(&time_knots);
        %let knots = %scan(&time_knots,1,%str( ));
        %do i = 2 %to &nknots ;
            %let knots = &knots , %scan(&time_knots,&i,%str( ));
        %end;
     %end; 


     %let get_cov = 0 ;
     %if %eval(&bootstrap) = 0 %then %do ;
          %let nboot = 0;
          %let bootstart = 0 ;
          %let bootend = 0 ; 
          %let k_step_update = 0 ; /* do not use any update to beta */
     %end;
     %else %if %eval(&bootstrap) = 1 %then %do ;       
          %let MSM_var = 0;   /* non double robust analytic estimate of variance */          
          %if %bquote(&bootend)= %then %let bootend=&nboot ;
          %if &bootstart = 0 AND &bootend > 1 %then %let get_cov = 1;                        
     %end;

     %if %eval(&just_models) = 1 %then %do;       
          %let msm_var = 0;
          %let get_cov = 0;
     %end;

     %let var_calc = 0;
     %if %eval(&get_cov) = 1 | %eval(&msm_var) = 1  %then %let var_calc = 1;

      %if %eval(&truncate_weights)^= 0 %then %do;
          %local check lower_limit upper_limit ;
          %let check = %numargs(&user_defined_limits) ;
          %if &check = 2 %then %do;
               %let lower_limit = %qscan(&user_defined_limits,1,%str( )) ;
               %let upper_limit = %qscan(&user_defined_limits,2,%str( ));
               %let user_defined_limits = 1 ;
           %end;
           %else %do ;
               %let user_defined_limits = 0 ;
               %let lower_limit = ;
               %let upper_limit = ;
           %end;
       %end;                       
                
     %setdata ;         
     %msmcox; 
     



     /* output results */
     %if %eval(&just_models) = 0 %then %do;       
 
             %if &bootstrap = 0 %then %do ;

                proc print data = _msmsw noobs ;
                %if &msm_var = 1 %then title "Results for marginal structural model with stabilized weights for &outcMSM using analytic variance" ;
                %if &msm_var = 0  %then %do;
                     %if &use_genmod=0 %then title "Results for marginal structural model with stabilized weights for &outcMSM using proc logistic " ;;
                     %if &use_genmod = 1 %then title "Results for marginal structural model with stabilized weights for &outcMSM using robust variance" ;;
                %end;
                %if &truncate_weights = 1 %then title2 "Truncating the weights so that &lower_limit <= sw <= &upper_limit ";; 
                var name estimate std lb ub z p_value ;
                run;

                
             %end;
             
 
              %if &get_cov = 1 & &bootstrap = 1 %then %do;

                
                   %bootstrap_results(
                       bootlib = &bootlib ,
                       bootname= &bootname,
                       modeltype = sw ,
                       numparts = 1,
                       samplestart = 0,
                       numboot = &bootend );

                
                    %if &use_natural_course=1 %then %do;    

                         %bootstrap_results(
                                  bootlib = &bootlib ,
                                  bootname= &bootname,
                                  modeltype = nc ,
                                  numparts = 1,
                                  samplestart = 0,
                                  numboot = &bootend,
                                  time = &time,
                                  outcMSM = &outcMSM ,
                                  use_natural_course = 1 );

 
                    %end;
                    %if %quote(&competing_risk)^= %then %do;
                          %bootstrap_results(
                                 bootlib = &bootlib ,
                                 bootname= &bootname,
                                 modeltype = cr ,
                                 numparts = 1,
                                 samplestart = 0,
                                 numboot = &bootend,
                                 time = &time ,
                                 outcMSM = &outcMSM,
                                 competing_risk = &competing_risk 
                                 );                         
                    %end;
                    %if &survival_curves = 1 %then %do;
                           %bootstrap_results(
                                               bootlib=&bootlib,
                                               bootname=&bootname,
                                               modeltype =surv ,
                                               numparts=1,
                                               samplestart = 0 ,
                                               numboot= &bootend ,
                                               time = &time,
                                               outcMSM = &outcMSM,
                                               use_natural_course = &use_natural_course,
                                               riskdiff1 = &riskdiff1 ,
                                               riskdiff0 = &riskdiff0  
                                                  );

                        
          
                    %end;
               %end; /* get_cov = 1 */

               %if &survival_curves = 1 %then %do;
			   
			            %let _resdata = _onegraph_;
					    %if &boostrap = 1 %then %let _resdata = _bgraphsr ;

          			    %plot_graphs(datain =  &_resdata ,
                              time = &time ,    
                              outcMSM= &outcMSM ,
                              natural_course = &use_natural_course  , 
                              sixgraphs=0 ,
                              gfilename=fourgraphs.pdf ,
                              title1=,title2=,title3=,titledata=,tsize=0.75 ,
                              frombootstrap = &bootstrap ) ;

               %end;
 

               %if %eval(&save_results) = 1 %then %do;
    
                    data &libsurv..&results_name;                     
                    set _msmsw  ;                    
                    run;
                 
                    
                    data &libsurv..betas_covariance  ;
                    %if &bootstrap = 0 %then %do;
                           set varbetas_sw ;
                    %end;
                    %else %if &bootstrap = 1 AND &survival_curves = 0 %then %do;
                        set _cov_sw (where = (_TYPE_="COV") rename = (_NAME_ = Variable));
                        drop _TYPE_ ;
                    %end ;
                    run;
                  
               %end; 
               

                %if &save_weights = 1 %then %do;
                     data &libsurv..&weights_name;
                     set _weights;
                     run;
                %end;
                %if %eval(&debug) = 0 %then %do;

                     proc datasets library = work nolist;
                     delete _cero0_ _msmsw_all  _weights  
                     _idholders_ _idsamples %if &use_genmod =0 %then  _varnames_  _for_names ;                      
                      ;
                    quit;
               %end;
          %end;
          title;
     
%mend msm ;




%macro setdata;

  
   
     proc sort data= &data out=_cero_ sortsize=256M ; 
     by &id &time; 
     run;
     
    

    
     %let class_inter_MSM = ;
     %let class_inter_time = ;
     %if %bquote(&inter_MSM)^= 
                  %then %create_class_list(cov_list = &inter_MSM,
                                        parent_list = &classMSMbh &classMSMbv,
                                        class_list = class_inter_MSM );
     %if %bquote(&inter_time) ^= 
                 %then %create_class_list(cov_list  = &inter_time,
                                        parent_list = &classMSMbv ,
                                        class_list  = class_inter_time );
          
     %let nclass = 13 ;
     %let class_var = %nrstr(  &classAd &classAn &classCd &classCn &classC2d &classC2n &classC3d &classC3n 
                                    &classMSMbh &classMSMbv &class_inter_MSM &class_inter_time &amsm_class);
     %let cov_var   = %nrstr(  &covAd &covAn &covCd &covCn &covC2d &covC2n &covC3d &covC3n 
                                    &covMSMbh &covMSMbv &inter_MSM &inter_time &amsm );
     %let aaaa = %unquote(&class_var) ;
     %if %bquote(&aaaa) ^=  %then %do;
                 %remove_extraneous;
                 %convert_class;
     %end;
              
    
/********** create macro variables needed for any interactions included in weighted models ***/

     %let _interactions = ;
     %let for_emppest = ;
     %let amsm_inter_def = ;
         
               
   
     %local ntmp vars vartmp num_treat _nvar ntreat treat_tmp inter_tmp ;               
               
     %let num_treat = %numargs(&AMSM) ;
     %let _nvar = %numargs(&inter_MSM) ;

     %let ndef_amsm = %eval(&num_treat * &_nvar) ;
     %let inter_count = 0 ;

     %do ntreat = 1 %to %eval(&num_treat);
         %let treat_tmp = %scan(&AMSM,%eval(&ntreat),%str(' '));               
         %let for_emppest = &for_emppest &treat_tmp ;
             %do inter = 1 %to %eval(&_nvar);
                   %let prod_tmp = ;
                   %let def_tmp = ;
                   %let varname_tmp = ;
                   %let inter_count = %eval(&inter_count + 1);
                   %let inter_tmp = %scan(&inter_MSM,%eval(&inter),%str(' '));
                   %let _interactions = &_interactions amsm_inter&inter_count ;
                   %let for_emppest = &for_emppest &treat_tmp.&inter_tmp;
                   %let prod_tmp = &treat_tmp * &inter_tmp ;                      
                   %let def_tmp = amsm_inter&inter_count = &prod_tmp ;
                   %if &inter_count = 1 %then %let amsm_inter_def = &def_tmp  ;
                   %else %let amsm_inter_def = &amsm_inter_def : &def_tmp ; 
              %end;
      %end;
                          

     %if &use_splines = 1 %then %do;
        %let covMSMbh_new = ;
        %let  covMSMbh_orig = &covMSMbh ;
        %let ncovMSMbh_orig = %numargs(&covMSMbh);
        %do itime = 1 %to &ncovMSMbh_orig ;
             %let word = %scan(&covMSMbh,%eval(&itime),%str( ));           
             %let covMSMbh_new = &covMSMbh_new &word ;
             %do iknot = 1 %to %eval(&nknots-2) ;
                 %let covMSMbh_new = &covMSMbh_new &word.&iknot ;
             %end;                                        
         %end;
         %let covMSMbh = &covMSMbh_new ;
     %end;

  
              
     %let num_bh = %numargs(&covMSMbh) ;
     %let _nvar2 = %numargs(&inter_time);
     %let _interactions_bh = ;
     %let for_emppest_bh= ;
     %let covbh_inter_def = ;
     %let ndef_time = %eval(&num_bh * &_nvar2) ;
     %let inter_count = 0 ;
     %do ntime = 1 %to %eval(&num_bh);
          %let time_tmp = %scan(&covMSMbh,%eval(&ntime),%str(' '));                 
          %let for_emppest_bh = &for_emppest_bh &time_tmp ;
          %do inter = 1 %to %eval(&_nvar2);
                %let def_tmp = ;
                %let prod_tmp = ;
                %let varname_tmp = ;
                %let inter_count = %eval(&inter_count + 1);
                %let inter_tmp = %scan(&inter_time,%eval(&inter),%str(' '));
                %let prod_tmp = &time_tmp * &inter_tmp ;
                %let _interactions_bh = &_interactions_bh covbh_inter&inter_count ;
                %let for_emppest_bh = &for_emppest_bh &time_tmp.&inter_tmp;
                %let def_tmp = covbh_inter&inter_count = &prod_tmp ;
                %if &inter_count = 1 %then %let covbh_inter_def = &def_tmp ;
                %else %let covbh_inter_def = &covbh_inter_def : &def_tmp ;
           %end;
     %end;
           

     %let cov_betas = Intercept &amsm &_interactions &covMSMbh &_interactions_bh  &covMSMbv  ;               
     %let n_betas = %numargs(&cov_betas) ;
     


 
     %let vars = &covMSMbh &covMSMbv ;
     %let ntmp = %numargs(&vars );
      
  
    
     data _cero0_ (index = (bsample));
     set  _cero_  end = _end_;
     by &id ;
     retain _newid_ ;
     if _n_ = 1 then _newid_ = 0 ;
     if first.&id then do;
          _newid_ = _newid_ + 1;                   
     end;       
     bsample = 0;     
     if _end_ then  do; 
          call symput('individuals',_newid_);           
     end;
     
     %if &use_splines = 1 %then %do;        
        %do itime = 1 %to &ncovMSMbh_orig ;
           %let word = %scan(&covMSMbh_orig,%eval(&itime),%str( ));
           %rcspline(&word,&knots);                                               
        %end;
     %end;

    


     %if &ndef_amsm > 0 %then %do;

         %do ninter = 1 %to &ndef_amsm ;
            %let word =  %scan(&amsm_inter_def,%eval(&ninter),%str(:)) ;
            &word ;
            label amsm_inter&ninter="%left(%qscan(&word,2,%str(=)))";

         %end;
     %end;
     %if &ndef_time > 0 %then %do ;
         %do ninter = 1 %to &ndef_time ;
           %let word =  %scan(&covbh_inter_def,%eval(&ninter),%str(:)) ;
           &word ;
           label covbh_inter&ninter="%left(%qscan(&word,2,%str(=)))";
         %end;
     %end;

     %if &user_function_of_time = 1 %then %user_function_of_time ;
        
     run;       

     proc datasets library = work nolist;
     modify _cero0_;
         label   
            %do ntreat = 1 %to %eval(&num_treat);
                %let vartmp = %scan(&AMSM,%eval(&ntreat),%str(' '));
                &vartmp = "&vartmp"
            %end;                         
            %do _i = 1 %to %eval(&ntmp);
                   %let vartmp = %scan(&vars,%eval(&_i),%str(' '));
                       &vartmp = "&vartmp"
                   %end;
                       ;
                  
     quit;          

     %let numids = &individuals ;
     
     %let id = _newid_; /* sequential id starting a 1 */

     


  
     data _idholders_ (index = (bsample));
     do bsample = 0 to &nboot ;
         do _newid_ = 1 to &numids;
           output ;
         end;
     end;
     run;

     proc surveyselect data= _idholders_ 
         method = urs
         n= &individuals
         seed = 1232  
         out = _idsamples (keep = bsample _newid_  numberhits  ) 
         outall    
         %if &bootstrap = 0 %then noprint               ;
           ;
      strata bsample ;
      run;

      data _idsamples ; 
      set _idsamples ;
      if bsample = 0 then numberhits = 1 ;
      run;
                 
/*  results data sets */
     data _msmsw_all; 
     run;

     %if &bootstrap = 1 %then %do;
    
         * for bootstraps save individual sample results in possible permanent data sets ;

         data &bootlib..&bootname._sw_&bootstart._&bootend ;
         run;

         %if &survival_curves = 1 %then %do;                          

              data &bootlib..&bootname._surv_&bootstart._&bootend ;
              run;


              
             %if %bquote(&competing_risk) ^= %then %do;
                 data &bootlib..&bootname._cr_&bootstart._&bootend ;
                 run;
             %end;
             %if &use_natural_course = 1 %then %do;
                 data &bootlib..&bootname._nc_&bootstart._&bootend ;
                 run;
             %end;
         %end;    
     %end;

     
%mend;




%macro msmcox; 

     title 'Marginal Structural Cox Model, Stabilized Weights';
    
   
     %do bsample = %eval(&bootstart) %to %eval(&bootend);
           %if &print_boot = 0 AND &bsample > 0 %then %do;
                     ods listing close ;
                     ods results off ;
           %end;
           %put RUNNING BOOTSTRAP SAMPLE &bsample ; 
          data _cero_ ;
          merge _cero0_ _idsamples (keep = _newid_ numberhits bsample 
                                   where = (bsample = %eval(&bsample))
                                    );
          by _newid_  ;
          drop bsample ;
          run;
          

          %if &user_defined_weights = 0 %then %calculate_weights ;
          %else %do;

                data _cero_ ;
                set  _cero_  end = _end_ ;
                by &id &time ;
                if first.&id then _occasion_ = -1;
                _occasion_ + 1 ;
                %if (&msm_var = 0 ) %then  if &outcMSM ne .;;  
                if &outcMSM = 1 then _events_used+1 ;
                if _end_ then call symput('events_used',trim(left(_events_used)));         
                run;
          %end;       


          %if %eval(&just_models) = 0 %then %do;
               %if %eval(&msm_var) = 1 %then %initialize_msm_var ;             
                            
               %if &use_genmod = 0 %then %do;
                                                     
                    proc logistic data=  _cero_ descending outest= emppest_sw (drop = _LINK_   _NAME_ _LNLIKE_   ) covout ;
                    ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
                    ods output ParameterEstimates = pe_sw;

                    model &outcMSM = &amsm  &_interactions &covMSMbh  &_interactions_bh  &covMSMbv 
                               %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
                    freq numberhits ;
                    weight sw;
                     
                    title "Weighted model for outcome &outcMSM";
                    %if &msm_var = 1 %then output out = _formsmvar_ (keep = &id &time hazard) p = hazard ;;
                    run;

                    %checksw(datain = emppest_sw) ;

                    %if &use_natural_course = 1 %then %do;
                          proc logistic data=  _cero_ descending outest= emppest_nc (drop = _LINK_   _NAME_ _LNLIKE_  ) covout ;
                          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
                          ods output ParameterEstimates = pe_nc ;

                          model &outcMSM =&covMSMbh  &_interactions_bh  &covMSMbv  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
                          weight natural_wt;
                          freq numberhits ;
                         
                          title "Weighted model for &outcMSM under the natural course.";                           
                          run;

                           %checksw(datain = emppest_nc) ;

                           %if %bquote(&competing_risk)^= %then %do;

                                proc logistic data=  _cero_ descending outest= emppest_nc_cr (drop = _LINK_   _NAME_ _LNLIKE_  ) covout ;
                                ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
                                ods output ParameterEstimates = pe_nc_cr ;

                                model &competing_risk =   &covMSMbh &_interactions_bh  &covMSMbv %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
                                weight natural_wt;
                                freq numberhits ;
                                
                                title "Weighted model for competing risk &competing_risk under the natural course";                           
                                run;

                                 %checksw(datain = emppest_nc_cr) ;
                          %end;
                    %end;


                    %if %bquote(&competing_risk)^= %then %do;

                            proc logistic data=  _cero_ descending outest= emppest_cr (drop = _LINK_   _NAME_ _LNLIKE_  ) covout  ;
                            ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
                            ods output ParameterEstimates = pe_cr ;

                            model &competing_risk = &amsm  &_interactions &covMSMbh &_interactions_bh  &covMSMbv %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
                            weight sw;
                            freq numberhits ;
                           
                            title "Weighted model for competing risk &competing_risk ";                           
                            run;

                             %checksw(datain = emppest_cr) ;
                    %end;

                    %if %eval(&bootstrap) = 0 & %eval(&msm_var) = 0 %then %do;
                           
                           proc contents data = work.emppest_sw (drop = _TYPE_ _STATUS_ ) 
                              out = _varnames_ (keep = varnum name)
                              /* noprint */;
                           run;
 

 
                            proc sql  noprint  ;
                              select  name
                              into :varnametmp separated by ' '
                              from _varnames_
                              order by varnum  ;
                           quit;

                        
                           data varbetas_sw;
                           set emppest_sw (where = (_TYPE_ = 'COV') drop =  _STATUS_  ) ;
                           drop _TYPE_ ;
                           array cov &varnametmp    ;
                           var = cov[_N_] ;
                           std = sqrt(cov[_N_] );
                           _row_ = _N_ ;
                           run;
                    %end;

                    data _NULL_;
                    set emppest_sw (keep = _STATUS_ );
                    good = substr(_STATUS_,1,1);                
                    call symput('msm_good',good);               
                    run;

              
                    %if %eval(&msm_good) ^= 0 %then %do;                       
                         %let msm_var = 0;
                         %put setting dr = 0 and msm_var = 0 due to problems in msm model ;
                    %end;

                   data emppest_sw;
                   set emppest_sw (drop = _STATUS_);
                   run;

                   data _for_names ;
                   set emppest_sw;
                   run;
                %end;
                %else %do;   

                        proc genmod data = _cero_ descending ;
                        ods exclude GEERCov ;
                        ods output GEEEmpPEst=emppest_g ConvergenceStatus = _g_status_ 
                                      GEERCov = varbetas_sw
                         ;
                        class &id ;
                        model  &outcMSM = &amsm &_interactions   &covMSMbh &_interactions_bh  &covMSMbv / dist = bin link = logit  ;
                        weight sw;
                        freq numberhits ;
                         repeated subject = &id / type = ind sorted  ecovb; 
                        run;
                        ods select all ;

                        data _null_;
                        set _g_status_;
                        if _n_ = 1 then call symput('msm_good',status);
                        run;

                        data emppest_sw;
                        set emppest_g (keep = Parm estimate LowerCL UpperCL ProbZ Stderr Z
                           rename = (parm = name LowerCL = lb UpperCL = ub ProbZ = p_value Stderr = std)) ;
                        _order_ = _n_ ;
                        run;

                        proc contents data = _cero_ (keep = &amsm &_interactions &_interactions_bh  &covMSMbh &covMSMbv)
                          out = eee_cont  (keep = name  label)  noprint;
                        run;

                        proc sql noprint ;
                          create table _kk_(drop = _order_ rename = (label=Notes)) as 
                          select * from emppest_sw 
                          full join  eee_cont
                          on emppest_sw.name=eee_cont.name
                          order _order_ 
                                  ;
                          select Notes into :labels separated by " : " from _kk_ ;
                       quit;

                       data _kk_ ;
                       set _kk_ ;
                       _row_ = _N_ ;
                       run;

                       data varbetas_sw ;
                       merge varbetas_sw (rename = (%do _parm = 1 %to &n_betas ;                                
                                                       prm&_parm = %scan(&cov_betas,%eval(&_parm),%str( ))
                                                     %end; ))
                               _kk_ (keep = name Notes _row_ );
                       by _row_ ;
                       drop rowname ;
                       
                       label %do _parm = 2 %to &n_betas ;
                                    %let word1 =   %scan(&cov_betas,%eval(&_parm),%str( )) ;
                                    %let word2 =   %qscan(&labels,%eval(&_parm-1),%str(":")); 
                                    &word1 =  "&word2"
                              %end;
                              ;
                       if index(Notes,"*")> 0 then name = Notes ; 
                       run;
                                                                     

                %end;

                  
                %if &use_genmod = 0   %then %do;
                          data empsw;
                          set emppest_sw ;
                          run;

                          proc contents data = empsw out=eee_cont(keep = varnum name type label where = (type ne 2)) varnum noprint ;
                          run;

                          proc sort data = eee_cont out=eee_cont(keep = NAME LABEL rename = ( LABEL =Notes )) ;
                          by varnum ;
                          run;               
              
                          data eee_cont ;
                          set eee_cont ;
                          _row_ = _N_ ;
                          run;

                         %if &bootstrap = 1 %then %do;
                            data _onesample_sw  ;
                            set emppest_sw ( where = (_TYPE_="PARMS"));
                            bsample = &bsample ;
                            drop _TYPE_ ;
                            run;

                            data &bootlib..&bootname._sw_&bootstart._&bootend ;
                            set &bootlib..&bootname._sw_&bootstart._&bootend  _onesample_sw;
                            run;

                            %if &use_natural_course = 1 %then %do;
                                 data _onesample_nc  ;
                                 set emppest_nc (drop = _STATUS_ where = (_TYPE_="PARMS"));
                                 bsample = &bsample ;
                                 drop _TYPE_ ;
                                 run;

                                 data &bootlib..&bootname._nc_&bootstart._&bootend  ;
                                 set &bootlib..&bootname._nc_&bootstart._&bootend  _onesample_nc;
                                 run;
                             %end;
                             %if %bquote(&competing_risk)^= %then %do;
                                 data _onesample_cr  ;
                                 set emppest_cr (drop = _STATUS_ where = (_TYPE_="PARMS"));
                                 bsample = &bsample ;
                                 drop _TYPE_ ;
                                 run;

                                 data &bootlib..&bootname._cr_&bootstart._&bootend  ;
                                 set &bootlib..&bootname._cr_&bootstart._&bootend  _onesample_cr;
                                 run;
                             %end;

                             proc datasets library = work nolist ;
                             delete _onesample_sw 
                                %if &use_natural_course = 1 %then _onesample_nc ;
                                %if %bquote(&competing_risk)^= %then _onesample_cr ;
                                ;
                              quit;

                        %end;

                        proc transpose data =  emppest_sw out = emppest_sw (drop = _LABEL_ rename = (col1 = estimate _NAME_ = name));                    
                        run;
                  
                        data emppest_sw ;
                        set emppest_sw ;
                        _row_ = _N_ ;
                        run;

                        %if %eval(&bootstrap) = 0 & %eval(&msm_var) = 0  %then %do;
                              data _kk_ ;
                              merge emppest_sw varbetas_sw( keep = var std _row_ ) ;
                              by _row_ ;
                              label p_value = 'Pr > |Z|';
                              format p_value PVALUE8.4 ;
                              if std = 0 | std = . then do;
                                  z = . ;
                                  p_value = . ;
                              end;
                              else do ;
                                   z = estimate / std ;
                                   p_value = 2*(1-probnorm(abs(z)));
                              end;
                              lb = estimate - 1.96 * std ;
                              ub = estimate + 1.96 * std ;  
                              
                              run;

                            data varbetas_sw ;
                            merge _kk_ (keep = name _row_) varbetas_sw(drop = var std  ) ;
                            by _row_ ;
                            drop _row_ ;
                            run;
                        %end;            
                %end; /* genmond = 0 */
                
                %if  %eval(&msm_var)  = 1 %then %calculate_msm_var;
 

                %if &bootstrap = 0 %then %do;
                     


                    data _msmsw ;
                    %if &use_genmod = 0 %then %do;
                        merge _kk_ eee_cont (keep = Notes _row_) ;
                        by _row_ ;
                    %end; 
                    %else set _kk_ ;;
                    if index(Notes,"*")> 0 then NAME = Notes ; 
                    if bsample = . then bsample = &bsample ;
                    drop _row_ ;
                    run;
                %end;

                %if &survival_curves = 1 %then %do;
                    %curves ;
                    %if &bootstrap = 1 %then %do;
                            data  &bootlib..&bootname._surv_&bootstart._&bootend  ;
                            set  &bootlib..&bootname._surv_&bootstart._&bootend  _onegraph_ ;
                            run ;
                    %end;
                %end;                                           
          %end; /* basic_models = 0 */
          
          %if &bootstrap = 1 %then %let just_models = &just_models_holder ;

          %if &debug = 0 %then %do;
               proc datasets library=work nolist;
               delete  _cero_   emppest_sw %if &bootstrap = 0 %then _kk_ ; eee_cont
                    %if %eval(&just_models) = 0 %then  %do;                         
                         %if &use_genmod = 1 %then  emppest_g _g_status_ ;
                         %if &use_genmod = 0 %then %do;
                                pe_sw  empsw  
                                %if &msm_var = 1 %then  nuisance ;                 
                                %if &bootstrap = 1 %then %do;
                                   
                                    %if &survival_curves = 1 %then _onegraph_ ;
                               %end;
                        %end;  
                    %end; 
                 ;
               quit;
          %end;
     ods listing ;
     ods results on ;
     %end; /* end bsample loop */

     title;

%mend msmcox;

%macro calculate_weights;

    %let modelAd = ;    %let modelAn = ;
    %let modelEd = ;    %let modelEn = ;
    %let modelCd = ;    %let modelCn = ;
    %let modelC2d = ;   %let modelC2n = ;
    %let modelC3d = ;   %let modelC3n = ;
    %let modelCdnc= ;   %let modelCnnc = ;
    %let modelC2dnc= ;   %let modelC2nnc = ;
    %let modelC3dnc= ;   %let modelC3nnc = ; 
    
     
    %if %eval(&treatment_weights) = 1 %then %do;
          proc logistic data = _cero_ descending outest = statAd(keep = _STATUS_ )    ; 
          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
          %if %bquote(&eligible)^= %then where &eligible = 1;;
          model &A = &covAd  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
          freq numberhits;
          output out = _modelAd_ (keep=   &id &time pA_d) p = pA_d;
          title "Model for denominator of stabilized weights for &A";
          run;
          %let modelAd = _modelAd_;

          proc logistic data = _cero_  descending outest = statAn(keep = _STATUS_ )     ;
          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
          %if %bquote(&eligible)^= %then where &eligible = 1;;
          model &A = &covAn  %if %bquote(&logistic_options)^= %then / &logistic_options ;;
          freq numberhits ;
          output out = _modelAn_ (keep=  &id &time pA_n) p = pA_n;
          title "Model for numerator of stabilized weights for &A";
          run;
          %let modelAn = _modelAn_;

    %end;

    %if %bquote(&Cens)^= %then %do;

        
          %let _where_ = ;
          %if %bquote(&eligible_cens) ^= %then %let _where_ = &eligible_cens = 1 ;
           
          %if %bquote(&_where_) ^= %then %let _where_ = where ( &cens ne . ) AND ( &_where_)  ;
          %else %let _where_ = where &cens ne . ;
           

          proc logistic data =  _cero_    outest = statcd(keep = _STATUS_ )    ;           
          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;       
          &_where_ ;
          model &Cens =  %if %bquote(&A_censoring)^= %then &A_censoring ; %else &A ; &covCd  %if %bquote(&logistic_options)^= %then / &logistic_options ;  ;
          freq numberhits ;
          output out = _modelCd_ (keep=  &id &time pC_d  ) p = pC_d;
          title "Model for denominator of stabilized weights for &cens";
          run;

      
          %let modelCd = _modelCd_;

          proc logistic data = _cero_   outest = statcn(keep = _STATUS_ )    ;
                
          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;          
          &_where_ ;
          model &Cens = %if %bquote(&A_censoring)^= %then &A_censoring ; %else &A ; &covCn  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
          freq numberhits ;
          output out = _modelCn_ (keep=  &id &time pC_n   ) p = pC_n;
          title "Model for numerator of stabilized weights for &cens";
          run;
         
          %let modelCn = _modelCn_;
          %let numC = numC;
      
          
     %end;

    %if %bquote(&Cens2)^= %then %do;

    

          %let _where_ = ;
          %if %bquote(&eligible_cens2) ^= %then %let _where_ = &eligible_cens2 = 1 ;
          
          %if %bquote(&_where_) ^= %then %let _where_ = where ( &cens2 ne . ) AND ( &_where_)  ;
          %else %let _where_ = where &cens2 ne . ;
           
          proc logistic data = _cero_   outest = statc2d(keep = _STATUS_ )   ;              
          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
          &_where_ ;        
          model &Cens2 =  %if %bquote(&A_censoring)^= %then &A_censoring ; %else &A ; &covC2d  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
          freq numberhits ;
          output out = _modelC2d_ (keep=  &id &time pC2_d  ) p = pC2_d;
          title "Model for denominator of stabilized weights for &cens2";
          run;

    
          %let modelC2d = _modelC2d_;
   
          proc logistic data = _cero_   outest = statc2n(keep = _STATUS_ )     ;                
          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
          &_where_ ; 
          model &Cens2 = %if %bquote(&A_censoring)^= %then &A_censoring ; %else &A ;  &covC2n  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
          freq numberhits;
          output out = _modelC2n_ (keep=  &id &time pC2_n   ) p = pC2_n;
          title "Model for numerator of stabilized weights for &cens2";
          run;
     
          %let modelC2n = _modelC2n_;
       

     %end;

    %if %bquote(&Cens3)^= %then %do;

       
          %let _where_ = ;
          %if %bquote(&eligible_cens3) ^= %then %let _where_ = &eligible_cens3 = 1 ;
          
          %if %bquote(&_where_) ^= %then %let _where_ = where ( &cens3 ne . ) AND ( &_where_)  ;
          %else %let _where_ = where &cens3 ne . ;
           

          proc logistic data = _cero_  outest = statc3d(keep = _STATUS_ )     ;
          
          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
          &_where_ ;
      
          model &Cens3 =  %if %bquote(&A_censoring)^= %then &A_censoring ; %else &A ; &covC3d  %if %bquote(&logistic_options)^= %then / &logistic_options ;  ;
          freq numberhits ;
          output out = _modelC3d_ (keep=  &id &time pC3_d  ) p = pC3_d;
          title "Model for denominator of stabilized weights for &cens3";
          run;
      
        %let modelC3d = _modelC3d_;

          proc logistic data = _cero_  (where = (&cens3 ne .))  outest = statc3n(keep = _STATUS_ )    ;         
          ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
          &_where_ ;   
          model &Cens3 =  %if %bquote(&A_censoring)^= %then &A_censoring ; %else &A ; &covC3n  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
          freq numberhits;
          output out = _modelC3n_ (keep=  &id &time pC3_n   ) p = pC3_n;
          title "Model for numerator of stabilized weights for &cens3";
          run;
      
          %let modelC3n = _modelC3n_;
       

     %end;


     %if &use_natural_course = 1 %then %do;

         %if %bquote(&Cens)^= %then %do;

        
               %let _where_ = ;
               %if %bquote(&eligible_cens) ^= %then %let _where_ = &eligible_cens = 1 ;
           
               %if %bquote(&_where_) ^= %then %let _where_ = where ( &cens ne . ) AND ( &_where_)  ;
               %else %let _where_ = where &cens ne . ;
           

               proc logistic data =  _cero_    outest = statcdnc(keep = _STATUS_ )     ;
                
               ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
            
               &_where_ ;
               model &Cens = &A_censoring   &covCd  %if %bquote(&logistic_options)^= %then / &logistic_options ;  ;
               freq numberhits ;
               output out = _modelCdnc_ (keep=  &id &time pCnc_d  ) p = pCnc_d;
               title "Model for denominator of stabilized weights for &cens under natural course";
               run;

      
               %let modelCdnc = _modelCdnc_;

               proc logistic data = _cero_   outest = statcnnc(keep = _STATUS_ )     ;
                
               ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;          
               &_where_ ;
               model &Cens =  &covCn  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
               freq numberhits ;
               output out = _modelCnnc_ (keep=  &id &time pCnc_n   ) p = pCnc_n;
               title "Model for numerator of stabilized weights for &cens under natural course.";
               run;
         
               %let modelCnnc = _modelCnnc_;
          
          %end;

          %if %bquote(&Cens2)^= %then %do;

               %let _where_ = ;
               %if %bquote(&eligible_cens2) ^= %then %let _where_ = &eligible_cens2 = 1 ;
          
               %if %bquote(&_where_) ^= %then %let _where_ = where ( &cens2 ne . ) AND ( &_where_)  ;
               %else %let _where_ = where &cens2 ne . ;
           
               proc logistic data = _cero_   outest = statc2d(keep = _STATUS_ )    ;              
               ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
               &_where_ ;        
               model &Cens2 = &A_censoring   &covC2d  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
               freq numberhits ;
               output out = _modelC2dnc_ (keep=  &id &time pC2nc_d  ) p = pC2nc_d;
               title "Model for denominator of stabilized weights for &cens2 under the natural course.";
               run;

               %let modelC2dnc = _modelC2dnc_;
   
               proc logistic data = _cero_   outest = statc2n(keep = _STATUS_ )    ;                
               ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
               &_where_ ;         
               model &Cens2 =   &covC2n  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
               freq numberhits;
               output out = _modelC2nnc_ (keep=  &id &time pC2nc_n   ) p = pC2nc_n;
               title "Model for numerator of stabilized weights for &cens2 under the natural course.";
               run;
     
               %let modelC2nnc = _modelC2nnc_;
       
          %end;

          %if %bquote(&Cens3)^= %then %do;
       
               %let _where_ = ;
               %if %bquote(&eligible_cens3) ^= %then %let _where_ = &eligible_cens3 = 1 ;
          
              %if %bquote(&_where_) ^= %then %let _where_ = where ( &cens3 ne . ) AND ( &_where_)  ;
              %else %let _where_ = where &cens3 ne . ;
           
               proc logistic data = _cero_  outest = statc3d(keep = _STATUS_ )    ;          
               ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
               &_where_ ;
               model &Cens3 = &A_censoring   &covC3d  %if %bquote(&logistic_options)^= %then / &logistic_options ;  ;
               freq numberhits ;
               output out = _modelC3dnc_ (keep=  &id &time pC3nc_d  ) p = pC3nc_d;
               title "Model for denominator of stabilized weights for &cens3 under the natural course";
               run;
      
               %let modelC3dnc = _modelC3dnc_;

               proc logistic data = _cero_  (where = (&cens3 ne .))  outest = statc3n(keep = _STATUS_ )     ;         
               ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
               &_where_ ;    
               model &Cens3 =   &covC3n  %if %bquote(&logistic_options)^= %then / &logistic_options ; ;
               freq numberhits;
               output out = _modelC3nnc_ (keep=  &id &time pC3nc_n   ) p = pC3nc_n;
               title "Model for numerator of stabilized weights for &cens3 under the natural course.";
               run;
      
               %let modelC3nnc = _modelC3nnc_;
       
          %end;
     %end;



     title;

     data  _cero_ ;
     merge _cero_ &modelAd &modelAn &modelCd  &modelCn &modelC2d  &modelC2n &modelC3d  &modelC3n 
                &modelCdnc &modelCnnc   &modelC2dnc &modelC2nnc  &modelC3dnc &modelC3nnc 
                end = end;
     by    &id &time;
     if first.&id then do; 
          numA=1; denA=1;  numC=1; denC=1; numC2 = 1; denC2 = 1; numC3 = 1; denC3 = 1;
          %if &msm_var = 1 %then %do;
               den = 1;
               num = 1;
          %end;
                     
          _occasion_ = -1;
     end;
     retain  numA denA   numC denC numC2 denC2 numC3 denC3 ;
     /*  we have two cases
           1) msm_var = 0 and macro %msm_var_initialization is not run and the proc logisitc or genmod will not use the observations
                      when the outcome is missing
           2) msm_var = 1 and the macro msm_var_initialization is run and the we will need to include the scores for the last observation
                       when the observation may have contributed to the censoring and treatment models, but not the final model. In the 
                       initialization macro the observation will be removed once the scores have been calculated.
       

     ***/
   
     %if %eval(&treatment_weights) = 0 %then %do;
          pA_d = 1;
          pA_n = 1 ;                
     %end;
     %else %do;
          if pA_d = . then pA_d = 1;
          if pA_n = . then pA_n = 1 ;
     %end;

     prob_holder_d = pA_d ;
     prob_holder_n = pA_n ;       

     %if %eval(&treatment_weights) = 1 %then %do;
          %if %bquote(&eligible)=  %then %do;
               if &A = 0 then do;
                    prob_holder_d = 1-pA_d;
                    prob_holder_n = 1-pA_n;
               end;
          %end;
          %else %if %bquote(&eligible)^= %then %do;
               if &eligible = 1 AND &A = 0 then do;
                    prob_holder_d = 1-pA_d;
                    prob_holder_n = 1-pA_n;
               end;
               if &eligible = 0 then do;
                    prob_holder_d = 1;
                    prob_holder_n = 1;
               end;
          %end;
     %end;
    
    numA = numA * prob_holder_n;
    denA = denA * prob_holder_d ;



     %if %bquote(&Cens)^= %then %do;
          %if %bquote(&eligible_cens)^= %then %do;
               if &eligible_cens = 0 then do;
                    pC_n = 1; /* overwrite any value that may have come from the proc logistic */
                    pC_d = 1;
               end;
          %end;
          if pC_d = . then pC_d = 1;
          if pC_n = . then pC_n = 1;

           
            denC = denC * pC_d;
            numC = numC * pC_n;
     %end;

     %if %bquote(&Cens2)^= %then %do;
           %if %bquote(&eligible_cens2)^= %then %do;
               if &eligible_cens2 = 0 then do;
                    pC2_n = 1;
                    pC2_d = 1;
               end;
          %end;
          if pC2_d = . then pC2_d = 1;
          if pC2_n = . then pC2_n = 1;
            
          numC2 = numC2* pC2_n;
          denC2 = denC2* pC2_d;
     %end;

     %if %bquote(&Cens3)^= %then %do;
          %if %bquote(&eligible_cens3)^= %then %do;
               if &eligible_cens3 = 0 then do;
                    pC3_n = 1;
                    pC3_d = 1;
               end;
          %end;
          if pC3_d = . then pC3_d = 1;
          if pC3_n = . then pC3_n = 1;
            
          numC3 = numC3* pC3_n;
          denC3 = denC3* pC3_d;
     %end;

  

     swnum = numA  ;
     %if %bquote(&Cens)^= %then swnum = swnum * numC ;;
     %if %bquote(&cens2)^= %then swnum = swnum * numC2 ;;
     %if %bquote(&cens3)^= %then swnum = swnum * numC3 ;;

     swdenum = denA   ;
     %if %bquote(&Cens)^= %then swdenum = swdenum * denC ;;
     %if %bquote(&cens2)^= %then swdenum = swdenum * denC2 ;;
     %if %bquote(&cens3)^= %then swdenum = swdenum * denC3 ;;




     %if &use_stabilized_wts = 1 %then %do ;
          sw = swnum / swdenum ;      
     %end;
     %else %do;
          sw = 1.0 / swdenum ;
     %end;

     wtA = numA / denA ;
     %if %bquote(&cens)^= %then wtC = numC / denC ;;
     %if %bquote(&cens2)^= %then wtC2 = numC2 / denC2 ;;
     %if %bquote(&cens3)^= %then wtC3 = numC3 / denC3 ;;

     %if &use_natural_course = 1 %then %do;
          retain natural_wt ;
          
          if first.&id then do ;
             natural_wt = 1.0 ;
          end;

          %if %bquote(&Cens)^= %then %do;
             %if %bquote(&eligible_cens)^= %then %do;
               if &eligible_cens = 0 then do;
                    pCnc_n = 1; /* overwrite any value that may have come from the proc logistic */
                    pCnc_d = 1;
               end;
             %end;
          
             if pCnc_d = . then pCnc_d = 1;
             if pCnc_n = . then pCnc_n = 1;

               wtCnc = pCnc_n / pCnc_d ;
          %end ;

          %if %bquote(&Cens2)^= %then %do;
             %if %bquote(&eligible_cens2)^= %then %do;
                if &eligible_cens2 = 0 then do;
                    pC2nc_n = 1; /* overwrite any value that may have come from the proc logistic */
                    pCn2c_d = 1;
                end;
             %end;
          
       
             if pC2nc_d = . then pC2nc_d = 1;
             if pC2nc_n = . then pC2nc_n = 1;
             wtC2nc = pC2nc_n / pC2nc_d ;
      %end;
     
          %if %bquote(&Cens3)^= %then %do;
             %if %bquote(&eligible_cens3)^= %then %do;
               if &eligible_cens3 = 0 then do;
                    pC3nc_n = 1; /* overwrite any value that may have come from the proc logistic */
                    pC3nc_d = 1;
               end;
             %end;
          
             if pC3nc_d = . then pC3nc_d = 1;
             if pC3nc_n = . then pC3nc_n = 1;
              wtC3nc = pC3nc_n / pC3nc_d ;
          %end;

        
                
         %if %bquote(&cens)^=  %then natural_wt = natural_wt * wtCnc ;;
         %if %bquote(&cens2)^= %then natural_wt = natural_wt * wtC2nc ;;
         %if %bquote(&cens3)^= %then natural_wt = natural_wt * wtC3nc  ;;



     %end;

     

     _occasion_ + 1;
     if &outcMSM = 1 then _events_used + 1 ;
     if end then  call symput('events_used',_events_used) ;

     keep numberhits &id &idname &covMSMbh &classMSMbh &covMSMbv &classMSMbv &inter_MSM &outcMSM &A &AMSM  sw _occasion_
     &time &_interactions &_interactions_bh 
     
     %if &msm_var = 1 %then %do;
        &covAd &covCd   &cens &cens2 &cens3 &eligible &eligible_cens &eligible_cens2 &eligible_cens3
     %end;
     %if %eval(&treatment_weights) = 1 %then wtA  pA_d pA_n;
     %if %bquote(&cens)^= %then wtC pC_n pC_d ;
     %if %bquote(&cens2)^= %then wtC2 pC2_n pC2_d ;
     %if %bquote(&cens3)^= %then wtC3 pC3_n pC3_d ;
     %if &use_natural_course = 1 %then natural_wt ;
     %if %bquote(&competing_risk) ^= %then &competing_risk ;
     
     ;  
     run;

     %if  %eval(&bsample) = 0 %then %do;
          data _weights;
          set _cero_ (keep = &idname &time sw numberhits
                   %if %eval(&treatment_weights) = 1 %then wtA pA_d pA_n ;
                   %if %bquote(&cens)^= %then wtC pC_n pC_d  ;
                   %if %bquote(&cens2)^= %then wtC2 pC2_n pC2_d ;
                   %if %bquote(&cens3)^= %then wtC3 pC3_n pC3_d ;
                   %if &use_natural_course = 1 %then natural_wt ;
                   ) ;             
          run;
     %end;

    
 
     %if %eval(&no_listing) = 0 %then %do;
          %let _weight_list = &A &cens &cens2 &cens3 ;
     
          %let _n_weight = %numargs(&_weight_list);
          %let _strt = 1;      
          %do _wght = &_strt %to %eval(&_n_weight);
          
               %if &_wght = 1    %then %do;
                    %if &treatment_weights = 1 %then %do;
                        %let var = wtA;
                        %let nametmp = &A;
                    %end;
                    %else %do;
                        %let nametmp = %scan(&_weight_list, %eval(&_wght),%str(' '));
                        %let var = wtC ;
                    %end; 
               %end;
               %else %do;
                    %if &treatment_weights = 1 %then %do;
                        %let nametmp = %scan(&_weight_list, %eval(&_wght),%str(' '));
                        %let ii = %eval(&_wght - 1);                     
                        %if &ii = 1 %then %let var = wtC;
                        %else %let var = wtC&ii;
                    %end;
                    %else %do;
                       %let nametmp = %scan(&_weight_list, %eval(&_wght),%str(' '));
                       %let var = wtc&_wght ;
                    %end;
               %end;
               proc univariate data = _cero_ (keep = &id  &var &outcMSM numberhits )
                                        %if %eval(&bsample) > 0 %then noprint ; mu0 = 1  ;
               title  "Distribution of stabilized weights contribution from  &nametmp "; 
               id &id ;
               where &outcMSM ^= . ;
               var  &var ;
               freq numberhits ;
               run;
               title ;
          %end;
     %end;


     proc univariate data = _cero_ (keep = &id sw &outcMSM numberhits )  mu0 = 1 ;
     title  'Distribution of stabilized weights'; 
     where &outcMSM ^= . ;
     id &id ;
     var sw   ;
     freq numberhits ;
     output out = _weight_levels   pctlpts= &truncate_weights %eval(100 - &truncate_weights) pctlpre=level_  
                     pctlname=lower upper ;

 
     run;

     %if %eval(&truncate_weights) > 0 %then %do;
       
          %if &user_defined_limits = 0 %then %do;
               data _null_;
               set _weight_levels ;
               call symput('lower_limit',trim(left(level_lower)));
               call symput('upper_limit',trim(left(level_upper)));
               run;
          %end;

          data _cero_ ;
          set _cero_  ;
                              
          if sw < &lower_limit then sw = &lower_limit;
          if sw > &upper_limit then sw = &upper_limit;           

          %if &use_natural_course = 1 %then %do;
              if natural_wt < &lower_limit then natural_wt = &lower_limit ;
              if natural_wt > &upper_limit then natural_wt = &upper_limit ;
          %end;

          drop wtA  %if %bquote(&cens)^= %then wtC;
                    %if %bquote(&cens2)^= %then wtC2;
                    %if %bquote(&cens3)^= %then wtC3 ; ;
          run;         

          proc univariate data = _cero_ mu0 = 1.0 ;
          title "Analysis of truncated weights with &lower_limit <= sw <= &upper_limit";
           where &outcMSM ^= . ;
           var sw ;
           freq numberhits ;
           run;
           title ;
     %end;

     data _status_;
     set 
      %if %eval(&treatment_weights) = 1 %then statAd statAn;
      %if %bquote(&cens)^= %then statCd statCn ;
      %if %bquote(&cens2)^= %then statC2d statCn ;
      %if %bquote(&cens3)^= %then statC3d statC3n ;
      %if &use_natural_course = 1 %then %do;
          %if %bquote(&cens)^= %then statCdnc statCnnc ;
          %if %bquote(&cens2)^= %then statC2dnc statCnnc ;
          %if %bquote(&cens3)^= %then statC3dnc statC3nnc ;
      %end;
      end = _end_ ;
      if _n_ = 1 then good = 1;
      holder = compress(substr(_STATUS_,1,1));
      if holder ^= '0' then good = 0;
      retain good;
      if _end_ = 1 then call symput('weight_good',good);
      run;

      %if %eval(&weight_good) = 0 %then %do;
          %let just_models_tmp = &just_models ; /* original value, do not want to change if was originally 1 */
          %let just_models = 1 ;
          %put Setting just_models = 1 due to problems with weight models ;
          proc print data = _status_;
          run;
          %if &override_warnings = 1 %then %do;
              %put resetting just_models = &just_models_tmp due to override_warnings = 1 ;
              %let just_models = &just_models_tmp ;
          %end;

     %end;



    %if &debug = 0 %then %do;
          proc datasets library=work nolist;
          delete &modelAd &modelAn   
                 &modelCd &modelCn  &modelC2d &modelC2n &modelC3d &modelC3n
                 &modelCdnc &modelCnnc  &modelC2dnc &modelC2nnc &modelC3dnc &modelC3nnc
                  _weight_levels
                  %if %eval(&treatment_weights) = 1 %then statAd statAn;
                  %if %bquote(&cens)^= %then statCd statCn ;
                  %if %bquote(&cens2)^= %then statC2d statCn ;
                  %if %bquote(&cens3)^= %then statC3d statC3n ;
                  %if %bquote(&competing_risk)^= %then %do;
                     %if %bquote(&cens)^= %then statCdnc statCnnc ;
                     %if %bquote(&cens2)^= %then statC2dnc statCnnc ;
                     %if %bquote(&cens3)^= %then statC3dnc statC3nnc ;
                  %end;
                  %if &treatment_weights = 1 or %bquote(&cens)^= OR   
                      %bquote(&cens2)^= OR %bquote(&cens3)^= %then _status_ ;
                 
                  ;     
          quit;
   %end;
%mend calculate_weights;

%macro initialize_msm_var;

     %if %eval(&treatment_weights) = 1 %then %let n_covAd = %numargs(&covAd);  
     %else %let n_covAd = 0;

     
     %let covMSM = &covMSMbh &covMSMbv;
     %let n_covMSM = %numargs(&covMSM);
    
     %if   %eval(&treatment_weights) = 1 %then %let scoreAd = scoreA_int;  
     %else %let scoreAd = ;
     %do i = 1 %to &n_covAd;
          %let vartmp = %scan(&covAd, %eval(&i));
          %let scoreAd = &scoreAd scoreA_&vartmp ;
     %end;

     %let scores = &scoreAd ;
     %if %eval(&treatment_weights) = 1 %then %let n_scores = %eval(&n_covAd + 1);
     %else %let n_scores = 0 ;

     %if %bquote(&Cens)^= %then %do;
          %let n_covCd = %numargs(&covCd);  
          %let scoreCd = scoreC_int;
          %do i = 1 %to &n_covCd;
               %let vartmp = %scan(&covCd, %eval(&i));
               %let scoreCd = &scoreCd scoreC_&vartmp ;
          %end;
          %let scores = &scores &scoreCd ;
          %let n_scores = %eval(&n_scores + &n_covCd + 1 );
     %end;
     %if %bquote(&Cens2)^= %then %do;
          %let n_covC2= %numargs(&covC2d);  
          %let scoreC2 = scoreC2_int;
          %do i = 1 %to &n_covC2;
               %let vartmp = %scan(&covC2d, %eval(&i));
               %let scoreC2 = &scoreC2 scoreC2_&vartmp ;
          %end;
          %let scores = &scores &scoreC2 ;
          %let n_scores = %eval(&n_scores + &n_covC2 + 1 );
     %end;
     %if %bquote(&Cens3)^= %then %do;
          %let n_covC3= %numargs(&covC3d); 
          %let scoreC3 = scoreC3_int;
          %do i = 1 %to &n_covC3;
               %let vartmp = %scan(&covC3d, %eval(&i));
               %let scoreC3 = &scoreC3 scoreC3_&vartmp ;
          %end;
          %let scores = &scores &scoreC3 ;
          %let n_scores = %eval(&n_scores + &n_covC3 + 1 );
     %end;
 


     data _for_follow;
     set  _cero_ end = _end_;
     by &id   &time;
     _keep = 0;
     if first.&id then _subjects + 1;   
     if last.&id then do;
          followup = _occasion_ ;
          if &outcMSM = 1 then _dies = 1;
          else _dies = 0;
          _keep = 1;
     end;
     if _end_ then call symput('individuals',_subjects);
     if _keep = 0 then delete;
     keep &id followup _dies;
     run;

     data _cero_;
     merge _cero_ _for_follow;
     by &id;
     derien = 1;
     run;



     data nuisance; 
     set  _cero_; 
     by &id &time; 
     retain  &scores ; 
        
     if first.&id  then do;                          
          %do _i = 1 %to  &n_scores;
               %let vartmp = %scan(&scores, %eval(&_i));
               &vartmp = 0;
          %end;
     end;
   /************************
         need to take into account the possiblity that the eligibility criterion is not monotonic that A_n is not
         necessarily once treated always treated so A=1 and pA_d = 1 when &eligible = 0 and the value of _tmp = 0. If
         &eligible = 0 and A_n = 0 and pA_d = 1 then _tmp is not 0 and there is a contribution to the scores.
      *****************/ 
     %if %eval(&treatment_weights) = 1 %then _tmp = &A - pA_d ;;

     %if %bquote(&eligible)^=  %then if &eligible = 0 then _tmp = 0 ;;

     %if %bquote(&Cens)^= %then %do;
          if &cens = . then _tmp = 0;
     %end;

     %if %eval(&treatment_weights) = 1 %then scoreA_int = _tmp + scoreA_int;;

     
     %do _i = 1 %to %eval(&n_covAd);
          %let vartmp = %scan(&covAd, %eval(&_i));
          _tmp = &vartmp * (&A - pA_d) ;
          %if %bquote(&Cens)^= %then %do;
               if &cens = . then _tmp = 0;
          %end;

         %if %bquote(&Cens2)^= %then %do;
               if &cens2 = . then _tmp = 0;
         %end;

         %if %bquote(&Cens3)^= %then %do;
               if &cens3 = . then _tmp = 0;
         %end;
         scoreA_&vartmp =  _tmp + scoreA_&vartmp ;
     %end;

     %if %bquote(&Cens)^= %then %do; 
          if &cens > . then do;
                    %if %bquote(&eligible_cens)^= %then if &eligible_cens = 1 then  do; ;                              
                              scoreC_int = (&Cens - 1 + pC_d)+scoreC_int;   
                              %do _i = 1 %to &n_covcd;
                                   %let vartmp = %scan(&covcd ,%eval(&_i));
                                   scoreC_&vartmp = &vartmp * (&Cens - 1 + pC_d ) + scoreC_&vartmp;
                              %end;
                    %if %bquote(&eligible_cens)^= %then end;;    
                     
                   
           end;
     %end; 
     %if %bquote(&Cens2)^= %then %do; 
          if &cens2 > . then do;
                    %if %bquote(&eligible_cens2)^= %then if &eligible_cens2 = 1 then  do; ;
                         scoreC2_int = (&Cens2 - 1 + pC2_d)+scoreC2_int;   
                         %do _i = 1 %to &n_covc2;
                              %let vartmp = %scan(&covC2d ,%eval(&_i));
                              scoreC2_&vartmp = &vartmp * (&Cens2 - 1 + pC2_d ) + scoreC2_&vartmp;
                         %end;
                    %if %bquote(&eligible_cens2)^= %then end;; 
           end;
     %end; 

     %if %bquote(&Cens3)^= %then %do; 
          if &cens3 > . then do;
                    %if %bquote(&eligible_cens3)^= %then if &eligible_cens3 = 1 then  do; ;
                         scoreC3_int = (&Cens3 - 1 + pC3_d)+scoreC3_int;   
                         %do _i = 1 %to &n_covc3;
                              %let vartmp = %scan(&covC3d ,%eval(&_i));
                              scoreC3_&vartmp = &vartmp * (&Cens3 - 1 + pC3_d ) + scoreC3_&vartmp;
                         %end;
                    %if %bquote(&eligible_cens3)^= %then end;; 
           end;
     %end; 

     if _occasion_ = followup then do; 
          keep &id  &scores ; 
          output;
     end;
     run;



    data _cero_; 
     merge _cero_ nuisance; 
     by &id ;
    /** change this to be delete observation only if &outcmsm = . , that is whenever the subject 
        is censored **/
     if &outcMSM = . then delete ;
     rien = 1; 
    run;

%mend initialize_msm_var;

%macro calculate_msm_var ;

     %local i  npred word new_u ;
     proc freq data = _cero_  noprint;
     table &id / out = tcount (drop = percent);        
     run;

    %if &debug = 0 %then %do;
         proc datasets library = work nolist;
               delete    _for_follow  ; 
         quit;
     %end; 
 

    
    %let pred = one &amsm &_interactions &covMSMbh &_interactions_bh &covMSMbv ;
    %let npred = %numargs(&pred) ;
  
    %let new_u = ;
    %let dE_dbeta =  ;
    %do i = 1 %to &npred ;
         %let new_u = &new_u new_u&i ;
         %let dE_dbeta = &dE_dbeta dE_dbeta&i ;
    %end;

    data _formsmvar_ ;
    merge _formsmvar_ _cero_ (in = _keep_ keep = &id &time &amsm &_interactions &covMSMbh &_interactions_bh  
                                   &covMSMbv &outcmsm sw );
    by &id &time ;
    if _keep_ ;

    array dE_dbeta dE_dbeta1 - dE_dbeta&npred ;
    array pred &pred ;

    u = sw *(&outcMSM - hazard);
    htmp = hazard * (1.0 - hazard);
   
    one = 1 ;
    /** calculate new_u **/
    %do i = 1 %to &npred ;
       %let word = %qscan(&pred, %eval(&i) ,%str( ));
       new_u&i = u * &word ;
    %end;

    do i = 1 to %eval(&npred);
         dE_dbeta[i] = -1 * sw * htmp * pred[i] ;
    end;

    run;

    /* sum the values of new_u to create NUvec by id */
    proc means data = _formsmvar_ (keep = &id &new_u ) noprint ;
    var &new_u ;
    by &id ;
    output out = NUvec(drop = _TYPE_ _FREQ_) sum(&new_u) =  ;
    run;

    proc corr data = NUvec sscp outp = lambda0 (keep = &new_u _TYPE_ _NAME_ 
               where = (_TYPE_ = "SSCP" & _NAME_^="Intercept")) noprint;
    var &new_u ;
    run;

    proc means data = _formsmvar_ (keep = &dE_dbeta) noprint;
    var &dE_dbeta ;
    output out = new_deriv sum( &dE_dbeta ) = ;
    run;

    %do i = 2 %to &npred ;
          %let word = %qscan(&pred,%eval(&i),%str( ));
          proc means data = _formsmvar_ (keep = &word &dE_dbeta) noprint ;
          var &dE_dbeta ;
          weight &word ;
          output out = deriv_tmp sum( &dE_dbeta ) = ;
          run;

          proc append base = new_deriv data = deriv_tmp ;
          run;
     %end;

    /* for gamma0 we need the outer product of two different vectors. We can use proc corr with
       all the variables and pull out the desired submatrix using iml. This should be the upper right-hand
       submatrix */

    data _forgamma0_ ;
    merge NUvec _cero_ (keep = &id _occasion_ &scores where = (_occasion_ = 0));
    by &id ;
    run;

    proc corr data = _forgamma0_ sscp outp = gamma0 (keep = &new_u &scores  _TYPE_ _NAME_ 
               where = (_TYPE_ = "SSCP" & _NAME_^="Intercept")) noprint;
    var &new_u &scores ;
    run;

    proc corr data = nuisance sscp outp = omega0 (keep = &scores  _TYPE_ _NAME_ 
               where = (_TYPE_ = "SSCP" & _NAME_^="Intercept")) noprint;
    var &scores ;
    run;

    data omega0 ;
    set omega0(drop = _TYPE_ _NAME_);
    run;

    data gamma0 ;
    set gamma0 (drop = _TYPE_ _NAME_ );
    run;


    %let nscores = %numargs(&scores);

  
     proc iml ;

     use lambda0 ;
     read all var { &new_u} into lambda0 ;

     use gamma0 ;
     read all into pregamma0 ;

     use omega0 ;
     read all into omega0 ;

     gamma0 = pregamma0[1:&npred , (%eval(&npred+1):%eval(&npred + &nscores))] ;

     omega0_eigen = eigval(omega0);
     omega0_det = det(omega0);
     omega0_cond =  max(abs(omega0_eigen))/min(abs(omega0_eigen)) ;
     print omega0_det omega0_cond ;
     file log ;
     put @1 "det(omega0) = " omega0_det +5 "cond(omega0)= " omega0_cond ;
     omega0 = &inv_type (omega0) ;

     meat = lambda0 - gamma0 * omega0 * t(gamma0) ;

     use new_deriv ;
     read all var { &dE_dbeta} into new_deriv ;

     new_deriv_eigen = eigval(new_deriv);
     new_deriv_det = det(new_deriv);
     new_deriv_cond =  max(abs(new_deriv_eigen))/min(abs(new_deriv_eigen)) ;
     put @1 "det(new_deriv) = " new_deriv_det +5 "cond(new_deriv)= " new_deriv_cond ;
     ginv_deriv = &inv_type (new_deriv ) ;

     varbetas = ginv_deriv * meat * t(ginv_deriv );

     use emppest_sw;
     read all var {estimate} into betas;
     create varbetas_sw from varbetas  ;
     append from varbetas;
     sterrB=varbetas[2,2]##0.5;
     onec = "a" ;
     test0 = repeat(onec,1,100);
     nametmp = rowcatc(test0);
     create _kk_ var{nametmp betasw std lb ub z p_value};
     %do i = 1 %to &n_betas;

         sterrb = varbetas[&i,&i]##0.5 ;
         %let vartmp = %scan(&cov_betas, %eval(&i));
         nametmp =  "&vartmp" ;
         betasw = betas[&i];
         lb = betas[&i] - 1.96 * sterrb;
         ub = betas[&i] + 1.96 * sterrb;
         std = sterrb;
         if sterrb = 0 |  sterrb = . then do ;
             z = . ;
             p_value = . ;
         end;
         else do;
             z = betasw / sterrb;
             p_value = 2*(1.0 - probnorm(abs(z)));
         end;
             
         append;

     %end;

     quit;

    data _kk_;
    set  _kk_ (rename = (nametmp = NAME betasw = estimate))  ;
    bsample = &bsample; 
    _row_ = _N_ ;         
    run;

    data varbetas_sw ;
    set varbetas_sw ;
    _row_ = _N_ ;
    run;

    data eee_cont ;
    set eee_cont ;
    _row_ = _N_ ;
    keep _row_ Notes ;
    run;


    data varbetas_sw;
    length name $20;
    merge varbetas_sw eee_cont ;
    by _row_ ;
    rename  
        %do i = 1 %to %eval(&n_betas );
             %let vartmp = %scan(&cov_betas,%eval(&i));
             col&i = &vartmp 
         %end;;
    array names[&n_betas] $20 names1-names&n_betas ;

    %do i = 1 %to %eval(&n_betas );
         %let vartmp = %scan(&cov_betas, %eval(&i));
         names&i =  "&vartmp " ;
    %end;
    name = names[_n_ ];  
     
    drop names1-names&n_betas;
    run;

    proc datasets library = work nolist;
    delete tcount gamma0 lambda0 new_deriv deriv_tmp 
           nuvec omega0 _forgamma0_ _formsmvar_;
    quit;

%mend calculate_msm_var;
%macro curves ;

/** three possible sets of curves 
     1) a) always treated vs never treated using original weights (uses empest data set )
        b) always treated vs never treated using original weights along with a competing risk model (empest_cr)
        c) with modeling for censoring due to lost to follow-up
     2) natural course model with no treatment variables (uses empest_nc)
 ****/


            
          %let numbh = %numargs(&covMSMbh) ;
          %let numbh_c = %numargs(&classMSMbh) ;
          %let numAMSM = %numargs(&AMSM) ;
          %let numAMSM_cats = %numargs(&amsm_knots);
          
                        
           proc sql  noprint  ;
              select variable into :varname_sw separated by ' ' from pe_sw ;
              select estimate into :estimate_sw separated by ' ' from pe_sw ;
              %if &use_natural_course = 1 %then %do;
                   select variable into :varname_nc separated by ' '  from pe_nc ;                   
                   select estimate into :estimate_nc separated by ' ' from pe_nc ;
                   %if %bquote(&competing_risk)^= %then %do;
                      select estimate into :estimate_nc_cr separated by ' ' from pe_nc_cr ;
                   %end;
              %end;
              %if %bquote(&competing_risk)^= %then %do;
                  select estimate into :estimate_cr separated by ' ' from pe_cr ;
              %end;
           quit;

                   
           %let numsw = %numargs(&estimate_sw) ;
           %if &use_natural_course = 1 %then %let numnc = %numargs(&estimate_nc);
 
           data curves ;
           set _cero_ (keep = &id &time  &covMSMbh &covMSMbv   numberhits );
           by &id &time ;
           if first.&id;
           array est_sw{&numsw} _TEMPORARY_ ( &estimate_sw ) ;
           array var_sw &varname_sw ;

           array ci_treated ci_treated&time_start-ci_treated&time_end ;
           array ci_nontreated ci_nontreated&time_start-ci_nontreated&time_end ;
           array km_treated km_treated&time_start-km_treated&time_end ;
           array km_nontreated km_nontreated&time_start - km_nontreated&time_end ;

           %if &use_natural_course = 1 %then %do;
               array ci_nc ci_nc&time_start - ci_nc&time_end;
               array km_nc km_nc&time_start - km_nc&time_end;
               array est_nc{&numnc} _TEMPORARY_ ( &estimate_nc ) ;
               array var_nc &varname_nc ;
               %if %bquote(&competing_risk)^= %then %do;
                  array est_nc_cr {&numnc} _TEMPORARY_ ( &estimate_nc_cr );
                  array var_nc_cr &varname_nc ;                  
               %end;
           %end;
           %if %bquote(&competing_risk)^= %then %do;
               array est_cr {&numsw} _TEMPORARY_ ( &estimate_cr );
               array var_cr &varname_sw ;
              
           %end;
           Intercept = 1.0 ;
           * initialize variables for calculating linear predictor for just baseline
              variables. This assumes that AMSM = 0 when treatment is also 0 for all amsm_type ;
           %do ii = 1 %to &numAMSM ;
              * treatment variables ;
              %let word = %scan(&AMSM,&ii,%str( ));
               &word = 0 ; 
           %end;


          /* interactions with amsm */
           %if &ndef_amsm > 0 %then %do;
                * any interactions with treatment variables ;
                %do ninter = 1 %to &ndef_amsm ;
                     %let word =  %scan(&amsm_inter_def,%eval(&ninter),%str(:)) ;
                     &word ;
                %end;
           %end;
           %do ii = 1 %to &numbh ;
                   * time variables ;
                   %let word = %qscan(&covMSMbh, &ii, %str( ));
                   &word = 0 ;
               %end ;

           /* interactions with covMSMbh */
           %if &ndef_time > 0 %then %do ;
              *interactions with time variables ;
              %do ninter = 1 %to &ndef_time ;
                  %let word =  %scan(&covbh_inter_def,%eval(&ninter),%str(:)) ;
                  &word ;          
              %end;
           %end;

           /* any user defined variables go here, only for initializing variables used in linear predictor : for example
            variables used in interactions but not used directly in weighted model. */

           %extra_variables_def ;
        
                 
           cuminc_treated = 0 ;
           cuminc_nontreated = 0 ;
           surv_treated = 1 ;
           surv_nontreated = 1 ;
           
           xbeta_base_sw = 0.0 ;
           do i = 1 to dim(est_sw);
               xbeta_base_sw = xbeta_base_sw + est_sw[i] * var_sw[i];
           end;

           %if &use_natural_course = 1 %then %do;
               surv_nc = 1.0 ;
               cuminc_nc = 0.0 ;                         
               xbeta_base_nc = 0 ;
               do i = 1 to dim(est_nc);
                   xbeta_base_nc = xbeta_base_nc + est_nc[i] * var_nc[i] ;
               end;
               %if %bquote(&competing_risk)^= %then %do;
                   surv_nc_cr = 1.0 ;
                   risk_nc_cr = 0.0 ;
                   xbeta_base_nc_cr = 0 ;
                   do i = 1 to dim(est_nc_cr);
                       xbeta_base_nc_cr = xbeta_base_nc_cr + est_nc_cr[i] * var_nc_cr[i];
                   end;
               %end;
           %end ;

           %if %bquote(&competing_risk)^= %then %do;
               surv_treated_cr    = 1.0 ;
               surv_nontreated_cr = 1.0 ;
               risk_treated_cr    = 0.0 ;
               risk_nontreated_cr = 0.0 ;
               xbeta_base_cr = 0 ;
               do i = 1 to dim(est_cr);
                   xbeta_base_cr = xbeta_base_cr + est_cr[i] * var_cr[i];
               end;
           %end;
           
           do &time = &time_start to &time_end by &time_step;
               
                  * calcuated time variables ;

                  %if &user_function_of_time = 1 %then %user_function_of_time ;
                  %else %do ;
                      %if &contMSMbh=1 and &use_splines = 1 %then %do;
                          %rcspline(&time,&knots );
                      %end;

                
         

                      %if &contMSMbh = 0 %then %do;
                           %do ii = 1 %to &numbh ;
                               %let word = %qscan(&covMSMbh, &ii, %str( ));
                               &word = 0 ;
                           
                               %if %upcase(&class_ref)=FIRST %then %do;
                                   %let _iq1 = %str(<);
                                   %let _iq2 = %str(<=);
                               %end;
                               %else %if %upcase(&class_ref)=LAST %then %do;
                                   %let _iq1 = %str(<=);
                                   %let _iq2 = %str(<) ;
                               %end;
                               %if &ii <= &numbh %then %do; 
                                  %let level0 = %scan(&time_knots,&ii,%str( ));
                                  %let level1 = %scan(&time_knots,%eval(&ii + 1),%str( ));
                                  if &level0 &_iq1 &time &_iq2 &level1 then &word = 1;
                               %end;
      
                           %end ;
                          
                     %end;
                 %end; 
                 * create any interaction terms with time variables ;
                  %if &ndef_time > 0 %then %do ;
                      %do ninter = 1 %to &ndef_time ;
                         %let word =  %scan(&covbh_inter_def,%eval(&ninter),%str(:)) ;
                         &word ;
          
                      %end;
                 %end;

                        
                
                 
                  /* ordering of variables is : intercept  amsm  [interactions with amsm ]  covMSMbh [interactions with covMSMbh ] covMSMbv */

                 %if &use_natural_course = 1 %then %do; /* no amsm variables */
                       * add time and interactions with time into natural course ;
                        xbeta_nc = xbeta_base_nc ;
                        %let timevars = %eval(&numbh + %numargs(&_interactions_bh)) ;
                        %if %bquote(&competing_risk)^= %then xbeta_nc_cr = xbeta_base_nc_cr ;;
                        do j = 2 to %eval(2+ &timevars /* &numbh */ -1);
                             xbeta_nc = xbeta_nc + est_nc[j] * var_nc[j] ;                       
                             %if %bquote(&competing_risk)^= %then xbeta_nc_cr = xbeta_nc_cr + est_nc_cr[j]*var_nc_cr[j];;
                        end; 

                      
                         
                                                                                                                                                          
                  %end;

                  /* ASSUMING THAT AMSM VARIABLES ARE 0 WHEN TREATMENT IS 0 (can change this for comparing two time
                     varying regimes (?) by using the amsm_type = 3 otptions and two separate sub-macros */


                  %if &amsm_type ^= 3 %then %do;
                       %do ii = 1 %to &numAMSM ;
                          %let word = %scan(&AMSM,&ii,%str( ));
                           &word = 0 ; 
                       %end;
                  %end;
                  %else %do;
                       %amsm_user_def_macro0 ;                       
                  %end;
                  
                  %if &ndef_amsm > 0 %then %do;

                     %do ninter = 1 %to &ndef_amsm ;
                         %let word =  %scan(&amsm_inter_def,%eval(&ninter),%str(:)) ;
                         &word ;
                     %end;
                  %end;



              
                 * for linear predictor for treated and never treated first add in time variables and interactions with time
                   to xbeta_base_sw   ;
                  

                  xbeta_sw = xbeta_base_sw ;

                  %let _startpos = %eval(1 + &numAMSM + %numargs(&_interactions) + 1);
                  %let _lastpos = %eval(&_startpos + &numbh + %numargs(&_interactions_bh) -1);

                  do j = &_startpos to &_lastpos ;
                        xbeta_sw = xbeta_sw + est_sw[j] * var_sw[j] ;
                  end;

                   
                                                                   
                  %if %bquote(&competing_risk)^= %then %do;
                       xbeta_cr = xbeta_base_cr ;                     
                       do j = &_startpos to &_lastpos ;
                            xbeta_cr = xbeta_cr + est_cr[j] * var_cr[j] ;
                       end;
                       

                 %end;
               
                 * for never treated, assume contribution is 0 so only need xbeta_sw 
                     except in possible case of amsm_type = 3 ;
                 
                  xbeta_sw0 = xbeta_sw ; 
                  %if %bquote(&competing_risk)^= %then xbeta_cr0 = xbeta_cr ;;

                  %if &amsm_type = 3 %then %do;
                        do j = 2 to %eval(1+&numAMSM + &ndef_amsm) ;
                            xbeta_sw0 = xbeta_sw0 + est_sw[j] * var_sw[j]; /* 7-2014 fixed loop for xbeta_sw to xbeta_sw0 */
                        end;
                        /* 7-2014 added following for general treatment type other than being constantly 0 */
                        %if %bquote(&competing_risk)^= %then %do;
                            do j = 2 to %eval(1+&numAMSM + &ndef_amsm) ;
                               xbeta_cr0 = xbeta_cr0 + est_cr[j] * var_cr[j]; 
                            end;
                        %end;
                  %end;

                 %if %bquote(&competing_risk)^= %then xbeta_cr0 = xbeta_cr ;;
                                           
                /* add in part for always treated */

                  * start linear predicter for always treated . only need to add in amsm variables ;
                 /* 7-2014 rwl changed to xbeta_sw  since for general case never-treated 
                         may change to form of delayed treatment so the contribution may be non-zero , do not want to include
                         these values in starting xbeta for treated. Only want time and baseline values. */
                  xbeta_sw1 = xbeta_sw ; 
                  %if %bquote(&competing_risk)^= %then xbeta_cr1 = xbeta_cr ;;
                  %if &amsm_type = 0 %then %do;
                      &amsm = 1 ;
                      xbeta_sw1 = xbeta_sw1 + var_sw[2] * est_sw[2] ;
                      %if %bquote(&competing_risk)^= %then %do;
                              xbeta_cr1 = xbeta_cr1 + est_cr[2] * var_cr[2] ;                              
                      %end;
                  %end;
                  %if &amsm_type = 1 %then %do;
                       %if &numAMSM = 1 %then %do;
                            if &time = &time_start then &AMSM = 0 ;
                            else if &time > &time_start then &AMSM = 1 ;
                            xbeta_sw1 = xbeta_sw1 + est_sw[2]* &AMSM ;
                            %if %bquote(&competing_risk)^= %then %do;
                              xbeta_cr1 = xbeta_cr1 + est_cr[2] * var_cr[2] ;                              
                            %end;
                       %end;
                       %else %if &numAMSM = 2 %then %do;
                            %scan(&AMSM,1,%str( )) = 1 ;
                            if &time = &time_start then %scan(&AMSM,2,%str( )) = 0;
                            else %scan(&AMSM,2,%str( )) = 1 ;
                            xbeta_sw1 = xbeta_sw1 + est_sw[2]* var_sw[2] + est_sw[3]*var_sw[3] ;
                            %if %bquote(&competing_risk)^= %then %do;
                              xbeta_cr1 = xbeta_cr1 + est_cr[2]*var_cr[2] + est_cr[3]*var_cr[3] ;                              
                            %end;

                       %end ;

                  %end;
                  %if &amsm_type = 2 %then %do;
                      /* amsm is a categorical variable depending length of time treated */
                       %do ii = 1 %to &numAMSM_cats ;
                          %scan(&AMSM_class,&ii,%str( )) = 0 ;
                           

                       %end;
                       %do ii = 1 %to %eval(&numAMSM_cats - 1) ;
                           if %scan(&amsm_knots,&ii,%str( )) < &time < %scan(&amsm_knots,%eval(&ii+1),%str( )) 
                                  then %scan(&amsm_class,&ii,%str( )) = 1 ;
                           xbeta_sw1 = xbeta_sw1 + est_sw[%eval(&ii+1)] * var_sw[%eval(&ii+1)] ;
                           %if %bquote(&competing_risk)^= %then xbeta_cr1 = xbeta_cr1 + 
                                 est_cr[%eval(&ii+1)]*var_cr[%eval(&ii+1)] ;;
                       %end;
                       if %scan(&amsm_knots,&numAMSM_cats,%str( )) <= &time 
                                then %scan(&amsm_class,&numAMSM_cats,%str( )) = 1;
                       xbeta_sw1 = xbeta_sw1 + est_sw[%eval(&numAMSM_cats+1)] * var_sw[%eval(&numAMSM_cats + 1)];
                       %if %bquote(&competing_risk)^= %then xbeta_cr1 = xbeta_cr1 + 
                            est_cr[%eval(&numAMSM_cats+1)] * var_cr[%eval(&numAMSM_cats + 1)];;

                       
                  %end;
                  %else %if &amsm_type = 3 %then %do;
                       %Amsm_user_def_macro1 ; 
                       do j = 2 to %eval(&numAMSM + 1);
                            xbeta_sw1 = xbeta_sw1 + est_sw[j] * var_sw[j] ;
                       end;
                          /* 7-2014 added following since it was missing in this case*/
                        %if %bquote(&competing_risk)^= %then %do;
                            do j = 2 to %eval(1+&numAMSM ) ;
                               xbeta_cr1 = xbeta_cr1 + est_cr[j] * var_cr[j]; /* 7-2014 fixed loop for xbeta_sw to xbeta_sw0 */
                            end;
                        %end;
                  %end;


                  
                 * interactions with amsm ;
                 %if &ndef_amsm > 0 %then %do;
                      %do ninter = 1 %to &ndef_amsm ;
                         %let word =  %scan(&amsm_inter_def,%eval(&ninter),%str(:)) ;
                         &word ;
                      %end;
                      %let start0 = %eval(1 + &numAMSM + 1);

                       do j = &start0 to %eval(&_startpos - 1) ;
                            xbeta_sw1 = xbeta_sw1 + est_sw[j] * var_sw[j];
                       end;

                       %if %bquote(&competing_risk)^= %then %do;
                            do j = &start0 to %eval(&_startpos - 1);
                                xbeta_cr1 = xbeta_cr1 + est_cr[j] * var_cr[j];
                            end;
                       %end;
                %end;
                                           
                  /* calculate the probabilities */

                  pr_sw0 = 1.0/(1.0+exp(-xbeta_sw0));
                  pr_sw1 = 1.0/(1.0+exp(-xbeta_sw1));
                  s_sw0 = 1.0 - pr_sw0 ;
                  s_sw1 = 1.0 - pr_sw1 ;


      
                  %if %bquote(&competing_risk)^= %then %do;
                     pr_cr0 = 1.0/(1.0+exp(-xbeta_cr0));
                     pr_cr1 = 1.0/(1.0+exp(-xbeta_cr1));
                     s_cr0 = 1.0 - pr_cr0;
                     s_cr1 = 1.0 - pr_cr1 ;                     
                  %end;
                  %else %do;
                     pr_cr0 = 0 ;
                     pr_cr1 = 0 ;
                     s_cr0 = 1 ;
                     s_cr1 = 1 ;
                  %end;


                  /* for cumulative incidence use KM_12 at previous time point */
                  /* 7-2014 rwl changed to match gformula changes. add in factors for competing risk s_cr 
                     cumsurv * p(outcome) * (1-p(competing risk) */
                  inc_treated         = surv_treated * pr_sw1 * s_cr1 ; 
                  inc_nontreated      = surv_nontreated * pr_sw0 * s_cr0;

                  cuminc_treated      = cuminc_treated + inc_treated ;
                  cuminc_nontreated   = cuminc_nontreated + inc_nontreated ;

                  /* cumulative surv */
                  surv_nontreated     = surv_nontreated * s_sw0 * s_cr0 ;
                  surv_treated        = surv_treated * s_sw1 * s_cr1 ;

                  /* same values at each time */
                 %if &time_start = 0 %then %do;
                     _index_ = &time + 1 ;
                 %end;
                 %else %do;
                     _index_ = &time - &time_start + 1;
                 %end;   
                  ci_treated[_index_] = cuminc_treated ;
                  ci_nontreated[_index_] = cuminc_nontreated ;

                  km_treated[_index_] = surv_treated ;
                  km_nontreated[_index_] = surv_nontreated ;



                  %if &use_natural_course=1 %then %do;
                      pr_nc = 1.0/(1.0+exp(-xbeta_nc));
                      s_nc = 1.0 - pr_nc;
                      pr_nc_cr = 0 ;
                      s_nc_cr = 1 ;
                      %if %bquote(&competing_risk)^= %then pr_nc_cr = 1.0/(1.0+exp(-xbeta_nc_cr)) ;;
                      %if %bquote(&competing_risk)^= %then s_nc_cr = 1.0 - pr_nc_cr ;;
                      inc_nc = surv_nc * pr_nc * s_nc_cr ; /* 7-2014 rwl add in factor for competing risk to match gformula */
                      cuminc_nc = cuminc_nc + inc_nc ;
                      surv_nc = surv_nc * s_nc * s_nc_cr ;

                      km_nc[_index_] = surv_nc ;
                      ci_nc[_index_] = cuminc_nc ;


                  %end;

                        
           end;
           run;
 
           proc means data = curves noprint  ;
           var km_treated&time_start-km_treated&time_end;
           freq numberhits ;         
           output out=_km_treat(drop = _TYPE_ _FREQ_ ) mean= ;
           run;

           proc transpose data = _km_treat out=_km_1 ;
           run;

           data _km_1 ;
           set _km_1 ;
           &time = input(substr(_NAME_,11),8.);
           drop _NAME_ ;
           run;

           proc means data = curves noprint ;
           var  km_nontreated&time_start-km_nontreated&time_end ;
           freq numberhits ;
           output out = _km_nontreat(drop = _TYPE_ _FREQ_) mean= ;
           run;

           proc transpose data = _km_nontreat out = _km_0 ;
           run;

           data _km_0 ;
           set _km_0 ;
           &time = input(substr(_NAME_,14),8.);
           drop _NAME_ ;
           run;

           proc means data = curves noprint ;
           var ci_treated&time_start-ci_treated&time_end ;
           freq numberhits ;
           output out =_ci_treat(drop = _TYPE_ _FREQ_) mean= ;
           run;

           proc transpose data = _ci_treat out=_ci_1 ;
           run;

           data _ci_1 ;
           set _ci_1 ;
           &time = input(substr(_NAME_,11),8.);
           drop _NAME_ ;
           run;

           proc means data = curves noprint ;
           var ci_nontreated&time_start-ci_nontreated&time_end  ;
           freq numberhits ;
           output out=_ci_nontreat(drop = _TYPE_ _FREQ_ ) mean= ;
           run;

           proc transpose data = _ci_nontreat out=_ci_0 ;
           run;

           data _ci_0 ;
           set _ci_0 ;
           &time = input(substr(_NAME_,14),8.);
           drop _NAME_ ;
           run;

           %if &use_natural_course = 1 %then %do;
                 proc means data = curves noprint ;
                 var ci_nc&time_start-ci_nc&time_end  ;
                 freq numberhits ;
                 output out=_ci_nc(drop = _TYPE_ _FREQ_ ) mean= ;
                 run;

                 proc transpose data = _ci_nc out=_ci_nc_t ;
                 run;

                 data _ci_nc_t ;
                 set _ci_nc_t ;
                 &time = input(substr(_NAME_,6),8.);
                 drop _NAME_ ;
                 run;

                 proc means data = curves noprint ;
                 var km_nc&time_start-km_nc&time_end  ;
                 freq numberhits ;
                 output out=_km_nc(drop = _TYPE_ _FREQ_ ) mean= ;
                 run;

                 proc transpose data = _km_nc out=_km_nc_t ;
                 run;

                 data _km_nc_t ;
                 set _km_nc_t ;
                 &time = input(substr(_NAME_,6),8.);
                 drop _NAME_ ;
                 run;

           %end;

           data _onegraph_;
           merge _km_0 ( rename = (col1=km_nontreat ) )
                 _km_1 ( rename = (col1 = km_treat  ))
                 _ci_0 ( rename = (col1 = ci_nontreat ))
                 _ci_1 ( rename = (col1 = ci_treat  ))
                 %if &use_natural_course = 1 %then %do;
                    _ci_nc_t (rename = (col1 = ci_nc ))
                    _km_nc_t (rename = (col1 = km_nc))
                 %end;
 
              ;

           by &time ;

           /* carry forward values to fill in possible gaps when stepsize is not 1 */
           %if &time_step > 1 %then %do;
               lkm0 = lag(km_nontreat);
               lkm1 = lag(km_treat);
               lci0 = lag(ci_nontreat);
               lci1 = lag(ci_treat);
           
               if km_nontreat = . then km_nontreat = lkm0 ;
               if km_treat    = . then km_treat = lkm1 ;
               if ci_nontreat = . then ci_nontreat = lci0 ;
               if ci_treat    = . then ci_treat = lci1 ;

               %if &use_natural_course = 1 %then %do;
                    lkmnc = lag(km_nc);
                    lcinc = lag(ci_nc);
                    if km_nc = . then km_nc = lkmnc ;
                    if ci_nc = . then ci_nc = lcinc;
                    drop lkmnc lcinc ;
               %end;
               drop lkm0 lkm1 lci0 lci1 ;
           %end;
           risk_diff = &riskdiff1 - &riskdiff0 /*ci_nontreat - ci_treat */ ;
           bsample = &bsample ;
           run;

           %if &bsample = 0 %then %do;
               
                proc print data = _onegraph_ ;
                %if %bquote(&competing_risk)= %then title "Survival probabilities for &outcMSM . ";;
                %if %bquote(&competing_risk)^= %then title "Survival probabilities for &outcMSM with competing risk ( &competing_risk ). ";;
                var &time km_treat km_nontreat  
                    %if &use_natural_course = 1 %then km_nc  ;
                    %if &print_ci = 1 %then ci_treat ci_nontreat 
                          %if &use_natural_course = 1 %then  ci_nc ; 
                           
                    risk_diff  ;
                run;
           %end;

               
           %if &debug = 0 %then %do;
                  proc datasets library = work nolist ;
                  delete curves _ci_0 _ci_1 _km_0 _km_1  _ci_treat _ci_nontreat _km_treat _km_nontreat
                      %if &use_natural_course = 1 %then %do;
                            pe_nc emppest_nc  _km_nc _ci_nc  _ci_nc_t _km_nc_t
                            %if %bquote(&competing_risk)^= %then  emppest_nc_cr pe_nc_cr ;
                      %end;
                      %if %bquote(&competing_risk)^= %then emppest_cr pe_cr ;
                      ;
                  quit;


           %end;


%mend;
%macro extra_variables_def ;
   /* any user defined variables used in constructing variables used in weighted model. */
   /* by default this macro is empty */

%mend ;

%macro bootstrap_results(
       bootlib = ,
       bootname= ,
       modeltype = sw ,  
       numparts = 1,
       samplestart = 0,
       numboot = 200,
       time = ,
       use_natural_course = 1,
       outcMSM = ,
       competing_risk = ,
       riskdiff1 = ci_treat ,
       riskdiff0 = ci_nc ,
       print_ci = 0 
       );

     %local datalist i j part bstart0 bstart1 bend1 betas ;

     %if &numparts = 1 %then %do;
         %let datalist = &bootlib..&bootname._&modeltype._0_&numboot ;
     %end;
     %else %do;
         %let bstart0 = %scan(&samplestart,1,%str( ));
         %let bstart1 = %scan(&samplestart,2,%str( )); 
         %let bend1 = %eval(&bstart1 - 1);
         %let datalist = &bootlib..&bootname._&modeltype._&bstart0._&bend1 ;
         %do part = 2 %to &numparts ;
              %let bstart0 = %scan(&samplestart,%eval(&part),%str( ));
              %if &part < &numparts %then %do;
                   %let bstart1 = %scan(&samplestart,%eval(&part+1),%str( ));
                   %let  bend1 = %eval(&bstart - 1);
              %end;
              %else %let bend1 = &numboot ;
              %let datalist = &datalist &bootlib..&bootname._&modeltype._&bstart0._&bend1 ;
         %end;
     %end;

     data _tmp ;
     set &datalist ;
     by bsample ;
     run;

     data _sample0;
     set _tmp (where = (bsample = 0));
     run;

     %if %upcase(&modeltype) ^= SURV %then %do;

          proc transpose data = _sample0(drop = bsample) out = _sample0( rename = (_NAME_ = name _LABEL_=name0 col1=estimate));
          run;

          proc sql noprint ;
            select name into :betas separated by ' ' from _sample0;
          quit;

          data _sample0 ;
          set _sample0 ;
          num = _n_ ;
          name = upcase(name);
          run;

          proc sort data = _sample0 ;
          by name ;
          run;

          proc corr data = _tmp( where = ( bsample > 0))  cov nocorr nosimple outp = _cov_&modeltype noprint ;
          var  &betas ;                                                                            
          run; 

 

          data _myout2;
          set _cov_&modeltype (drop = _name_ );
          if _type_ = 'MEAN' or _type_ = 'STD' ;                                          
          run;

          proc transpose data = _myout2 out = _myout2t (rename = ( _name_ = name col1 = mean col2 = std));                     
          run;

          data _myout2t;
          set _myout2t;
          name = upcase(name);
          run;

          proc sort data = _myout2t sortsize= 256M; 
          by name ; 
          run;
 

         data _sample0 ;
         merge _sample0 _myout2t;
         by name ;
         run;

         proc sort data = _sample0 ;
         by num ;
         run;


         data _msm&modeltype;
         set _sample0;       
         label p_value = 'Pr > |Z|'
                name0="Variable"
                estimate="Estimate"
                bbias = "Bias"
                std="Standard Error"
                z="Z-score"
                lb="95% lower bound"
                ub="95% upper bound" 
           ;
         format p_value PVALUE8.4 ;                  
         lb = estimate - 1.96 * std;
         ub = estimate + 1.96 * std; 
         if std = 0 | std = . then do ;
             z = . ;
             p_value = . ;
         end;
         else do ;
             z = estimate / std;
             p_value = 2*(1-probnorm(abs(z))); 
         end;
         bbias = estimate - mean; ;
         if _n_ = 1 then name0 = "Intercept";
         drop num ; 
         run;   

       
         proc print data = _msm&modeltype noobs label;
         %if %upcase(&modeltype)=SW %then   
                title1 "Results for marginal structural model with stabilized weights for &outcMSM";;
         %if %upcase(&modeltype)=CR %then   
                title1 "Results for marginal structural model with stabilized weights for &competing_risk";;
         %if %upcase(&modeltype)=NC %then   
                title1 "Results for marginal structural model with stabilized weights for &outcMSM under natural course";;
  
        title2 "with &numboot bootstrap samples";
        var name0  estimate bbias std  lb ub p_value ;
        run;


        proc datasets library = work nolist ;
        delete _tmp _sample0 _myout2 _myout2t ;
        quit ;

    %end;
    %else %do;

       /* for graphs */

       data _bgraphs ;
       set _tmp ;
       run;

       proc means data = _bgraphs (where = (bsample > 0))mean std clm noprint ;
       var km_nontreat km_treat ci_nontreat ci_treat
          %if &use_natural_course = 1 %then  ci_nc km_nc ;
           risk_diff;
       class &time ;
       types &time  ;
       output out = _bgraphsm mean= std= /autoname ;
       run;
   

       data _bgraphsr ;
       merge _bgraphs (where = (bsample = 0) ) _bgraphsm ;
       by &time ;

       label km_treat = "Survival Probability"


             km_nontreat = "Survival probability"
             km_nontreat_bias = "Bias"
             km_nontreat_StdDev = "Standard error"
             km_nontreat_lb = "95% CI lower bound"
             km_nontreat_ub = "95% CI upper bound"

             km_treat = "Survival probability"
             km_treat_bias = "Bias"
             km_treat_StdDev = "Standard error"
             km_treat_lb = "95% CI lower bound"
             km_treat_ub = "95% CI upper bound"

             ci_nontreat = "Cummulative incidence"
             ci_nontreat_bias = "Bias"
             ci_nontreat_StdDev = "Standard Error"
             ci_nontreat_lb = "95% CI lower bound"
             ci_nontreat_ub = "95% CI upper bound"

             ci_treat = "Cummulative incidence"
             ci_treat_bias = "Bias"
             ci_treat_StdDev = "Standard error"
             ci_treat_lb = "95% CI lower bound"
             ci_treat_ub = "95% CI upper bound"

             risk_diff = "Risk difference"
             risk_diff_bias = "Bias"
             risk_diff_StdDev = "Standard error"
             risk_diff_lb = "95% CI lower bound"
             risk_diff_ub = "95% CI upper bound"

             %if &use_natural_course = 1 %then %do;
                 km_nc = "Survival probability"
                 km_nc_bias = "Bias"
                 km_nc_StdDev = "Standard error"
                 km_nc_lb = "95% CI lower bound"
                 km_nc_ub = "95% CI upper bound"

                 ci_nc = "Cummulative incidence"
                 ci_nc_bias = "Bias"                     
                 ci_nc_StdDev = "Standard error" 
                 ci_nc_lb = "95% CI lower bound"              
                 ci_nc_ub = "95% CI upper bound" 

             %end;
             ;

             
       km_treat_lb = km_treat - 1.96*km_treat_StdDev ;
       km_treat_ub = km_treat + 1.96*km_treat_StdDev ;
       km_treat_bias = km_treat - km_treat_mean ;

       km_nontreat_lb = km_nontreat - 1.96*km_nontreat_StdDev ;
       km_nontreat_ub = km_nontreat + 1.96*km_nontreat_StdDev ;
       km_nontreat_bias = km_nontreat - km_nontreat_mean ;


       ci_treat_lb = ci_treat - 1.96*ci_treat_StdDev ;
       ci_treat_ub = ci_treat + 1.96*ci_treat_StdDev ;
       ci_treat_bias = ci_treat - ci_treat_mean ;

       ci_nontreat_lb = ci_nontreat - 1.96*ci_nontreat_StdDev ;
       ci_nontreat_ub = ci_nontreat + 1.96*ci_nontreat_StdDev ;
       ci_nontreat_bias = ci_nontreat - ci_nontreat_mean ;

       %if &use_natural_course = 1 %then %do;
           ci_nc_lb = ci_nc - 1.96*ci_nc_StdDev ;
           ci_nc_ub = ci_nc + 1.96*ci_nc_StdDev ;
           ci_nc_bias = ci_nc - ci_nc_mean ;

           km_nc_lb = km_nc - 1.96*km_nc_StdDev ;
           km_nc_ub = km_nc + 1.96*km_nc_StdDev ;
           km_nc_bias = km_nc - km_nc_mean ;
       %end;

       risk_diff_lb = risk_diff - 1.96* risk_diff_StdDev ;
       risk_diff_ub = risk_diff + 1.96 * risk_diff_StdDev ;
       risk_diff_bias = risk_diff - risk_diff_mean ;

       run;

       %if &print_ci = 1 %then %do;     
           proc print data = _bgraphsr noobs label ;
           title "Cummulatice incidence for &outcMSM under  never treated using &numboot samples. ";
           var &time  ci_nontreat ci_nontreat_bias ci_nontreat_StdDev ci_nontreat_lb ci_nontreat_ub ;
           run;

       
           proc print data = _bgraphsr noobs label ;
           title "Cummulatice incidence for &outcMSM under always treated using &numboot samples. ";
           var &time ci_treat ci_treat_bias ci_treat_StdDev ci_treat_lb ci_treat_ub ;
           run;

           %if &use_natural_course = 1  %then %do;
        
               proc print data = _bgraphsr noobs label ;
               title "Cummulatice incidence  for &outcMSM under natural course using &numboot samples. ";
               var &time  ci_nc ci_nc_bias ci_nc_StdDev ci_nc_lb ci_nc_ub ;
               run;
           %end;

       %end;
     
       proc print data = _bgraphsr noobs label ;
       title "Survival probabilities  for &outcMSM under  never treated using &numboot samples. ";
       var &time  km_nontreat km_nontreat_bias km_nontreat_StdDev km_nontreat_lb km_nontreat_ub ;
       run;

      
       proc print data = _bgraphsr noobs label ;
       title "Survival probabilities for &outcMSM under always treated using &numboot samples. ";
       var &time km_treat km_treat_bias km_treat_StdDev km_treat_lb km_treat_ub ;
       run;

       %if &use_natural_course = 1 %then %do;
          
           proc print data = _bgraphsr noobs label ;
           title "Survival probabilities for &outcMSM under then natural course using &numboot samples. ";
           var &time  km_nc km_nc_bias km_nc_StdDev km_nc_lb km_nc_ub ;
           run;
       %end;
      
       proc print data = _bgraphsr noobs label ;
       title "Risk difference of &riskdiff1 - &riskdiff0 using &numboot samples.";
       var &time  risk_diff risk_diff_bias risk_diff_StdDev risk_diff_lb risk_diff_ub ;
       run;

       proc datasets library = work nolist ;
       delete _tmp _sample0 _bgraphs _bgraphsm  ;
       quit;
    %end;




%mend ;

%macro checksw(datain = );
                  %local dsid chk rc ;
                   %let dsid = %sysfunc(open(&datain));
                   %let chk = %sysfunc(varnum(&dsid,_ESTTYPE_));
                   %let rc = %sysfunc(close(&dsid));
                  
                   %if &chk > 0 %then %do;
                          data &datain ;
                          set  &datain (drop = _ESTTYPE_) ;
                          run;
                    %end;  
              

%mend ;


%macro numargs(arg);

     %let n = 1;
     %if %bquote(&arg)^= %then %do;
          %do %until (%qscan(&arg,%eval(&n),%str( ))=%str());
           
               %let word = %qscan(&arg,&n);
               %let n = %eval(&n+1);

          %end;

     %end;
     %eval(&n-1) /* there is no ; here since it will be used as %let a = %numargs(&b) ;
                      and the ; is included at the end of this line  */                               
                                                
%mend numargs;



%macro remove_class;

     /* remove entries in list2 from list1 */
     %local i j  test word ;
     %let i = 1;
     %let _nonclass= ;
     %do %until(%scan(&_allcov,&i,%str( )) = %str() );

       %let test = %scan(&_allcov,&i);
       %let j = 1;
       %let add = 1;
       %do %until(%scan(&_class,&j,%str( ))=%str());
          %let word = %scan(&_class,&j);
           %if &test = &word %then %let add = 0 ;
           %let j = %eval(&j + 1);
       %end;
       %if %eval(&add) = 1 %then %let _nonclass = &_nonclass &test ;;
       %let i = %eval(&i + 1);

     %end;

%mend remove_class;

/*****/
%macro create_class_list(cov_list = , parent_list = , class_list = );
     %local tmp_list  n_parent  ncov ;
     %let tmp_list = ;
     %let n_parent = %numargs(&parent_list);
     %let ncov = %numargs(&cov_list);
     %do i = 1 %to %eval(&ncov);
          %let vartmp0 = %scan(&cov_list,%eval(&i),%str(' '));
          %do j = 1 %to %eval(&n_parent);
               %let vartmp1 = %scan(&parent_list,%eval(&j),%str(' '));
               %if %upcase(&vartmp1) = %upcase(&vartmp0) %then %do;
                    %let j = %eval(&n_parent + 1) ;
                    %let tmp_list = &tmp_list &vartmp0 ;
               %end;
          %end;
     %end;
     %let &class_list = &tmp_list ;
%mend;
/****/
%macro remove_extraneous;
  
     %do i = 1 %to &nclass;
          /* class var is a list of macro variables */
          %let _var  = %qscan(&class_var, &i,%str( ));
          %let _cvar = %qscan(&cov_var,&i,%str( ));
          %let _class  = %unquote(&_var);
          %let _allcov = %unquote(&_cvar);
          %if %bquote(&_class)^= %then %do; 

               /* remove entries in list1 that are not in list2 */
            %local   test word ;
            %let ii = 1;
            %let _class_used= ;
            %do %until(%scan(&_class,&ii,%str( )) = %str() );

                %let test = %scan(&_class,&ii);
                %let j = 1;
                %let remove = 1;
                %do %until(%scan(&_allcov,&j,%str( ))=%str());
                      %let word = %scan(&_allcov,&j);
                      %if &test = &word %then %let remove = 0 ;
                      %let j = %eval(&j + 1);
                %end;
                %if %eval(&remove) = 0 %then %let _class_used = &_class_used &test ;;
                %let ii = %eval(&ii + 1);

            %end;

           %let %substr(&_var,2) = &_class_used;

        %end;  
             
     %end;

%mend remove_extraneous;

/*************************************************************************/


 
%macro convert_class ;

     %local cnt uniquecnt newlist _nonclass i in test testword add tmp keyfld ;
     %let keyfld =  %unquote(&class_var) ;

     /* create list of unique elements of all the class variables */
     %let i = 1;

     %do %until(%scan(&keyfld,&i,%str( ) ) = %str() );
          %let var&i = %scan(&keyfld,&i,%str( ));
          %let i = %eval(&i + 1);
     %end;

     %let cnt = %eval(&i - 1);
     %let newlist = &var1 ;
     %let in = 1;
     %let test = 2;
 
     %do %while(%eval( &test) <= %eval(&cnt));

          %let testword = &&var&test;
          %let add = 1;

          %do i = 1 %to &in;
               %let tmp = %scan(&newlist,&i,%str( ));
               %if &testword = &tmp  %then %do;
                    %let add = 0;
                    %let i = %eval(&in + 1);
               %end;
          %end;

          %if &add = 1 %then %do;
               %let newlist = &newlist &testword ;
               %let in = %eval(&in + 1);
          %end;
    
          %let test = %eval(&test + 1);

     %end;
 
     %let uniquecnt = &in;


     %let newlist2= ;

   proc glmmod data=_cero_(keep = &id &time &newlist)  outparm=_Parm_ outdesign=_Design_ ;
   ods exclude ClassLevels NObs Parameters DesignPoints;
   ods output ClassLevels=_cl_
      ;
      class &newlist;
      model &id &time =  &newlist /noint;
   run;



   proc contents data = _design_ varnum out=_test_ (keep=name varnum label) noprint;
   run;

   data _test_ ;
   set _test_ ;
   label=compress(label);
   run;

   %local label cat_levels cat_values  total num_cat _newlist word ;

   proc sql noprint;
   select label into :label 
     separated by ' ' from _test_;
   select name into :names
     separated by ' ' from _test_;
   select levels into :cat_levels 
     separated by ' : ' from _cl_ ;
   select values into :cat_values 
     separated by ' : ' from _cl_ ;
   quit;

   proc means data = _cl_ sum n  noprint  ;
   var levels ;
   output out = _num_bin_ sum=total n = num_cat;
   run;
 
   data _null_;
   set _num_bin_;
   if _n_ = 1 then do;
         call symput("total",trim(left(total)));
         call symput("num_cat",trim(left(num_cat)));
   end;      
   run;

 

  proc datasets library=work nolist ;
  modify _design_;
  rename %do i = 1 %to &total ;
      %let word = %scan(&label,&i,%str( ));
      col&i = &word
      %end;
      ;
  quit;

  %local word1 val lev j first last word2 _newlist ;
  %let _newlist = ;
/*
  %put old variables =  &newlist ;
  %put old categorical values = &cat_values ;
  %put number of levels per variable =  &cat_levels ;
*/

  /*****
  %do i = 1 %to &num_cat ;
     %let word1 = %scan(&newlist,&i,%str( ));
     %let val = %scan(&cat_values,&i,%str(:));
     %let lev = %scan(&cat_levels,&i,%str(:)); 
     %let j = %eval(&lev) ;

     %put lev = &lev ;
     %put val = &val ; 
     %if %upcase(&class_ref)=FIRST %then %do;
        %let first=2;
        %let last = %eval(&j);
     %end;
     %else %if %upcase(&class_ref)=LAST %then %do;
         %let first = 1 ;
         %let last = %eval(&j - 1) ;
     %end;
     %do k = &first %to &last ;
          %let word2 = %scan(&val,&k,%str( ));
          %let _newlist = &_newlist &word1.&word2 ;
         
     %end;
    
 
 
  %end;
***/
   %let num_old = %numargs(&newlist) ;
   %do i = 1 %to &nclass;


          %let var = %qscan(&class_var, &i,%str( ));
          %let cvar = %qscan(&cov_var,&i,%str( ));
          %let _class = %unquote(&var);
          %let _allcov = %unquote(&cvar);
          %if %bquote(&_class)^= %then %do;
             %remove_class;
             %let j = 1;
             %let in = 1;
             %let tmp = ;
             %let ctmp= ;
             %let _new_class_tmp = ;
            %do %until(%scan(%unquote(&var),&j,%str( ) ) = %str() );

                 %let vartmp = %scan(%unquote(&var),&j,%str( ));
                 %let j = %eval(&j + 1);
                 

                 /* find position in list of old categorical variables */

                 %let found = 0 ;
                 %let position = 1 ;

                 %do %while (&found = 0 and &position <= &num_old ) ;
                      %let word = %scan(&newlist,&position,%str( ));
                      %if %upcase(&vartmp) = %upcase(&word) %then %do;
                           %let found = &position ;
                      %end;
                      %else %do;
                          %let position = %eval(&position + 1);
                      %end;
                 /*     %put word = &word vartmp = &vartmp found = &found position = %eval(&position-1) ;*/
                 %end;
                 %if &found > 0 %then %do;
                      %let levels_tmp = %scan(&cat_levels,&position,%str(:));
                      %let values_tmp = %scan(&cat_values,&position,%str(:));

                      %if %upcase(&class_ref) = FIRST %then %let k0 = 2 ;
                      %else %let k0 = 1 ;
                      %if %upcase(&class_ref) = LAST %then %let k1 = %eval(&levels_tmp - 1);
                      %else %let k1 = &levels_tmp ;
                 
                           
             
              /*   %let new_levels = %eval(&&_&vartmp.levels - 1);*/
                     %do k = &k0 %to &k1;
                         %let value = %scan(&values_tmp,&k,%str( ));
                         %let tmp = &tmp &&vartmp.&value ;
                         %let _new_class_tmp = &_new_class_tmp &&vartmp.&value ; 
                     %end;
                     
                 %end;
            %end;
  
            %let %substr(&cvar,2) = &_nonclass &tmp;
            %let %substr(&var,2) = &_new_class_tmp ;

        %end;
/**
            %put %substr(&cvar,2) = &_nonclass &tmp;
            %put %substr(&var,2) = &_new_class_tmp ;
**/
     %end;
  
    data _cero_ ;
    merge _cero_ _design_;
    by &id &time ;
    run;

    proc datasets library = work nolist ;
    delete _num_bin_ _parm_ _test_ _cl_ _design_ ;
    quit; 
   
%mend;
/*****/
 /*MACRO RCSPLINE

   For a given variable named X and from 3-10 knot locations,
   generates SAS assignment statements to compute k-2 components
   of cubic spline function restricted to be linear before the
   first knot and after the last knot, where k is the number of
   knots given.  These component variables are named c1, c2, ...
   ck-2, where c is the first 7 letters of X.

   Usage:

   DATA; ....
   %RCSPLINE(x,knot1,knot2,...,norm=)   e.g. %RCSPLINE(x,-1.4,0,2,8)

        norm=0 : no normalization of constructed variables
        norm=1 : divide by cube of difference in last 2 knots
                 makes all variables unitless
        norm=2 : (default) divide by square of difference in outer knots
                 makes all variables in original units of x

   Reference:

   Devlin TF, Weeks BJ (1986): Spline functions for logistic regression
   modeling. Proc Eleventh Annual SAS Users Group International.
   Cary NC: SAS Institute, Inc., pp. 646-51.


   Author  : Frank E. Harrell Jr.
             Clinical Biostatistics, Duke University Medical Center
   Date    : 10 Apr 88
   Mod     : 22 Feb 91 - normalized as in S function rcspline.eval
             06 May 91 - added norm, with default= 22 Feb 91
             10 May 91 - fixed bug re precedence of <>

                                                                      */
%MACRO RCSPLINE(x,knot1,knot2,knot3,knot4,knot5,knot6,knot7,
                  knot8,knot9,knot10, norm=2);
%LOCAL j v7 k tk tk1 t k1 k2;
%LET v7=&x; %IF %LENGTH(&v7)=8 %THEN %LET v7=%SUBSTR(&v7,1,7);
  %*Get no. knots, last knot, next to last knot;
    %DO k=1 %TO 10;
    %IF %QUOTE(&&knot&k)=  %THEN %GOTO nomorek;
    %END;
%LET k=11;
%nomorek: %LET k=%EVAL(&k-1); %LET k1=%EVAL(&k-1); %LET k2=%EVAL(&k-2);
%IF &k<3 %THEN %PUT ERROR: <3 KNOTS GIVEN.  NO SPLINE VARIABLES CREATED.;
%ELSE %DO;
 %LET tk=&&knot&k;
 %LET tk1=&&knot&k1;
 DROP _kd_; _kd_=
 %IF &norm=0 %THEN 1;
 %ELSE %IF &norm=1 %THEN &tk - &tk1;
 %ELSE (&tk - &knot1)**.666666666666; ;
    %DO j=1 %TO &k2;
    %LET t=&&knot&j;
    &v7&j=max((&x-&t)/_kd_,0)**3+((&tk1-&t)*max((&x-&tk)/_kd_,0)**3
        -(&tk-&t)*max((&x-&tk1)/_kd_,0)**3)/(&tk-&tk1)%STR(;);
    %END;
 %END;
%MEND;



%macro plot_graphs (datain = _bgraphsr,
                    time = time ,    
                    outcMSM= ,
                    natural_course = 0 , 
                    sixgraphs=0 ,
                    gfilename=fourgraphs.pdf ,
                    title1=,title2=,title3=,titledata=,tsize=0.75 ,
                    frombootstrap = 1) ;
 



 proc greplay igout = GSEG nofs ;
 delete _all_ ;
 run;
 quit;

%local nperiods mintime mydate mywork  outcome ;
%let outcome = &outcMSM ;



%let mydate = &sysdate9._%scan(&systime,1,':')_%scan(&systime,2,':') ;
%let mywork = %sysfunc(getoption(work)) ;       
     filename mygraphs "&mywork./graphs.pdf";
       
 
 

    %let maxnc = -1000000 ; /* need a default value when natural_course = 0 */
    proc sql noprint ;
    select count(&time)      as mtime into :nperiods from &datain ;
    select min(&time)        as mintime into :mintime from &datain ;   
    select max(&time)     as maxtime into :maxtime from &datain ;
    %if &frombootstrap = 1 %then %do;
        select min(risk_diff_lb) as mindiff into :mindiff from &datain ;
        select max(risk_diff_ub) as maxdiff into :maxdiff from &datain ; 
    %end;
    %else %do;
         select min(risk_diff) as mindiff into :mindiff from &datain ;
         select max(risk_diff) as maxdiff into :maxdiff from &datain ;
    %end;
    select min(ci_treat) as minctreat into :minctreat from &datain ;
    select min(ci_nontreat) as minnontreat into :minnontreat from &datain ;
    %if &natural_course = 1 %then %do;
         select min(ci_nc) as minnc into :minnc from &datain ;
         select max(ci_nc) as maxnc into :maxnc from &datain ;
    %end ;
    select max(ci_treat) as maxctreat into :maxctreat from &datain ;
    select max(ci_nontreat) as maxnontreat into :maxnontreat from &datain ;

    quit;
    
    %let mindiff = %sysevalf(%sysfunc(round(&mindiff,0.01))-0.005);
    %let maxdiff = %sysevalf(%sysfunc(round(&maxdiff,0.01))+0.005);
    %let mydiff = %sysfunc(round((&maxdiff - &mindiff)/5,0.001));
     
        
    %let nperiods = %sysfunc(compress(&nperiods)) ;
    %let mintime = %sysfunc(compress(&mintime)) ;
    %let timediff = %sysfunc(round((&maxtime - &mintime)/6,1.0));

    data _null_ ;
    maxci =  max( &maxctreat , &maxnontreat , &maxnc )  ;
    maxci = round(maxci+0.01,0.01) ;
    minci = 0 ;
     
    put maxci= minci= ;
    mydiffci = round(maxci /6 ,0.01);
    call symput("maxci",trim(left(maxci)));
    call symput("minci",trim(left(minci)));
    call symput("mydiffci",trim(left(mydiffci)));
    run;


 goptions reset=all noborder device=pdf gsfname=mygraphs  ;
 proc greplay tc=work.tempcat nofs;
 tdef newtemp des="my two panel template"
     1/llx=20   lly=52
       ulx=20   uly=90
       urx=83   ury=90
       lrx=83   lry=52
  
           

     2/llx=20   lly=10
       ulx=20   uly=48
       urx=83   ury=48
       lrx=83   lry=10
 

        ;
        tdef newtemp2 des="my five panel template"
     1/llx=0   lly=52
       ulx=0   uly=90
       urx=83  ury=90
       lrx=83  lry=52
        

     2/llx=0   lly=10
       ulx=0   uly=48
       urx=48  ury=48
       lrx=48  lry=10
       

     3/llx=52 lly=52
       ulx=52 uly=90
       urx=100 ury=90
       lrx=100 lry=52
       

     4/llx=52 lly=10
       ulx=52 uly=48
       urx=100 ury=48
       lrx=100 lry=10
      
      5/llx=0   lly = 0  
        ulx=0   uly=100
        urx=100 ury=100
        lrx=100 lry=0
        ;

  template newtemp;
 list template;
 quit;


 %if %bquote(&title1)= %then %do;

  %let period = ;
  %if &frombootstrap = 0 %then %let period= . ;
        title1 height= &tsize  'Left column: observed (solid line), natural course (dotted line) estimates by follow-up period.';
        title2 height= &tsize   "Right column: differences between observed and natural course estimates (solid lines)&period  ";
        %if &frombootstrap = 1 %then title3 height=&tsize  'and 95% pointwise confidence intervals (dotted lines).';;
%end;
%else %do;
        title1 height= &tsize &title1 ;
        title2 height= &tsize &title2 ;
        title3 height= &tsize &title3 ;
%end; 

proc gslide gout=work.gseg ;
run;
quit;

title ;



/* extract data for graphs 1) cumulative incidence , 2) survival curves, 3) risk difference  
    any natural course will be included in the corresponding data sets  */

data cuminc ;
set &datain (keep= &time ci_treat  %if &frombootstrap = 1 %then ci_treat_lb ci_treat_ub ;
                       ci_nontreat  %if &frombootstrap = 1 %then ci_nontreat_lb ci_nontreat_ub ;
                       %if &natural_course = 1 %then ci_nc ; %if &frombootstrap = 1 and &natural_course = 1  %then ci_nc_lb ci_nc_ub  ;
                       );
label ci_treat="Treated"
      ci_nontreat = "Non-treated"
      %if &natural_course = 1 %then ci_nc="Natural Course";
      ;
run;

data surv ;
set &datain (keep = &time km_treat %if &frombootstrap = 1 %then km_treat_lb km_treat_ub ;
                          km_nontreat %if &frombootstrap = 1 %then km_nontreat_lb km_nontreat_ub ;
                            %if &natural_course = 1 %then km_nc ; %if &frombootstrap = 1 and &natural_course = 1  %then km_nc_lb km_nc_ub ;

);

run;

data riskdiff ;
set &datain (keep = &time risk_diff %if &frombootstrap = 1 %then risk_diff_lb risk_diff_ub  ; );
run;
                       
 
 

  
    goptions reset = all  /* display */ hsize=8 in vsize = 10 in device=pdf  gsfname=mygraphs  ;

    symbol1 line = 1 width = 3 interpol = line  color = red ;
    symbol2 line = 2 width = 3 interpol = line color = blue ;
    symbol3 line = 2 width = 3 interpol = line color = black ;

    axis1 order = &mintime to &maxtime by &timediff  minor = none label= ( h = 3 "&time" )value=(h=2) ;
    axis2 order = &minci to &maxci by &mydiffci minor = none value=(h=2) label = (h = 3  justify=center angle=90 "Cumulative &outcome incidence" ) ;
    axis3 minor = none value=(h=2) label= none  ;
    axis4 minor = none value=(h=2) label = (h = 3  justify=center angle=90 "Probability of &outcome  survival" ) ;
    axis5 minor = none value=(h=2) label = (h = 3  justify=center angle=90 "Risk difference" ) ;
     
    legend1 label= none
        
       
       %if &natural_course = 1 %then  value=(h=1.15 "Treated" "Non-treated" "Natural Course") ;
       %else value=(h=1.15 "Treated" "Non-treated" ) ;
        shape=symbol(4,2)
        position=(bottom center inside)
        mode=share
        down=1;

         proc gplot data = cuminc  ;
         plot   ci_treat * &time ci_nontreat * &time %if &natural_course = 1 %then ci_nc * &time ; / overlay   
              vaxis = axis2 haxis = axis1 name="outc" noframe legend=legend1;            
         run;
         quit;
           
         proc gplot data = surv  ;
         plot   km_treat * &time km_nontreat * &time %if &natural_course = 1 %then km_nc * &time ; / overlay   vaxis = axis4 haxis = axis1 name="outkm" noframe nolegend;            
         run;
         quit;
           
         symbol1 line = 1 width = 3 interpol = line  color = red ;
         symbol2 line = 2 width = 3 interpol = line color = black ;
         symbol3 line = 2 width = 3 interpol = line color = black ;
         legend2 label= none
        
       
            value=(h=1.15 "Risk difference" %if &frombootstrap = 1 %then "95% confidence bounds" ;)
            shape=symbol(4,2)
            position=(bottom left inside)
            mode=share
            down = 1

         ;

 axis5 minor = none value=(h=2) order = &mindiff to &maxdiff by &mydiff label = (h = 3  justify=center angle=90 "Risk difference" )  ;
 axis6 minor = none value=(h=2) order = &mindiff to &maxdiff by &mydiff label = (h = 3  justify=center angle=90 "Risk difference" )    ;

         proc gplot data = riskdiff  ;         
         plot  risk_diff * &time  %if &frombootstrap = 1 %then risk_diff_ub * &time  ; / overlay  
            haxis = axis1 vaxis= axis5 name = "rdiffu" noframe 
            legend = legend2;   
         
         run;
         quit;
           
    symbol1 line = 2 width = 3 interpol = line color = black ;

     %if &frombootstrap = 1 %then %do;
         proc gplot data = riskdiff  ;         
           
         plot risk_diff_lb * &time / overlay  
            haxis = axis1 vaxis= axis6 name = "rdiffl" noframe nolegend  ;
         run;
         quit;

    
      %end;
      * create a page of graphs, 4 per page = 2 covariates per page ;

     
      

      
       ods pdf body="&gfilename"  ;

       
       %let graphlist =  1:outc  2:rdiffu   2:rdiffl   ;
       %if &frombootstrap = 0 %then %let graphlist = 1:outc 2:rdiffu ;   

       proc greplay igout = GSEG nofs tc=work.tempcat ;
       template= newtemp ;
       treplay &graphlist ;
       run;
       quit;

      

  ods pdf close  ;

  
 %mend ;













