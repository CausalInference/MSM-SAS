
options linesize=132 pagesize=54;
 
*%include "msm.sas";
%include 'msm2_24_2015.sas';
**SAMPLEDATA Data;
data sampledata(drop = i j  ahbp1-ahbp8 abmi1-abmi8);

    call streaminit(5027);
    
    do i=1 to 1000;
        baseage = int( 35 + 25*rand('uniform'));
        
        array ahbp(8);
        array abmi(8);
       


        do j=1 to 8;
            ahbp(j) = (0.4>rand('uniform'));
            abmi(j) = log(25+5*(rand('normal')));
        end;
        
        do j=3 to 8 until (dia or censlost or censdead);
            id=i;
            time=j-3;
            


            hbp     = ahbp(j);
            hbp_l1  = ahbp(j-1);
            hbp_l2  = ahbp(j-2);
            hbp_b   = ahbp(2);

            bmi     = abmi(j);
            bmi_l1  = abmi(j-1);
            bmi_l2  = abmi(j-2);
            bmi_b   = abmi(2);

            dia = ( (j/500) >rand('uniform'));
            censlost  = (0.05>rand('uniform'));
            censdead  = (0.05>rand('uniform'));

            if bmi<=3.22 then A=1;
                         else A=0;
            if bmi_l1 <= 3.22 then A_l1=1;
                              else A_l1=0;
            if time=0 then do;
                if A_l1=1 then do;
                    Apreb=1;
                    A_l1=0;
                end;
                else do;
                    Apreb=0;
                end;
            end;    
                
            hbp0 = ahbp(1); 
            output;

            end;
       end;
   
run;

proc freq data = sampledata (where = (dia = 1)) ;
table time ;
run;

proc means data=sampledata;
title 'Means of SAMPLEDATA data';
run;

* Exampe 1 : Normal MSM with analytic variance (msm_var = 1). The resulting data sets when save_results = 1 are
             contained in the library indicated by the libsurv variable (default = work). The covariance matrix is 
             in the data set betas_covariance and the table of results is contained in the data set contained in the
             results_name variable (default = _results). The results data set and covariance matrix are also in the 
             non-permanents data set _msmsw and varbetas_sw. ;
* In the normal case with msm_var = 1, the variable use_genmod is set to 0. ;
* When the macro is run with msm_var = 0, the output is saved are  from the last proc logistic. When save_results = 1, the 
  saved results are in the data sets results_name and betas_covariance. 
* When running with msm_var = 0 and use_genmod = 1 the robust variance is used.
* When the macro is run with msm_var = 1 the resulting data sets contain the analytic covariance  (beta_covariance) and the 
  resulting table of final results with the confidence bounds (in results_name ). ;
* When running with bootstrap = 1 the macro will set msm_var = 0 and the variance will be based on the bootstrap samples. The 
  calculated covariance matrix is contained in the temporary data set _cov_sw and in the permanent data set beta_covariance, if
  desired.  ;
* When including interactions with covMSMbh and/or AMSM the macro will create the numeric variables internaly and name them 
  covbh_interX and amsm_interX where X is a number of the interaction. In the call below, when including inter_MSM=hbp_b, the
  created variable will be amsm_inter1. This is the name of the variable that will be used in the weighted model and also 
  when using the survival_cures = 1 option. This is to ensure that the definitions are consistent in the weighted models and
  the calculation of the survival probabilities.;

options nonotes nomprint ;
options notes mprint  ;
%msm(
    data= sampledata,
    id=id,
    time=time,
    time_knots = 0 1 2 3 4 5,
     
    outcMSM = dia,
    AMSM = A,
    covMSMbh = time ,
    classMSMbh = time ,
    contMSMbh = 0,
    covMSMbv = baseage hbp_b bmi_b, 
    inter_MSM =   hbp_b    ,
    inter_time = A baseage  ,
    msm_var = 1,
     
    A_censoring = A ,
    cens  = censlost,
    covCd =  time baseage hbp_b bmi_b hbp  , 
    classCd = time ,
     
    covCn = time baseage hbp_b bmi_b ,
    classCn = time ,
   

    
    A= A,
    covAd=time baseage hbp_b bmi_b   hbp bmi_l1 bmi_l2,
    classAd = time ,
    covAn = time baseage hbp_b bmi_b ,
    classAn = time ,
    
    use_genmod = 0 ,

    bootstrap = 1 ,
    nboot = 25,
    bootstart = 0,
    bootend = 25,

   

    save_results = 0 ,
    save_weights = 0 ,
    libsurv = work,
    results_name = _results ,
    weights_name = _weights ,

   survival_curves = 1,
   /* use_natural_course = 0,
    competing_risk = censdead , */
    time_start = 0,
    time_end = 5 ,



    debug = 0
  
    
    );



* Example 2 uses the same models as in example 1, but now calculates the probabily of
  survival using a competing risk of censdead. Here we also calculate the natural course and
  calculate the risk difference to be the default difference of the cumulative incidence between never 
  treated and alsways treated. ;
* The order of the desired output is that the probabilities are first and then come the results for the
  weighted model for outcMSM. ;
* When a single sample is run (bootstrap = 0) the results are contained in the temporary data set _onegraph_. ;
 
options nonotes nomprint ;
*options notes mprint ;
 *%msm(
    data= sampledata,
    id=id,
    time=time,
    time_knots = 0 1 2 3 4 5 ,
     
    outcMSM = dia,
    AMSM = A,
    covMSMbh = time ,
    classMSMbh = time ,
    contMSMbh = 0 ,
    covMSMbv = baseage hbp_b bmi_b, 

    msm_var = 1,
    A_censoring = A ,
    cens  = censlost,
    covCd =  time baseage hbp_b bmi_b hbp  , 
    classCd = time ,
    
    covCn = time baseage hbp_b bmi_b ,
    classCn = time ,
       
    A= A,
    covAd=time baseage hbp_b bmi_b   hbp bmi_l1 bmi_l2,
    classAd = time ,
    covAn = time baseage hbp_b bmi_b ,
    classAn = time ,

    survival_curves = 1,
    use_natural_course = 0,
    competing_risk = censdead ,
    time_start = 0,
    time_end = 5 ,

    bootstrap = 1 ,
    nboot = 10,
    bootstart = 0 ,
    bootend = 10 ,

    save_results = 0,
    save_weights = 0,
    libsurv = work,
    results_name = _results,
    weights_name = _weights

    

    
    );

* Advanced use of msm macro using amsm_type = 3 to compare two general
  regimes of always treated verses delaying treatment until the third period;
* To perform this comparison the user must construct the two regimes by modifying the
  two sub-macros amsm_user_def_macro0 and amsm_user_def_macro1 as below ;
* Example 3 uses the same models as the previous examples, but now compares two general 
  regimes and calculates the bootstrap variance using 25 bootstrap samples. ;
* In the bootstrap analysis several data sets are created. Some are permement  and some are
  not. The permenent data sets contain the bootstrap results for each type of analysis that was
  performed : The permenent data sets are set in the call statement to hold the individual bootstrap 
  results for the weighted model for outcMSM (sw), competing risk (cr), and natural course (nc), results 
  each type of "survival probability" ( Kaplan-Meyer type curves ,km ), cumulative incidence (ci), and 
  natural course (nc) for both cases of always treated and never treated. These are in the _surv data sets.;
* The bootstrap convariance matrices are contained in the data sets _cov_sw, _cov_nc, and _cov_cr. ;
* The output tables included at the end of the program run are contained in the data sets _msmsw (outcMSM), 
  _msmcr (competing risk), and _msmnc (natural course). ;
* The dataset _bgraphsr contains the tables for the bootstrap results for the various survival probabilites.  ;

%macro amsm_user_def_macro0 ;
  * delay starting till after 2 period ;
  A = 0 ;
  if &time > 2 then A = 1 ;

%mend ;

%macro amsm_user_def_macro1 ;
 * start at baseline ;
   A = 1;
%mend ;

options nonotes nomprint ;
  *  %msm(
    data= sampledata,
    id=id,
    time=time,
    time_knots = 0 1 2 3 4 5 ,
     
    outcMSM = dia,
    AMSM = A,
    amsm_type = 3 ,
    covMSMbh = time ,
    classMSMbh = time ,
    contMSMbh = 0 ,
    covMSMbv = baseage hbp_b bmi_b, 

    msm_var = 1,
    A_censoring = A ,
    cens  = censlost,
    covCd =  time baseage hbp_b bmi_b hbp  , 
    classCd = time ,
    
    covCn = time baseage hbp_b bmi_b ,
    classCn = time ,
    

    
    A= A,
    covAd=time baseage hbp_b bmi_b   hbp bmi_l1 bmi_l2,
    classAd = time ,
    covAn = time baseage hbp_b bmi_b ,
    classAn = time ,

    survival_curves = 1,
    use_natural_course = 1,
    competing_risk = censdead ,
    time_start = 0,
    time_end = 5 ,

    bootstrap = 1,
    nboot = 25,
    bootstart = 0,
    bootend = 25 ,

    
    save_results = 0,
    save_weights = 0,
    libsurv = work,
    results_name = _results,
    weights_name = _weights

    
    );

* Example 4 : Running a bootstrap using two separate calls and putting the results back together. 
  In some situations the time needed to run a set of bootstraps may be too long and it may be desired
  to run the analysis in separate parts at the same time. This can be accomplished with modifying the 
  input parameters bootstart and bootend while keeping nboot to be the same for both runs. In running a
  bootstrap analysis the macro will first try to create a data set with the ids that will contribute to 
  each of the nboot samples. These ids are only determined by the nboot option so that the same collection 
  samples would be obtained by running all nboot samples in a single call to the macro or by splitting up
  the work into separate parts. Below we split the 25 samples into two runs. The first will run the samples 
  from 0 to 12 and the second will run the samples 13 to 25. The last step is to call the bootstrap_results
  macro to put the individual results back together to obtain the various output tables for the complete 25
  samples.

  Note: Due to the way the macro works for creating the output, the portion of the bootstrap runs that 
        contains sample 0 (original data), the msm macro will create output tables for this shorter run. All
        other sets of bootstrap samples will not contain any output.
    ;

option nonotes nomprint ;
      *%msm(
    data= sampledata,
    id=id,
    time=time,
    time_knots = 0 1 2 3 4 5 ,
     
    outcMSM = dia,
    AMSM = A,
    amsm_type = 3 ,
    covMSMbh = time ,
    classMSMbh = time ,
    contMSMbh = 0 ,
    covMSMbv = baseage hbp_b bmi_b, 

    msm_var = 1,
     
    A_censoring = A ,
    cens  = censlost,
    covCd =  time baseage hbp_b bmi_b hbp  , 
    classCd = time ,    
    covCn = time baseage hbp_b bmi_b ,
    classCn = time ,
    

    
    A= A,
    covAd=time baseage hbp_b bmi_b   hbp bmi_l1 bmi_l2,
    classAd = time ,
    covAn = time baseage hbp_b bmi_b ,
    classAn = time ,

    survival_curves = 1,
    use_natural_course = 1,
    competing_risk = censdead ,
    time_start = 0,
    time_end = 5 ,

    bootstrap = 1,
    nboot = 25,
    bootstart = 0,
    bootend = 12 
    
    );

     *%msm(
    data= sampledata,
    id=id,
    time=time,
    time_knots = 0 1 2 3 4 5 ,
     
    outcMSM = dia,
    AMSM = A,
    amsm_type = 3 ,
    covMSMbh = time ,
    classMSMbh = time ,
    contMSMbh = 0 ,
    covMSMbv = baseage hbp_b bmi_b, 

    msm_var = 1,
     
    A_censoring = A ,
    cens  = censlost,
    covCd =  time baseage hbp_b bmi_b hbp  , 
    classCd = time ,
    
    covCn = time baseage hbp_b bmi_b ,
    classCn = time ,
        
    A= A,
    covAd=time baseage hbp_b bmi_b   hbp bmi_l1 bmi_l2,
    classAd = time ,
    covAn = time baseage hbp_b bmi_b ,
    classAn = time ,

    survival_curves = 1,
    use_natural_course = 1,
    competing_risk = censdead ,
    time_start = 0,
    time_end = 5 ,

    bootstrap = 1,
    nboot = 25,
    bootstart = 13,
    bootend = 25 
    
    );



    * Call bootstrap_results macro to create tables of final output. Some of the variables used need
      to match what was used in the original calls to the msm macro. ;
    * The variable modeltype can be one of : sw, nc, cr, surv for the four different tables for the
      bootstrap results. In the call below, we are interested in the results when the survival_curves
      option was used. ;
    * The variables are defined in the documentation for the msm macro. Just note that since there were
      two runs the variables numparts and samplestart need to be defined correctly to determine the names
      of the data sets to put together to obtain a data set containing all nboot + 1 samples.;

   * %bootstrap_results(
       bootlib = work ,
       bootname= bootstrap ,
       modeltype = surv ,
       numparts = 2,
       samplestart = 0 13,
       numboot = 25,
       time = time ,
       use_natural_course = 1,
       outcMSM = dia,
       competing_risk = censdead,
       riskdiff1 = ci_treat ,
       riskdiff0 = ci_nontreat ,
       print_ci = 0 
       );

