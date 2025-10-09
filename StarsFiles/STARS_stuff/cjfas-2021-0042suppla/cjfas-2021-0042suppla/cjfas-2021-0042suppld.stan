functions{
  
  row_vector arrival_probs(int T, real mu, real sigma){
    vector[T] base;
    vector[T] tease = cumulative_sum(rep_vector(1, T)) - 0.5; // the t's

    base = -square(log(tease) - mu) / (2 * square(sigma)) - log(tease);
    
    return(softmax(base)');
  }  
  
  row_vector arrival_probs_re(int T, real mu, real sigma, vector eps){
    vector[T] base;
    vector[T] tease = cumulative_sum(rep_vector(1, T)) - 0.5; // the t's

    base = -square(log(tease) - mu) / (2 * square(sigma)) + eps - log(tease);

    return(softmax(base)');
  }
  
  int first_capture(int[] x_i) {
    for (k in 1:size(x_i))
      if (x_i[k] > 0)
        return k;
      return 0;
  }
  
  int next_capture(int[] x_i, int now, int last) {
    
    for (k in (now + 1):last)
      if (x_i[k] > 0)
        return k;
    
    return 0;
  }
  
  int last_capture(int[] x_i) {
    for (k_rev in 0:(size(x_i) - 1)) {
      int k;
      k =  size(x_i) - k_rev;
      if (x_i[k] > 0)
        return k;
    }
    return 0;
  }
  
}

data {
  int<lower = 0>               N_14;
  int<lower = 0>               N_15;
  int<lower = 0>               N_16;
  int<lower = 0>               N_17;
  int<lower = 0>               N_18;
  
  int<lower = 1>               max_T; //maximum observed T
  int<lower = 1>               T;

  int<lower = 1, upper = 2>    rls_grp_15[N_15];
  int<lower = 1, upper = T>    rls_day_15[N_15];
  int<lower = 1, upper = 2>    rls_grp_16[N_16];
  int<lower = 1, upper = T>    rls_day_16[N_16];  

  int<lower = 1, upper = 2>    rls_grp_18[N_18];
  int<lower = 1, upper = T>    rls_day_18[N_18];
  int                          ch_14[2, N_14, 11]; //capture history  
  int                          ch_15[2, N_15, 11]; //capture history
  int                          ch_16[2, N_16, 11]; //capture history
  int                          ch_17[2, N_17, 11]; //capture history
  int                          ch_18[2, N_18, 11]; //capture history
  
  matrix[T, 5]                 KL;  //KL Flow - Fremont WEir
  matrix[T, 5]                 SAC; //Freeport Flwo
  matrix[T, 5]                 RIO;  //rio vista flow
  matrix[T, 5]                 TMP;  //Freeport temp, proxy for all temps except Yolo
  matrix[T, 5]                 EI;  //import export ratio by
  
  matrix[T, 3]                 YOL;  //Yolo flow
  matrix[T, 3]                 YOL_O;  //Yolo overtop indicator
  matrix[T, 3]                 YOL_R;  //Stage height above Fremont Weir
  matrix[T, 3]                 YOL_T;  //Yolo temperatuire

  vector[T]                    tag_surv_prob_14;  
  vector[T]                    tag_surv_prob_15;
  vector[T]                    tag_surv_prob_16;
  vector[T]                    tag_surv_prob_17;
  vector[T]                    tag_surv_prob_18[2];
  
  real<lower = 0>              dist[5];
  
}

// If you're actually reading this code, congratulations you've have exceeded my expectations
// So let's be honest - this code is ugly, it's ugly as hell
// But it works
// There is probably a clever way to code this using 4 dimensional arrays, but there is 
// a ragged array problem
// My solution was to hardcode everything with named variables, it's kludgey and fragile
// and entirely dependent on the structure of the data I have
// It is what it is. Please feel free to email me and tell me how 
// absolutely horrified and embarassed you are for me: dhance@usgs.gov

transformed data{
  //Transform capture histories into M matrices and L vectors
  //Only capture history pairs observed in the data are calculated
  
  real   N_rls_grp_15[2];
  real   N_rls_grp_16[2];
  real   N_rls_grp_18[2];
  
  //M Matrices 2014---------------------------------------------------------------------------\\
  
  row_vector[T]       m_t0s1_t1s1_14 = rep_row_vector(0, T); //n first detected knight's landing
  row_vector[T]       m_t0s1_t4s1_14 = rep_row_vector(0, T); //n first detected feather
  row_vector[T]       m_t0s1_t5s1_14 = rep_row_vector(0, T); //n first detected sacramento

  matrix[T, T]        m_t1s1_t4s1_14 = rep_matrix(0, T, T); //n knight's landing -> feather
  matrix[T, T]        m_t1s1_t5s1_14 = rep_matrix(0, T, T); //n knight's landing -> sacramento

  matrix[T, T]        m_t4s1_t5s1_14 = rep_matrix(0, T, T); //n feather -> sacramento

  matrix[T, T]        m_t5s1_t6s1_14 = rep_matrix(0, T, T); //n sacramento -> freeport
  matrix[T, T]        m_t5s1_t7s1_14 = rep_matrix(0, T, T); //n sacramento -> sac blw sutter/steam
  matrix[T, T]        m_t5s1_t7s3_14 = rep_matrix(0, T, T); //n sacramento -> sutter
  matrix[T, T]        m_t5s1_t7s4_14 = rep_matrix(0, T, T); //n sacramento -> steam
  matrix[T, T]        m_t5s1_t9s1_14 = rep_matrix(0, T, T); //n sacramento -> rio vista

  matrix[T, T]        m_t6s1_t7s1_14 = rep_matrix(0, T, T); //n freeport -> sac blw sutter/steam
  matrix[T, T]        m_t6s1_t7s3_14 = rep_matrix(0, T, T); //n freeport -> sutter
  matrix[T, T]        m_t6s1_t7s4_14 = rep_matrix(0, T, T); //n freeport -> steam
  matrix[T, T]        m_t6s1_t9s1_14 = rep_matrix(0, T, T); //n freeport -> rio vista

  matrix[T, T]        m_t7s1_t8s1_14 = rep_matrix(0, T, T); //n sac blw sutter/steam -> sac blw georgiana
  matrix[T, T]        m_t7s1_t8s5_14 = rep_matrix(0, T, T); //n sac blw sutter/steam -> georgiana
  matrix[T, T]        m_t7s1_t9s1_14 = rep_matrix(0, T, T); //n sac blw sutter/steam -> rio vista

  matrix[T, T]        m_t7s3_t9s1_14 = rep_matrix(0, T, T); //n sutter -> rio vista
  matrix[T, T]        m_t7s4_t9s1_14 = rep_matrix(0, T, T); //n steam -> rio vista

  matrix[T, T]        m_t8s1_t9s1_14 = rep_matrix(0, T, T); //n sac blw georgiana -> rio vista
  matrix[T, T]        m_t8s1_t11s1_14 = rep_matrix(0, T, T); //n sac blw georgiana -> benecia

  matrix[T, T]        m_t9s1_t11s1_14 = rep_matrix(0, T, T); //n rio vista -> benecia
  
  //M Matrices 2015-----------------------------------------------------------------\\
  matrix[2, T]        m_t0s1_t1s1_15 = rep_matrix(0, 2, T); //n first detected knight's landing
  matrix[2, T]        m_t0s1_t4s1_15 = rep_matrix(0, 2, T); //n first detected feather
  matrix[2, T]        m_t0s1_t5s1_15 = rep_matrix(0, 2, T); //n first detected sacramento
  matrix[2, T]        m_t0s1_t7s1_15 = rep_matrix(0, 2, T); //n first detected sac blw sutter/steam
  matrix[2, T]        m_t0s1_t7s3_15 = rep_matrix(0, 2, T); //n first detected sutter
  matrix[2, T]        m_t0s1_t7s4_15 = rep_matrix(0, 2, T); //n first detected steam

  matrix[T, T]        m_t1s1_t4s1_15 = rep_matrix(0, T, T); //n knight's landing -> feather
  matrix[T, T]        m_t1s1_t5s1_15 = rep_matrix(0, T, T); //n knight's landing -> sacramento
  matrix[T, T]        m_t1s1_t7s1_15 = rep_matrix(0, T, T); //n knight's landing -> sac blw sutter/steam
  matrix[T, T]        m_t1s1_t7s3_15 = rep_matrix(0, T, T); //n knight's landing ->sutter
  matrix[T, T]        m_t1s1_t7s4_15 = rep_matrix(0, T, T); //n knight's landing -> steam
  matrix[T, T]        m_t1s1_t8s1_15 = rep_matrix(0, T, T); //n knight's landing -> sac blw georgiana

  matrix[T, T]        m_t4s1_t5s1_15 = rep_matrix(0, T, T); //n feather -> sacramento
  matrix[T, T]        m_t4s1_t6s1_15 = rep_matrix(0, T, T); //n feather -> freeport
  matrix[T, T]        m_t4s1_t7s1_15 = rep_matrix(0, T, T); //n feather -> sac blw sutter/steam
  matrix[T, T]        m_t4s1_t7s3_15 = rep_matrix(0, T, T); //n feather -> sutter
  matrix[T, T]        m_t4s1_t7s4_15 = rep_matrix(0, T, T); //n feather -> steam
  matrix[T, T]        m_t4s1_t8s1_15 = rep_matrix(0, T, T); //n feather -> sac blw georgiana
  matrix[T, T]        m_t4s1_t8s5_15 = rep_matrix(0, T, T); //n feather -> georgiana
  matrix[T, T]        m_t4s1_t9s1_15 = rep_matrix(0, T, T); //n feather -> rio vista

  matrix[T, T]        m_t5s1_t6s1_15 = rep_matrix(0, T, T); //n sacramento -> freeport
  matrix[T, T]        m_t5s1_t7s1_15 = rep_matrix(0, T, T); //n sacramento -> sac blw sutter/steam
  matrix[T, T]        m_t5s1_t7s3_15 = rep_matrix(0, T, T); //n sacramento -> sutter
  matrix[T, T]        m_t5s1_t7s4_15 = rep_matrix(0, T, T); //n sacramento -> steam
  matrix[T, T]        m_t5s1_t8s1_15 = rep_matrix(0, T, T); //n sacramento -> sac blw georgiana
  matrix[T, T]        m_t5s1_t8s5_15 = rep_matrix(0, T, T); //n sacramento -> georgiana
  matrix[T, T]        m_t5s1_t9s1_15 = rep_matrix(0, T, T); //n sacramento -> rio vista

  matrix[T, T]        m_t6s1_t7s1_15 = rep_matrix(0, T, T); //n freeport -> sac blw sutter/steam
  matrix[T, T]        m_t6s1_t7s3_15 = rep_matrix(0, T, T); //n freeport -> sutter
  matrix[T, T]        m_t6s1_t7s4_15 = rep_matrix(0, T, T); //n freeport -> steam

  matrix[T, T]        m_t7s1_t8s1_15 = rep_matrix(0, T, T); //n sac blw sutter/steam -> sac blw georgiana
  matrix[T, T]        m_t7s1_t8s5_15 = rep_matrix(0, T, T); //n sac blw sutter/steam -> georgiana
  matrix[T, T]        m_t7s1_t9s1_15 = rep_matrix(0, T, T); //n sac blw sutter/steam -> rio vista

  matrix[T, T]        m_t7s3_t9s1_15 = rep_matrix(0, T, T); //n sutter -> rio vista
  matrix[T, T]        m_t7s3_t11s1_15 = rep_matrix(0, T, T); //n sutter -> benecia
  matrix[T, T]        m_t7s4_t9s1_15 = rep_matrix(0, T, T); //n steam -> rio vista
  matrix[T, T]        m_t7s4_t11s1_15 = rep_matrix(0, T, T); //n steam -> benecia

  matrix[T, T]        m_t8s1_t9s1_15 = rep_matrix(0, T, T); //n sac blw georgiana -> rio vista
  matrix[T, T]        m_t8s1_t11s1_15 = rep_matrix(0, T, T); //n sac blw georgiana -> benecia
  matrix[T, T]        m_t8s5_t9s5_15 = rep_matrix(0, T, T); //n georgiana -> interior delta
  matrix[T, T]        m_t8s5_t11s1_15 = rep_matrix(0, T, T); //n georgiana -> benecia

  matrix[T, T]        m_t9s1_t11s1_15 = rep_matrix(0, T, T); //n rio vista -> benecia
  matrix[T, T]        m_t9s5_t11s1_15 = rep_matrix(0, T, T); //n interior delta -> benecia
  
  //M Matrices 2016-----------------------------------------------------------------\\
  matrix[2, T]        m_t0s1_t1s1_16 = rep_matrix(0, 2, T); //n first detected knight's landing
  matrix[2, T]        m_t0s1_t4s1_16 = rep_matrix(0, 2, T); //n first detected feather
  matrix[2, T]        m_t0s1_t4s2_16 = rep_matrix(0, 2, T); //n first detected toe drain
  matrix[2, T]        m_t0s1_t5s1_16 = rep_matrix(0, 2, T); //n first detected sacramento
  matrix[2, T]        m_t0s1_t7s1_16 = rep_matrix(0, 2, T); //n first detected sac blw sutter/steam
  matrix[2, T]        m_t0s1_t10s1_16 = rep_matrix(0, 2, T); //n first detected chipps

  matrix[T, T]        m_t1s1_t4s1_16 = rep_matrix(0, T, T); //n knight's landing -> feather
  matrix[T, T]        m_t1s1_t4s2_16 = rep_matrix(0, T, T); //n knight's landing -> toe drain
  matrix[T, T]        m_t1s1_t5s1_16 = rep_matrix(0, T, T); //n knight's landing -> sacramento
  matrix[T, T]        m_t1s1_t7s4_16 = rep_matrix(0, T, T); //n knight's landing -> steam

  matrix[T, T]        m_t4s1_t5s1_16 = rep_matrix(0, T, T); //n feather -> sacramento
  matrix[T, T]        m_t4s1_t6s1_16 = rep_matrix(0, T, T); //n feather -> freeport
  matrix[T, T]        m_t4s1_t7s1_16 = rep_matrix(0, T, T); //n feather -> sac blw sutter/steam

  matrix[T, T]        m_t4s2_t9s1_16 = rep_matrix(0, T, T); //n toe drain -> rio vista

  matrix[T, T]        m_t5s1_t6s1_16 = rep_matrix(0, T, T); //n sacramento -> freeport
  matrix[T, T]        m_t5s1_t7s4_16 = rep_matrix(0, T, T); //n sacramento -> steam
  matrix[T, T]        m_t5s1_t9s1_16 = rep_matrix(0, T, T); //n sacramento -> rio vista

  matrix[T, T]        m_t6s1_t7s1_16 = rep_matrix(0, T, T); //n freeport -> sac blw sutter/steam
  matrix[T, T]        m_t6s1_t7s3_16 = rep_matrix(0, T, T); //n freeport -> sutter
  matrix[T, T]        m_t6s1_t7s4_16 = rep_matrix(0, T, T); //n freeport -> steam
  matrix[T, T]        m_t6s1_t9s1_16 = rep_matrix(0, T, T); //n freeport -> rio vista
  
  matrix[T, T]        m_t7s1_t9s1_16 = rep_matrix(0, T, T); //n sac blw sutter/steam -> rio vista
  matrix[T, T]        m_t7s1_t9s5_16 = rep_matrix(0, T, T); //n sac blw sutter/steam -> interior delta
  matrix[T, T]        m_t7s1_t10s1_16 = rep_matrix(0, T, T); //n sac blw sutter/steam -> chipps

  matrix[T, T]        m_t7s3_t9s1_16 = rep_matrix(0, T, T); //n sutter -> rio vista
  matrix[T, T]        m_t7s3_t10s1_16 = rep_matrix(0, T, T); //n sutter -> chipps

  matrix[T, T]        m_t7s4_t9s1_16 = rep_matrix(0, T, T); //n steam -> rio vista
  matrix[T, T]        m_t7s4_t10s1_16 = rep_matrix(0, T, T); //n steam -> chipps

  matrix[T, T]        m_t9s1_t10s1_16 = rep_matrix(0, T, T); //n rio vista -> chipps
  matrix[T, T]        m_t9s1_t11s1_16 = rep_matrix(0, T, T); //n rio vista -> benecia

  matrix[T, T]        m_t9s5_t10s1_16 = rep_matrix(0, T, T); //n interior delta -> chipps
  matrix[T, T]        m_t9s5_t11s1_16 = rep_matrix(0, T, T); //n interior delta -> benecia

  matrix[T, T]        m_t10s1_t11s1_16 = rep_matrix(0, T, T); //n chipps -> benecia
  
  //M Matrices 2017------------------------------------------------------------------\\
  row_vector[T]       m_t0s1_t1s1_17 = rep_row_vector(0, T); //n first detected knight's landing
  row_vector[T]       m_t0s1_t3s1_17 = rep_row_vector(0, T); //n first detected blw fremont weir
  row_vector[T]       m_t0s1_t4s1_17 = rep_row_vector(0, T); //n first detected feather
  row_vector[T]       m_t0s1_t5s1_17 = rep_row_vector(0, T); //n first detected sacramento
  row_vector[T]       m_t0s1_t6s1_17 = rep_row_vector(0, T); //n first detected freeport
  row_vector[T]       m_t0s1_t7s1_17 = rep_row_vector(0, T); //n first detected sac blw sutter/steam
  row_vector[T]       m_t0s1_t8s1_17 = rep_row_vector(0, T); //n first detected sac blw georgiana
  row_vector[T]       m_t0s1_t9s1_17 = rep_row_vector(0, T); //n first detected rio vista
  row_vector[T]       m_t0s1_t10s1_17 = rep_row_vector(0, T); //n first detected chipps

  matrix[T, T]        m_t1s1_t2s1_17 = rep_matrix(0, T, T); //n knight's landing -> abv fremont weir
  matrix[T, T]        m_t1s1_t3s1_17 = rep_matrix(0, T, T); //n knight's landing -> blw fremont weir
  matrix[T, T]        m_t1s1_t4s2_17 = rep_matrix(0, T, T); //n knight's landing -> toe drain
  matrix[T, T]        m_t1s1_t9s1_17 = rep_matrix(0, T, T); //n knight's landing -> rio vista

  matrix[T, T]        m_t2s1_t3s1_17 = rep_matrix(0, T, T); //n abv fremont weir -> blw fremont weir
  matrix[T, T]        m_t2s1_t4s1_17 = rep_matrix(0, T, T); //n abv fremont weir -> feather
  matrix[T, T]        m_t2s1_t4s2_17 = rep_matrix(0, T, T); //n abv fremont weir -> toe drain
  matrix[T, T]        m_t2s1_t5s1_17 = rep_matrix(0, T, T); //n abv fremont weir -> sacramento
  matrix[T, T]        m_t2s1_t9s1_17 = rep_matrix(0, T, T); //n abv fremont weir -> rio vista
  matrix[T, T]        m_t2s1_t10s1_17 = rep_matrix(0, T, T); //n abv fremont weir -> chipps

  matrix[T, T]        m_t3s1_t4s1_17 = rep_matrix(0, T, T); //n blw fremont weir -> feather
  matrix[T, T]        m_t3s1_t5s1_17 = rep_matrix(0, T, T); //n blw fremont weir -> sacramento
  matrix[T, T]        m_t3s1_t6s1_17 = rep_matrix(0, T, T); //n blw fremont weir -> freeport
  matrix[T, T]        m_t3s1_t7s1_17 = rep_matrix(0, T, T); //n blw fremont weir -> sac blw sutter/steam
  matrix[T, T]        m_t3s1_t7s3_17 = rep_matrix(0, T, T); //n blw fremont weir -> sutter
  matrix[T, T]        m_t3s1_t7s4_17 = rep_matrix(0, T, T); //n blw fremont weir -> steam
  matrix[T, T]        m_t3s1_t8s1_17 = rep_matrix(0, T, T); //n blw fremont weir -> sac blw georgiana
  matrix[T, T]        m_t3s1_t9s1_17 = rep_matrix(0, T, T); //n blw fremont weir -> rio vista
  matrix[T, T]        m_t3s1_t9s5_17 = rep_matrix(0, T, T); //n blw fremont weir -> interior delta
  matrix[T, T]        m_t3s1_t10s1_17 = rep_matrix(0, T, T); //n blw fremont weir -> chipps

  matrix[T, T]        m_t4s1_t5s1_17 = rep_matrix(0, T, T); //n feather -> sacramento
  matrix[T, T]        m_t4s1_t6s1_17 = rep_matrix(0, T, T); //n feather -> freeport
  matrix[T, T]        m_t4s1_t7s1_17 = rep_matrix(0, T, T); //n feather -> sac blw sutter/steam

  matrix[T, T]        m_t4s2_t9s1_17 = rep_matrix(0, T, T); //n toe drain -> rio vista
  matrix[T, T]        m_t4s2_t10s1_17 = rep_matrix(0, T, T); //n toe drain -> chipps
  matrix[T, T]        m_t4s2_t11s1_17 = rep_matrix(0, T, T); //n toe drain -> benecia

  matrix[T, T]        m_t5s1_t6s1_17 = rep_matrix(0, T, T); //n sacramento -> freeport
  matrix[T, T]        m_t5s1_t7s1_17 = rep_matrix(0, T, T); //n sacramento -> sac blw sutter/steam
  matrix[T, T]        m_t5s1_t7s3_17 = rep_matrix(0, T, T); //n sacramento -> sutter
  matrix[T, T]        m_t5s1_t7s4_17 = rep_matrix(0, T, T); //n sacramento -> steam
  matrix[T, T]        m_t5s1_t8s1_17 = rep_matrix(0, T, T); //n sacramento -> sac blw georgiana
  matrix[T, T]        m_t5s1_t9s1_17 = rep_matrix(0, T, T); //n sacramento -> rio vista
  matrix[T, T]        m_t5s1_t9s5_17 = rep_matrix(0, T, T); //n sacramento -> interior delta
  matrix[T, T]        m_t5s1_t10s1_17 = rep_matrix(0, T, T); //n sacramento -> chipps

  matrix[T, T]        m_t6s1_t7s1_17 = rep_matrix(0, T, T); //n freeport -> sac blw sutter/steam
  matrix[T, T]        m_t6s1_t7s3_17 = rep_matrix(0, T, T); //n freeport -> sutter
  matrix[T, T]        m_t6s1_t7s4_17 = rep_matrix(0, T, T); //n freeport -> steam
  matrix[T, T]        m_t6s1_t8s1_17 = rep_matrix(0, T, T); //n freeport -> sac blw georgiana
  matrix[T, T]        m_t6s1_t9s1_17 = rep_matrix(0, T, T); //n freeport -> rio vista
  matrix[T, T]        m_t6s1_t10s1_17 = rep_matrix(0, T, T); //n freeport -> chipps

  matrix[T, T]        m_t7s1_t8s1_17 = rep_matrix(0, T, T); //n sac blw sutter/steam -> sac blw georgiana
  matrix[T, T]        m_t7s1_t9s1_17 = rep_matrix(0, T, T); //n sac blw sutter/steam -> rio vista
  matrix[T, T]        m_t7s1_t9s5_17 = rep_matrix(0, T, T); //n sac blw sutter/steam -> interior delta
  matrix[T, T]        m_t7s1_t10s1_17 = rep_matrix(0, T, T); //n sac blw sutter/steam -> chipps

  matrix[T, T]        m_t7s3_t9s1_17 = rep_matrix(0, T, T); //n sutter -> rio vista
  matrix[T, T]        m_t7s3_t10s1_17 = rep_matrix(0, T, T); //n sutter -> chipps
  matrix[T, T]        m_t7s3_t11s1_17 = rep_matrix(0, T, T); //n sutter -> benecia

  matrix[T, T]        m_t7s4_t9s1_17 = rep_matrix(0, T, T); //n steam -> rio vista
  matrix[T, T]        m_t7s4_t10s1_17 = rep_matrix(0, T, T); //n steam -> chipps
  matrix[T, T]        m_t7s4_t11s1_17 = rep_matrix(0, T, T); //n steam -> benecia

  matrix[T, T]        m_t8s1_t9s1_17 = rep_matrix(0, T, T); //n sac blw georgiana -> rio vista
  matrix[T, T]        m_t8s1_t10s1_17 = rep_matrix(0, T, T); //n sac blw georgiana -> chipps

  matrix[T, T]        m_t9s1_t10s1_17 = rep_matrix(0, T, T); //n rio vista -> chipps
  matrix[T, T]        m_t9s1_t11s1_17 = rep_matrix(0, T, T); //n rio vista -> benecia

  matrix[T, T]        m_t9s5_t10s1_17 = rep_matrix(0, T, T); //n interior delta -> chipps

  matrix[T, T]        m_t10s1_t11s1_17 = rep_matrix(0, T, T); //n chipps -> benecia
 //M Matrices 2018-------------------------------------------------------------------\\
  row_vector[T]       m_t0s1_t1s1_18[2]; //n first detected knight's landing

  matrix[T, T]        m_t1s1_t2s1_18[2]; //n knight's landing -> abv fremont weir
  matrix[T, T]        m_t1s1_t3s1_18[2]; //n knight's landing -> blw fremont weir
  matrix[T, T]        m_t1s1_t4s1_18[2]; //n knight's landing -> feather
  matrix[T, T]        m_t1s1_t5s1_18[2]; //n knight's landing -> sacramento
  matrix[T, T]        m_t1s1_t6s1_18[2]; //n knight's landing -> freeport
  matrix[T, T]        m_t1s1_t7s1_18[2]; //n knight's landing -> sac blw sutter/steam

  matrix[T, T]        m_t2s1_t3s1_18[2]; //n abv fremont weir -> blw fremont weir
  matrix[T, T]        m_t2s1_t4s1_18[2]; //n abv fremont weir -> feather
  matrix[T, T]        m_t2s1_t5s1_18[2]; //n abv fremont weir -> sacramento
  matrix[T, T]        m_t2s1_t6s1_18[2]; //n abv fremont weir -> freeport
  matrix[T, T]        m_t2s1_t7s1_18[2]; //n abv fremont weir -> sac blw sutter/steam
  matrix[T, T]        m_t2s1_t7s3_18[2]; //n abv fremont weir -> sutter
  matrix[T, T]        m_t2s1_t7s4_18[2]; //n abv fremont weir -> steam
  matrix[T, T]        m_t2s1_t9s1_18[2]; //n abv fremont weir -> rio vista

  matrix[T, T]        m_t3s1_t4s1_18[2]; //n blw fremont weir -> feather
  matrix[T, T]        m_t3s1_t5s1_18[2]; //n blw fremont weir -> sacramento
  matrix[T, T]        m_t3s1_t6s1_18[2]; //n blw fremont weir -> freeport
  matrix[T, T]        m_t3s1_t7s1_18[2]; //n blw fremont weir -> sac blw sutter/steam
  matrix[T, T]        m_t3s1_t9s1_18[2]; //n blw fremont weir -> rio vista

  matrix[T, T]        m_t4s1_t5s1_18[2]; //n feather -> sacramento
  matrix[T, T]        m_t4s1_t6s1_18[2]; //n feather -> freeport
  matrix[T, T]        m_t4s1_t7s1_18[2]; //n feather -> sac blw sutter/steam
  matrix[T, T]        m_t4s1_t7s3_18[2]; //n feather -> sutter
  matrix[T, T]        m_t4s1_t7s4_18[2]; //n feather -> steam

  matrix[T, T]        m_t5s1_t6s1_18[2]; //n sacramento -> freeport
  matrix[T, T]        m_t5s1_t7s1_18[2]; //n sacramento -> sac blw sutter/steam
  matrix[T, T]        m_t5s1_t7s3_18[2]; //n sacramento -> sutter
  matrix[T, T]        m_t5s1_t7s4_18[2]; //n sacramento -> steam

  matrix[T, T]        m_t6s1_t7s1_18[2]; //n freeport -> sac blw sutter/steam
  matrix[T, T]        m_t6s1_t7s3_18[2]; //n freeport -> sutter
  matrix[T, T]        m_t6s1_t7s4_18[2]; //n freeport -> steam
  matrix[T, T]        m_t6s1_t8s5_18[2]; //n freeport -> georgiana

  matrix[T, T]        m_t7s1_t8s1_18[2]; //n sac blw sutter/steam -> sac blw georgiana
  matrix[T, T]        m_t7s1_t8s5_18[2]; //n sac blw sutter/steam -> georgiana
  matrix[T, T]        m_t7s1_t9s1_18[2]; //n sac blw sutter/steam -> rio vista
  matrix[T, T]        m_t7s1_t10s1_18[2]; //n sac blw sutter/steam -> chipps

  matrix[T, T]        m_t7s3_t9s1_18[2]; //n sutter -> rio vista
  matrix[T, T]        m_t7s3_t10s1_18[2]; //n sutter -> chipps

  matrix[T, T]        m_t7s4_t9s1_18[2]; //n steam -> rio vista
  matrix[T, T]        m_t7s4_t10s1_18[2]; //n steam -> chipps

  matrix[T, T]        m_t8s1_t9s1_18[2]; //n sac blw georgiana -> rio vista
  matrix[T, T]        m_t8s1_t10s1_18[2]; //n sac blw georgiana -> chipps

  matrix[T, T]        m_t8s5_t9s5_18[2]; //n georgiana -> interior delta

  matrix[T, T]        m_t9s1_t10s1_18[2]; //n rio vista -> chipps
  matrix[T, T]        m_t9s1_t11s1_18[2]; //n rio vista -> benecia

  matrix[T, T]        m_t9s5_t10s1_18[2]; //n interior delta -> chipps

  matrix[T, T]        m_t10s1_t11s1_18[2]; //n chipps -> benecia  
  
  //L Vectors 2014----------------------------------------------------------\\
  vector[T]           l_9_s1_14  = rep_vector(0, T);
  vector[T]           l_9_s5_14  = rep_vector(0, T);
  vector[T]           l_8_s1_14  = rep_vector(0, T);
  vector[T]           l_8_s5_14  = rep_vector(0, T);
  vector[T]           l_7_s1_14  = rep_vector(0, T);
  vector[T]           l_7_s3_14  = rep_vector(0, T);
  vector[T]           l_7_s4_14  = rep_vector(0, T);
  vector[T]           l_6_s1_14  = rep_vector(0, T);
  vector[T]           l_5_s1_14  = rep_vector(0, T);
  vector[T]           l_4_s1_14  = rep_vector(0, T);
  vector[T]           l_1_s1_14  = rep_vector(0, T);
  real                l_0_s1_14  = 0;
  //L Vectors 2015----------------------------------------------------------\\
  vector[T]           l_9_s1_15  = rep_vector(0, T);
  vector[T]           l_9_s5_15  = rep_vector(0, T);
  vector[T]           l_8_s1_15  = rep_vector(0, T);
  vector[T]           l_8_s5_15  = rep_vector(0, T);
  vector[T]           l_7_s1_15  = rep_vector(0, T);
  vector[T]           l_7_s3_15  = rep_vector(0, T);
  vector[T]           l_7_s4_15  = rep_vector(0, T);
  vector[T]           l_6_s1_15  = rep_vector(0, T);
  vector[T]           l_5_s1_15  = rep_vector(0, T);
  vector[T]           l_4_s1_15  = rep_vector(0, T);
  vector[T]           l_1_s1_15  = rep_vector(0, T);
  vector[2]           l_0_s1_15  = rep_vector(0, 2);
  //L Vectors 2016----------------------------------------------------------\\
  vector[T]           l_10_s1_16  = rep_vector(0, T); //n last seen chipps
  vector[T]           l_9_s1_16  = rep_vector(0, T); //n last seen rio vista
  vector[T]           l_9_s5_16  = rep_vector(0, T); //n last seen interior delta
  vector[T]           l_7_s1_16  = rep_vector(0, T); //n last seen sac blw sutter/steam
  vector[T]           l_7_s3_16  = rep_vector(0, T); //n last seen sutter
  vector[T]           l_7_s4_16  = rep_vector(0, T); //n last seen steam
  vector[T]           l_6_s1_16  = rep_vector(0, T); //n last seen freeport
  vector[T]           l_5_s1_16  = rep_vector(0, T); //n last seen sacramento
  vector[T]           l_4_s2_16  = rep_vector(0, T); //n last seen toe drain
  vector[T]           l_4_s1_16  = rep_vector(0, T); //n last seen feather
  vector[T]           l_1_s1_16  = rep_vector(0, T); //n last seen knight's landing
  vector[2]           l_0_s1_16  = rep_vector(0, 2); //n last seen release
  //L Vectors 2017----------------------------------------------------------\\
  vector[T]           l_10_s1_17  = rep_vector(0, T); //n last seen chipps
  vector[T]           l_9_s1_17  = rep_vector(0, T); //n last seen rio vista
  vector[T]           l_9_s5_17  = rep_vector(0, T); //n last seen interior delta
  vector[T]           l_8_s1_17  = rep_vector(0, T); //n last seen sac blw georgiana
  vector[T]           l_8_s5_17  = rep_vector(0, T); //n last seen georgiana
  vector[T]           l_7_s1_17  = rep_vector(0, T); //n last seen sac blw sutter/steam
  vector[T]           l_7_s3_17  = rep_vector(0, T); //n last seen sutter
  vector[T]           l_7_s4_17  = rep_vector(0, T); //n last seen steam
  vector[T]           l_6_s1_17  = rep_vector(0, T); //n last seen freeport
  vector[T]           l_5_s1_17  = rep_vector(0, T); //n last seen sacramento
  vector[T]           l_4_s2_17  = rep_vector(0, T); //n last seen toe drain
  vector[T]           l_4_s1_17  = rep_vector(0, T); //n last seen feather
  vector[T]           l_3_s1_17  = rep_vector(0, T); //n last seen blw fremont weir
  vector[T]           l_2_s1_17  = rep_vector(0, T); //n last seen abv fremont weir
  vector[T]           l_1_s1_17  = rep_vector(0, T); //n last seen knight's landing
  real                l_0_s1_17  = 0; //n last seen release
  //L Vectors 2018----------------------------------------------------------\\
  vector[T]           l_10_s1_18[2]; //n last seen chipps
  vector[T]           l_9_s1_18[2]; //n last seen rio vista
  vector[T]           l_9_s5_18[2]; //n last seen interior delta
  vector[T]           l_8_s1_18[2]; //n last seen sac blw georgiana
  vector[T]           l_8_s5_18[2]; //n last seen georgiana
  vector[T]           l_7_s1_18[2]; //n last seen sac blw sutter/steam
  vector[T]           l_7_s3_18[2]; //n last seen sutter
  vector[T]           l_7_s4_18[2]; //n last seen steam
  vector[T]           l_6_s1_18[2]; //n last seen freeport
  vector[T]           l_5_s1_18[2]; //n last seen sacramento
  vector[T]           l_4_s2_18[2]; //n last seen toe drain
  vector[T]           l_4_s1_18[2]; //n last seen feather
  vector[T]           l_3_s1_18[2]; //n last seen blw fremont weir
  vector[T]           l_2_s1_18[2]; //n last seen abv fremont weir
  vector[T]           l_1_s1_18[2]; //n last seen knight's landing
  real                l_0_s1_18[2]; //n last seen release
  
  //Begin Transformed Data Calculations----------------------------------------------------\\
  N_rls_grp_15[1] = 0;
  N_rls_grp_15[2] = 0;
  N_rls_grp_16[1] = 0;
  N_rls_grp_16[2] = 0;
  N_rls_grp_18[1] = 0;
  N_rls_grp_18[2] = 0;
  
  //Tranformed Data 2014----\\ 
  for (i in 1:N_14){
    int first = first_capture(ch_14[2, i, ]);
    int last = last_capture(ch_14[2, i, ]);
    {
      if (first == 1)
        m_t0s1_t1s1_14[ch_14[1, i, first]] += 1;
      if (first == 4)
        m_t0s1_t4s1_14[ch_14[1, i, first]] += 1;
      if (first == 5)
        m_t0s1_t5s1_14[ ch_14[1, i, first]] += 1;
    }

    if (first > 0 && first != last){
      int now = first;
      int next = next_capture(ch_14[2, i, ], now, last);

      while (next < last){ // This while loop takes us to through penultimate capture 
        if (now == 1){
          if (next == 4)
            m_t1s1_t4s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          if (next == 5)
            m_t1s1_t5s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }

        if (now == 4) {
          if (next == 5)
            m_t4s1_t5s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }

        if (now == 5) {
          if (next == 6)
            m_t5s1_t6s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          if (next == 7){
            if (ch_14[2, i, next] == 1)
              m_t5s1_t7s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
            if (ch_14[2, i, next] == 3)
              m_t5s1_t7s3_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
            if (ch_14[2, i, next] == 4)
              m_t5s1_t7s4_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }

          if (next == 9){
            if (ch_14[2, i, next] == 1)
              m_t5s1_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }
        }

        if (now == 6) {
          if (next == 7){
            if (ch_14[2, i, next] == 1)
              m_t6s1_t7s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
            if (ch_14[2, i, next] == 3)
              m_t6s1_t7s3_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
            if (ch_14[2, i, next] == 4)
              m_t6s1_t7s4_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }

          if (next == 9){
            if (ch_14[2, i, next] == 1)
              m_t6s1_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }
        }

        if (now == 7) {
          if(ch_14[2, i, now] == 1){
            if (next == 8){
              if (ch_14[2, i, next] == 1)
                m_t7s1_t8s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
              if (ch_14[2, i, next] == 5)
                m_t7s1_t8s5_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
            }
            if (next == 9){
              if (ch_14[2, i, next] == 1)
                m_t7s1_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
            }
          }

          if(ch_14[2, i, now] == 3){
            if (next == 9)
                m_t7s3_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }

          if(ch_14[2, i, now] == 4){
            if (next == 9)
                m_t7s4_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }
        }

        if (now == 8) {
          if(ch_14[2, i, now] == 1){
            if (next == 9)
              m_t8s1_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }
        }

        now = next;
        next = next_capture(ch_14[2, i, ], now, last);
      }

      // from second-to-last capture to last capture

      if (now == 1){
        if (next == 4)
          m_t1s1_t4s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        if (next == 5)
          m_t1s1_t5s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
      }

      if (now == 4) {
        if (next == 5)
          m_t4s1_t5s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
      }

      if (now == 5) {
        if (next == 6)
          m_t5s1_t6s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        if (next == 7){
          if (ch_14[2, i, next] == 1)
            m_t5s1_t7s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          if (ch_14[2, i, next] == 3)
            m_t5s1_t7s3_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          if (ch_14[2, i, next] == 4)
            m_t5s1_t7s4_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }
        if (next == 9){
          if (ch_14[2, i, next] == 1)
            m_t5s1_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }
      }

      if (now == 6) {
        if (next == 7){
          if (ch_14[2, i, next] == 1)
            m_t6s1_t7s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          if (ch_14[2, i, next] == 3)
            m_t6s1_t7s3_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          if (ch_14[2, i, next] == 4)
            m_t6s1_t7s4_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }
        if (next == 9){
          if (ch_14[2, i, next] == 1)
            m_t6s1_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }
      }

      if (now == 7) {
        if(ch_14[2, i, now] == 1){
          if (next == 8){
            if (ch_14[2, i, next] == 1)
              m_t7s1_t8s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
            if (ch_14[2, i, next] == 5)
              m_t7s1_t8s5_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_14[2, i, next] == 1)
              m_t7s1_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          }
        }

        if(ch_14[2, i, now] == 3){
          if (next == 9)
              m_t7s3_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }

        if(ch_14[2, i, now] == 4){
          if (next == 9)
              m_t7s4_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }
      }

      if (now == 8) {
        if(ch_14[2, i, now] == 1){
          if (next == 9)
            m_t8s1_t9s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
          if (next == 11)
            m_t8s1_t11s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
        }
      }

      if (now == 9){
        if(ch_14[2, i, now] == 1)
          m_t9s1_t11s1_14[ch_14[1, i, now], ch_14[1, i, next]] += 1;
      }
    }

    {
      if (last == 0)
        l_0_s1_14 += 1;
      if (last == 1)
        l_1_s1_14[ch_14[1, i, last]] += 1;
      if (last == 4)
        l_4_s1_14[ch_14[1, i, last]] += 1;
      if (last == 5)
        l_5_s1_14[ch_14[1, i, last]] += 1;
      if (last == 6)
        l_6_s1_14[ch_14[1, i, last]] += 1;
      if (last == 7){
        if (ch_14[2, i, last] == 1)
          l_7_s1_14[ch_14[1, i, last]] += 1;
        if (ch_14[2, i, last] == 3)
          l_7_s3_14[ch_14[1, i, last]] += 1;
        if (ch_14[2, i, last] == 4)
          l_7_s4_14[ch_14[1, i, last]] += 1;
      }
      if (last == 8){
        if (ch_14[2, i, last] == 1)
          l_8_s1_14[ch_14[1, i, last]] += 1;
        if (ch_14[2, i, last] == 5)
          l_8_s5_14[ch_14[1, i, last]] += 1;
      }
      if (last == 9){
        if (ch_14[2, i, last] == 1)
          l_9_s1_14[ch_14[1, i, last]] += 1;
        if (ch_14[2, i, last] == 5)
          l_9_s5_14[ch_14[1, i, last]] += 1;
      }
    }
  }

  //Tranformed Data 2015----\\ 
  for (i in 1:N_15){
    int first = first_capture(ch_15[2, i, ]);
    int last = last_capture(ch_15[2, i, ]);
    
    if (rls_grp_15[i] == 1)
      N_rls_grp_15[1] += 1;
    else
      N_rls_grp_15[2] += 1;
    {
      if (first == 1)
        m_t0s1_t1s1_15[rls_grp_15[i], ch_15[1, i, first]] += 1;
      if (first == 4)
        m_t0s1_t4s1_15[rls_grp_15[i], ch_15[1, i, first]] += 1;
      if (first == 5)
        m_t0s1_t5s1_15[rls_grp_15[i], ch_15[1, i, first]] += 1;
      if (first == 7){
        if (ch_15[2, i, first] == 1)
          m_t0s1_t7s1_15[rls_grp_15[i], ch_15[1, i, first]] += 1;
        if (ch_15[2, i, first] == 3)
          m_t0s1_t7s3_15[rls_grp_15[i], ch_15[1, i, first]] += 1;
        if (ch_15[2, i, first] == 4)
          m_t0s1_t7s4_15[rls_grp_15[i], ch_15[1, i, first]] += 1;
      }
    }

    if (first > 0 && first != last){
      int now = first;
      int next = next_capture(ch_15[2, i, ], now, last);

      while (next < last){
        if (now == 1){
          if (next == 4)
            m_t1s1_t4s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 5)
            m_t1s1_t5s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 7){
            if (ch_15[2, i, next] == 1)
              m_t1s1_t7s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 3)
              m_t1s1_t7s3_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 4)
              m_t1s1_t7s4_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_15[2, i, next] == 1)
              m_t1s1_t8s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
        }

        if (now == 4) {
          if (next == 5)
            m_t4s1_t5s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 6)
            m_t4s1_t6s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 7){
            if (ch_15[2, i, next] == 1)
              m_t4s1_t7s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 3)
              m_t4s1_t7s3_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 4)
              m_t4s1_t7s4_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_15[2, i, next] == 1)
              m_t4s1_t8s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 5)
              m_t4s1_t8s5_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_15[2, i, next] == 1)
              m_t4s1_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }

        }

        if (now == 5) {
          if (next == 6)
            m_t5s1_t6s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 7){
            if (ch_15[2, i, next] == 1)
              m_t5s1_t7s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 3)
              m_t5s1_t7s3_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 4)
              m_t5s1_t7s4_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_15[2, i, next] == 1)
              m_t5s1_t8s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 5)
              m_t5s1_t8s5_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_15[2, i, next] == 1)
              m_t5s1_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
        }

        if (now == 6) {
          if (next == 7){
            if (ch_15[2, i, next] == 1)
              m_t6s1_t7s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 3)
              m_t6s1_t7s3_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 4)
              m_t6s1_t7s4_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
        }

        if (now == 7) {
          if(ch_15[2, i, now] == 1){
            if (next == 8){
              if (ch_15[2, i, next] == 1)
                m_t7s1_t8s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
              if (ch_15[2, i, next] == 5)
                m_t7s1_t8s5_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            }
            if (next == 9){
              if (ch_15[2, i, next] == 1)
                m_t7s1_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            }
          }

          if(ch_15[2, i, now] == 3){
            if (next == 9)
                m_t7s3_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }

          if(ch_15[2, i, now] == 4){
            if (next == 9)
                m_t7s4_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
        }

        if (now == 8) {
          if(ch_15[2, i, now] == 1){
            if (next == 9)
              m_t8s1_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
          if(ch_15[2, i, now] == 5){
            if (next == 9)
              m_t8s5_t9s5_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
        }

        now = next;
        next = next_capture(ch_15[2, i, ], now, last);
      }

      // from second-to-last capture to last capture

      if (now == 1){
        if (next == 4)
          m_t1s1_t4s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        if (next == 5)
          m_t1s1_t5s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        if (next == 7){
          if (ch_15[2, i, next] == 1)
            m_t1s1_t7s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 3)
            m_t1s1_t7s3_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 4)
            m_t1s1_t7s4_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
        if (next == 8){
          if (ch_15[2, i, next] == 1)
            m_t1s1_t8s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
      }

      if (now == 4) {
        if (next == 5)
          m_t4s1_t5s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        if (next == 6)
          m_t4s1_t6s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        if (next == 7){
          if (ch_15[2, i, next] == 1)
            m_t4s1_t7s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 3)
            m_t4s1_t7s3_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 4)
            m_t4s1_t7s4_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
        if (next == 8){
          if (ch_15[2, i, next] == 1)
            m_t4s1_t8s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 5)
            m_t4s1_t8s5_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
        if (next == 9){
          if (ch_15[2, i, next] == 1)
            m_t4s1_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
      }

      if (now == 5) {
        if (next == 6)
          m_t5s1_t6s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        if (next == 7){
          if (ch_15[2, i, next] == 1)
            m_t5s1_t7s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 3)
            m_t5s1_t7s3_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 4)
            m_t5s1_t7s4_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
        if (next == 8){
          if (ch_15[2, i, next] == 1)
            m_t5s1_t8s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 5)
            m_t5s1_t8s5_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
        if (next == 9){
          if (ch_15[2, i, next] == 1)
            m_t5s1_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
      }

      if (now == 6) {
        if (next == 7){
          if (ch_15[2, i, next] == 1)
            m_t6s1_t7s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 3)
            m_t6s1_t7s3_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (ch_15[2, i, next] == 4)
            m_t6s1_t7s4_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
      }

      if (now == 7) {
        if(ch_15[2, i, now] == 1){
          if (next == 8){
            if (ch_15[2, i, next] == 1)
              m_t7s1_t8s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
            if (ch_15[2, i, next] == 5)
              m_t7s1_t8s5_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_15[2, i, next] == 1)
              m_t7s1_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          }
        }

        if(ch_15[2, i, now] == 3){
          if (next == 9)
              m_t7s3_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 11)
            m_t7s3_t11s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }

        if(ch_15[2, i, now] == 4){
          if (next == 9)
              m_t7s4_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 11)
            m_t7s4_t11s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
      }

      if (now == 8) {
        if(ch_15[2, i, now] == 1){
          if (next == 9)
            m_t8s1_t9s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 11)
            m_t8s1_t11s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
        if(ch_15[2, i, now] == 5){
          if (next == 9)
            m_t8s5_t9s5_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
          if (next == 11)
            m_t8s5_t11s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        }
      }

      if (now == 9){
        if(ch_15[2, i, now] == 1)
          m_t9s1_t11s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
        if(ch_15[2, i, now] == 5)
          m_t9s5_t11s1_15[ch_15[1, i, now], ch_15[1, i, next]] += 1;
      }
    }

    {
      if (last == 0)
        l_0_s1_15[rls_grp_15[i]] += 1;
      if (last == 1)
        l_1_s1_15[ch_15[1, i, last]] += 1;
      if (last == 4)
        l_4_s1_15[ch_15[1, i, last]] += 1;
      if (last == 5)
        l_5_s1_15[ch_15[1, i, last]] += 1;
      if (last == 6)
        l_6_s1_15[ch_15[1, i, last]] += 1;
      if (last == 7){
        if (ch_15[2, i, last] == 1)
          l_7_s1_15[ch_15[1, i, last]] += 1;
        if (ch_15[2, i, last] == 3)
          l_7_s3_15[ch_15[1, i, last]] += 1;
        if (ch_15[2, i, last] == 4)
          l_7_s4_15[ch_15[1, i, last]] += 1;
      }
      if (last == 8){
        if (ch_15[2, i, last] == 1)
          l_8_s1_15[ch_15[1, i, last]] += 1;
        if (ch_15[2, i, last] == 5)
          l_8_s5_15[ch_15[1, i, last]] += 1;
      }
      if (last == 9){
        if (ch_15[2, i, last] == 1)
          l_9_s1_15[ch_15[1, i, last]] += 1;
        if (ch_15[2, i, last] == 5)
          l_9_s5_15[ch_15[1, i, last]] += 1;
      }
    }
  }

  //Tranformed Data 2016----\\ 
  for (i in 1:N_16){
    int first = first_capture(ch_16[2, i, ]);
    int last = last_capture(ch_16[2, i, ]);
    
    if (rls_grp_16[i] == 1)
      N_rls_grp_16[1] += 1;
    else
      N_rls_grp_16[2] += 1;
    
    {
      if (first == 1)
        m_t0s1_t1s1_16[rls_grp_16[i], ch_16[1, i, first]] += 1;
      if (first == 4){
        if (ch_16[2, i, first] == 1)
          m_t0s1_t4s1_16[rls_grp_16[i], ch_16[1, i, first]] += 1;
        if (ch_16[2, i, first] == 2)
          m_t0s1_t4s2_16[rls_grp_16[i], ch_16[1, i, first]] += 1;
      }
      if (first == 5)
        m_t0s1_t4s1_16[rls_grp_16[i], ch_16[1, i, first]] += 1;
      if (first == 7){
        if (ch_16[2, i, first] == 1)
          m_t0s1_t7s1_16[rls_grp_16[i], ch_16[1, i, first]] += 1;
      }
      if (first == 10)
        m_t0s1_t10s1_16[rls_grp_16[i], ch_16[1, i, first]] += 1;
    }

    if (first > 0 && first != last){
      int now = first;
      int next = next_capture(ch_16[2, i, ], now, last);

      while (next < last){
        if (now == 1){
          if (next == 4){
            if (ch_16[2, i, next] == 1)
            m_t1s1_t4s1_16[rls_grp_16[i], ch_16[1, i, next]] += 1;
            if (ch_16[2, i, next] == 2)
            m_t1s1_t4s2_16[rls_grp_16[i], ch_16[1, i, next]] += 1;
          }
          if (next == 5)
            m_t1s1_t5s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          if (next == 7){
            if (ch_16[2, i, next] == 4)
              m_t1s1_t7s4_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 4) {
          if(ch_16[2, i, now] == 1){
            if (next == 5)
              m_t4s1_t5s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 6)
              m_t4s1_t6s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 7){
              if (ch_16[2, i, next] == 1)
                m_t4s1_t7s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            }
          }
          
          if(ch_16[2, i, now] == 2){
            if (next == 9)
              m_t4s2_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 5) {
          if (next == 6)
            m_t5s1_t6s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          if (next == 7){
            if (ch_16[2, i, next] == 4)
              m_t5s1_t7s4_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_16[2, i, next] == 1)
              m_t5s1_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 6) {
          if (next == 7){
            if (ch_16[2, i, next] == 1)
              m_t6s1_t7s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (ch_16[2, i, next] == 3)
              m_t6s1_t7s3_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (ch_16[2, i, next] == 4)
              m_t6s1_t7s4_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_16[2, i, next] == 1)
              m_t6s1_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 7) {
          if(ch_16[2, i, now] == 1){
            if (next == 9){
              if (ch_16[2, i, next] == 1)
                m_t7s1_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
              if (ch_16[2, i, next] == 5)
                m_t7s1_t9s5_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            }
            if (next == 10)
              m_t7s1_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }

          if(ch_16[2, i, now] == 3){
            if (next == 9)
              m_t7s3_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 10)
              m_t7s3_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }

          if(ch_16[2, i, now] == 4){
            if (next == 9)
              m_t7s4_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 10)
              m_t7s4_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 9) {
          if(ch_16[2, i, now] == 1){
            if (next == 10)
              m_t9s1_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
          if(ch_16[2, i, now] == 5){
            if (next == 10)
              m_t9s5_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }
        now = next;
        next = next_capture(ch_16[2, i, ], now, last);
      }

      // from second-to-last capture to last capture

      {
        if (now == 1){
          if (next == 4){
            if (ch_16[2, i, next] == 1)
              m_t1s1_t4s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (ch_16[2, i, next] == 2)
              m_t1s1_t4s2_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
          if (next == 5)
            m_t1s1_t5s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          if (next == 7){
            if (ch_16[2, i, next] == 4)
              m_t1s1_t7s4_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 4) {
          if(ch_16[2, i, now] == 1){
            if (next == 5)
            m_t4s1_t5s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 6)
            m_t4s1_t6s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 7){
              if (ch_16[2, i, next] == 1)
                m_t4s1_t7s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            }
          }
          if(ch_16[2, i, now] == 2){
            if (next == 9)
            m_t4s2_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 5) {
          if (next == 6)
            m_t5s1_t6s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          if (next == 7){
            if (ch_16[2, i, next] == 4)
              m_t5s1_t7s4_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_16[2, i, next] == 1)
              m_t5s1_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 6) {
          if (next == 7){
            if (ch_16[2, i, next] == 1)
              m_t6s1_t7s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (ch_16[2, i, next] == 3)
              m_t6s1_t7s3_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (ch_16[2, i, next] == 4)
              m_t6s1_t7s4_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_16[2, i, next] == 1)
              m_t6s1_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 7) {
          if(ch_16[2, i, now] == 1){
            if (next == 9){
              if (ch_16[2, i, next] == 1)
                m_t7s1_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
              if (ch_16[2, i, next] == 5)
                m_t7s1_t9s5_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            }
            if (next == 10)
              m_t7s1_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }

          if(ch_16[2, i, now] == 3){
            if (next == 9)
              m_t7s3_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 10)
              m_t7s3_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }

          if(ch_16[2, i, now] == 4){
            if (next == 9)
              m_t7s4_t9s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 10)
              m_t7s4_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 9) {
          if(ch_16[2, i, now] == 1){
            if (next == 10)
              m_t9s1_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 11)
              m_t9s1_t11s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
          if(ch_16[2, i, now] == 5){
            if (next == 10)
              m_t9s5_t10s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
            if (next == 11)
              m_t9s5_t11s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;
          }
        }

        if (now == 10)
          if (next == 11)
            m_t10s1_t11s1_16[ch_16[1, i, now], ch_16[1, i, next]] += 1;

      }

    }

    {
      if (last == 0)
        l_0_s1_16[rls_grp_16[i]] += 1;
      if (last == 1)
        l_1_s1_16[ch_16[1, i, last]] += 1;
      if (last == 4) {
        if (ch_16[2, i, last] == 1)
          l_4_s1_16[ch_16[1, i, last]] += 1;
        if (ch_16[2, i, last] == 2)
          l_4_s2_16[ch_16[1, i, last]] += 1;
      }
      if (last == 5)
        l_5_s1_16[ch_16[1, i, last]] += 1;
      if (last == 6)
        l_6_s1_16[ch_16[1, i, last]] += 1;
      if (last == 7){
        if (ch_16[2, i, last] == 1)
          l_7_s1_16[ch_16[1, i, last]] += 1;
        if (ch_16[2, i, last] == 3)
          l_7_s3_16[ch_16[1, i, last]] += 1;
        if (ch_16[2, i, last] == 4)
          l_7_s4_16[ch_16[1, i, last]] += 1;
      }
      if (last == 9){
        if (ch_16[2, i, last] == 1)
          l_9_s1_16[ch_16[1, i, last]] += 1;
        if (ch_16[2, i, last] == 5)
          l_9_s5_16[ch_16[1, i, last]] += 1;
      }
      if (last == 10)
        l_10_s1_16[ch_16[1, i, last]] +=1;
    }
  }

  //Tranformed Data 2017----\\ 
  for (i in 1:N_17){
    int first = first_capture(ch_17[2, i, ]);
    int last = last_capture(ch_17[2, i, ]);
    {
      if (first == 1)
        m_t0s1_t1s1_17[ch_17[1, i, first]] += 1;
      if (first == 3)
          m_t0s1_t3s1_17[ch_17[1, i, first]] += 1;
      if (first == 4){
        if (ch_17[2, i, first] == 1)
          m_t0s1_t4s1_17[ch_17[1, i, first]] += 1;
      }
      if (first == 5)
        m_t0s1_t4s1_17[ch_17[1, i, first]] += 1;
      if (first == 6)
        m_t0s1_t4s1_17[ch_17[1, i, first]] += 1;
      if (first == 7){
        if (ch_17[2, i, first] == 1)
          m_t0s1_t7s1_17[ch_17[1, i, first]] += 1;
      }
      if (first == 8){
        if (ch_17[2, i, first] == 1)
          m_t0s1_t8s1_17[ch_17[1, i, first]] += 1;
      }
      if (first == 9){
        if (ch_17[2, i, first] == 1)
          m_t0s1_t9s1_17[ch_17[1, i, first]] += 1;
      }
      if (first == 10)
        m_t0s1_t10s1_17[ch_17[1, i, first]] += 1;
    }

    if (first > 0 && first != last){
      int now = first;
      int next = next_capture(ch_17[2, i, ], now, last);

      while (next < last){
        if (now == 1){
          if (next == 2)
            m_t1s1_t2s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 3)
            m_t1s1_t3s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 4){
            if (ch_17[2, i, next] == 2)
            m_t1s1_t4s2_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t1s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
        }

        if (now == 2) {
          if (next == 3)
            m_t2s1_t3s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 4){
            if (ch_17[2, i, next] == 1)
            m_t2s1_t4s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 2)
            m_t2s1_t4s2_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 5)
            m_t2s1_t5s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t2s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 10)
            m_t2s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 3) {
          if (next == 4)
            m_t3s1_t4s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 5)
            m_t3s1_t5s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 6)
            m_t3s1_t6s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 7){
            if (ch_17[2, i, next] == 1)
              m_t3s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 3)
              m_t3s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 4)
              m_t3s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_17[2, i, next] == 1)
              m_t3s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t3s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 5)
              m_t3s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 10)
            m_t3s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 4) {
          if(ch_17[2, i, now] == 1){
            if (next == 5)
              m_t4s1_t5s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 6)
              m_t4s1_t6s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 7){
              if (ch_17[2, i, next] == 1)
                m_t4s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            }
          }
          if(ch_17[2, i, now] == 2){
            if (next == 9)
              m_t4s2_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 10)
              m_t4s2_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
        }

        if (now == 5) {
          if (next == 6)
            m_t5s1_t6s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 7){
            if (ch_17[2, i, next] == 1)
              m_t5s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 3)
              m_t5s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 4)
              m_t5s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_17[2, i, next] == 1)
              m_t5s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t5s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 5)
              m_t5s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 10)
            m_t5s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 6) {
          if (next == 7){
            if (ch_17[2, i, next] == 1)
              m_t6s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 3)
              m_t6s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 4)
              m_t6s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_17[2, i, next] == 1)
              m_t6s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t6s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 10)
            m_t6s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 7) {
          if(ch_17[2, i, now] == 1){
            if (next == 8){
              if (ch_17[2, i, next] == 1)
                m_t7s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            }
            if (next == 9){
              if (ch_17[2, i, next] == 1)
                m_t7s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
              if (ch_17[2, i, next] == 5)
                m_t7s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            }
            if (next == 10)
              m_t7s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }

          if(ch_17[2, i, now] == 3){
            if (next == 9)
              m_t7s3_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 10)
              m_t7s3_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }

          if(ch_17[2, i, now] == 4){
            if (next == 9)
              m_t7s4_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 10)
              m_t7s4_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
        }

        if (now == 8) {
          if(ch_17[2, i, now] == 1){
            if (next == 9)
              m_t8s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 10)
              m_t8s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
        }

        if (now == 9) {
          if(ch_17[2, i, now] == 1){
            if (next == 10)
              m_t9s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if(ch_17[2, i, now] == 5){
            if (next == 10)
              m_t9s5_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
        }

        now = next;
        next = next_capture(ch_17[2, i, ], now, last);
      }

      // from second-to-last capture to last capture

    {
        if (now == 1){
          if (next == 2)
            m_t1s1_t2s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 3)
            m_t1s1_t3s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 4){
            // if (ch_17[2, i, next] == 1)
            //   m_t1s1_t4s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 2)
              m_t1s1_t4s2_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          // if (next == 5)
          //   m_t1s1_t5s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 6)
          //   m_t1s1_t6s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 7){
          //   if (ch_17[2, i, next] == 1)
          //     m_t1s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          //   if (ch_17[2, i, next] == 3)
          //     m_t1s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          //   if (ch_17[2, i, next] == 4)
          //     m_t1s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // }
          // if (next == 8){
          //   if (ch_17[2, i, next] == 1)
          //     m_t1s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          //   if (ch_17[2, i, next] == 5)
          //     m_t1s1_t8s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t1s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (ch_17[2, i, next] == 5)
            //   m_t1s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          // if (next == 10)
          //   m_t1s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 11)
          //   m_t1s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 2) {
          if (next == 3)
            m_t2s1_t3s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 4){
            if (ch_17[2, i, next] == 1)
              m_t2s1_t4s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 2)
              m_t2s1_t4s2_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 5)
            m_t2s1_t5s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 6)
          //   m_t2s1_t6s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 7){
          //   if (ch_17[2, i, next] == 1)
          //     m_t2s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          //   if (ch_17[2, i, next] == 3)
          //     m_t2s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          //   if (ch_17[2, i, next] == 4)
          //     m_t2s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // }
          // if (next == 8){
          //   if (ch_17[2, i, next] == 1)
          //     m_t2s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          //   if (ch_17[2, i, next] == 5)
          //     m_t2s1_t8s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t2s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (ch_17[2, i, next] == 5)
            //   m_t2s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 10)
            m_t2s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 11)
          //   m_t2s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 3) {
          if (next == 4)
            m_t3s1_t4s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 5)
            m_t3s1_t5s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 6)
            m_t3s1_t6s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 7){
            if (ch_17[2, i, next] == 1)
              m_t3s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 3)
              m_t3s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 4)
              m_t3s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_17[2, i, next] == 1)
              m_t3s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (ch_17[2, i, next] == 5)
            //   m_t3s1_t8s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t3s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 5)
              m_t3s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 10)
            m_t3s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 11)
          //   m_t3s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 4) {
          if(ch_17[2, i, now] == 1){
            if (next == 5)
              m_t4s1_t5s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 6)
              m_t4s1_t6s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 7){
              if (ch_17[2, i, next] == 1)
                m_t4s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
              // if (ch_17[2, i, next] == 3)
              //   m_t4s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
              // if (ch_17[2, i, next] == 4)
              //   m_t4s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            }
            // if (next == 8){
            //   if (ch_17[2, i, next] == 1)
            //     m_t4s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            //   if (ch_17[2, i, next] == 5)
            //     m_t4s1_t8s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // }
            // if (next == 9){
            //   if (ch_17[2, i, next] == 1)
            //     m_t4s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            //   if (ch_17[2, i, next] == 5)
            //     m_t4s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // }
            // if (next == 10)
            //   m_t4s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (next == 11)
            //   m_t4s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if(ch_17[2, i, now] == 2){
            if (next == 9)
              m_t4s2_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 10)
              m_t4s2_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 11)
              m_t4s2_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
        }

        if (now == 5) {
          if (next == 6)
            m_t5s1_t6s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          if (next == 7){
            if (ch_17[2, i, next] == 1)
              m_t5s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 3)
              m_t5s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 4)
              m_t5s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_17[2, i, next] == 1)
              m_t5s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (ch_17[2, i, next] == 5)
            //   m_t5s1_t8s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t5s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 5)
              m_t5s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 10)
            m_t5s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 11)
          //   m_t5s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 6) {
          if (next == 7){
            if (ch_17[2, i, next] == 1)
              m_t6s1_t7s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 3)
              m_t6s1_t7s3_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (ch_17[2, i, next] == 4)
              m_t6s1_t7s4_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_17[2, i, next] == 1)
              m_t6s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (ch_17[2, i, next] == 5)
            //   m_t6s1_t8s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_17[2, i, next] == 1)
              m_t6s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (ch_17[2, i, next] == 5)
            //   m_t6s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if (next == 10)
            m_t6s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // if (next == 11)
          //   m_t6s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
        }

        if (now == 7) {
          if(ch_17[2, i, now] == 1){
            if (next == 8){
              if (ch_17[2, i, next] == 1)
                m_t7s1_t8s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
              // if (ch_17[2, i, next] == 5)
              //   m_t7s1_t8s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            }
            if (next == 9){
              if (ch_17[2, i, next] == 1)
                m_t7s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
              if (ch_17[2, i, next] == 5)
                m_t7s1_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            }
            if (next == 10)
              m_t7s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (next == 11)
            //   m_t7s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }

          if(ch_17[2, i, now] == 3){
            if (next == 9)
              m_t7s3_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 10)
              m_t7s3_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 11)
              m_t7s3_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }

          if(ch_17[2, i, now] == 4){
            if (next == 9)
              m_t7s4_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 10)
              m_t7s4_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 11)
              m_t7s4_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
        }

        if (now == 8) {
          if(ch_17[2, i, now] == 1){
            if (next == 9)
              m_t8s1_t9s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 10)
              m_t8s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (next == 11)
            //   m_t8s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          // if(ch_17[2, i, now] == 5){
          //   if (next == 9)
          //     m_t8s5_t9s5_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          //   if (next == 10)
          //     m_t8s5_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          //   if (next == 11)
          //     m_t8s5_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          // }
        }

        if (now == 9) {
          if(ch_17[2, i, now] == 1){
            if (next == 10)
              m_t9s1_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            if (next == 11)
              m_t9s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
          if(ch_17[2, i, now] == 5){
            if (next == 10)
              m_t9s5_t10s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
            // if (next == 11)
            //   m_t9s5_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;
          }
        }
      if (now == 10)
         if (next == 11)
            m_t10s1_t11s1_17[ch_17[1, i, now], ch_17[1, i, next]] += 1;

    }

    }

    {
      if (last == 0)
        l_0_s1_17 += 1;
      if (last == 1)
        l_1_s1_17[ch_17[1, i, last]] += 1;
      if (last == 2)
        l_2_s1_17[ch_17[1, i, last]] += 1;
      if (last == 3)
        l_3_s1_17[ch_17[1, i, last]] += 1;
      if (last == 4) {
        if (ch_17[2, i, last] == 1)
          l_4_s1_17[ch_17[1, i, last]] += 1;
        if (ch_17[2, i, last] == 2)
          l_4_s2_17[ch_17[1, i, last]] += 1;
      }
      if (last == 5)
        l_5_s1_17[ch_17[1, i, last]] += 1;
      if (last == 6)
        l_6_s1_17[ch_17[1, i, last]] += 1;

      if (last == 7){
        if (ch_17[2, i, last] == 1)
          l_7_s1_17[ch_17[1, i, last]] += 1;
        if (ch_17[2, i, last] == 3)
          l_7_s3_17[ch_17[1, i, last]] += 1;
        if (ch_17[2, i, last] == 4)
          l_7_s4_17[ch_17[1, i, last]] += 1;
      }
      if (last == 8){
        if (ch_17[2, i, last] == 1)
          l_8_s1_17[ch_17[1, i, last]] += 1;
        if (ch_17[2, i, last] == 5)
          l_8_s5_17[ch_17[1, i, last]] += 1;
      }
      if (last == 9){
        if (ch_17[2, i, last] == 1)
          l_9_s1_17[ch_17[1, i, last]] += 1;
        if (ch_17[2, i, last] == 5)
          l_9_s5_17[ch_17[1, i, last]] += 1;
      }
      if (last == 10)
        l_10_s1_17[ch_17[1, i, last]] +=1;
    }
  }

  //Tranformed Data 2018----\\ 

  for (g in 1:2){
    l_10_s1_18[g] = rep_vector(0, T);
    l_9_s1_18[g] = rep_vector(0, T); 
    l_9_s5_18[g] = rep_vector(0, T); 
    l_8_s1_18[g] = rep_vector(0, T);
    l_8_s5_18[g] = rep_vector(0, T); 
    l_7_s1_18[g] = rep_vector(0, T); 
    l_7_s3_18[g] = rep_vector(0, T); 
    l_7_s4_18[g] = rep_vector(0, T); 
    l_6_s1_18[g] = rep_vector(0, T); 
    l_5_s1_18[g] = rep_vector(0, T); 
    l_4_s1_18[g] = rep_vector(0, T); 
    l_3_s1_18[g] = rep_vector(0, T); 
    l_2_s1_18[g] = rep_vector(0, T); 
    l_1_s1_18[g] = rep_vector(0, T); 
    l_0_s1_18[g] = 0;
    
    m_t0s1_t1s1_18[g] = rep_row_vector(0, T);

    m_t1s1_t2s1_18[g] = rep_matrix(0, T, T);
    m_t1s1_t3s1_18[g] = rep_matrix(0, T, T); 
    m_t1s1_t4s1_18[g] = rep_matrix(0, T, T); 
    m_t1s1_t5s1_18[g] = rep_matrix(0, T, T); 
    m_t1s1_t6s1_18[g] = rep_matrix(0, T, T); 
    m_t1s1_t7s1_18[g] = rep_matrix(0, T, T); 

    m_t2s1_t3s1_18[g] = rep_matrix(0, T, T); 
    m_t2s1_t4s1_18[g] = rep_matrix(0, T, T); 
    m_t2s1_t5s1_18[g] = rep_matrix(0, T, T); 
    m_t2s1_t6s1_18[g] = rep_matrix(0, T, T); 
    m_t2s1_t7s1_18[g] = rep_matrix(0, T, T); 
    m_t2s1_t7s3_18[g] = rep_matrix(0, T, T); 
    m_t2s1_t7s4_18[g] = rep_matrix(0, T, T); 
    m_t2s1_t9s1_18[g] = rep_matrix(0, T, T); 

    m_t3s1_t4s1_18[g] = rep_matrix(0, T, T); 
    m_t3s1_t5s1_18[g] = rep_matrix(0, T, T); 
    m_t3s1_t6s1_18[g] = rep_matrix(0, T, T); 
    m_t3s1_t7s1_18[g] = rep_matrix(0, T, T); 
    m_t3s1_t9s1_18[g] = rep_matrix(0, T, T); 

    m_t4s1_t5s1_18[g] = rep_matrix(0, T, T); 
    m_t4s1_t6s1_18[g] = rep_matrix(0, T, T); 
    m_t4s1_t7s1_18[g] = rep_matrix(0, T, T);
    m_t4s1_t7s3_18[g] = rep_matrix(0, T, T); 
    m_t4s1_t7s4_18[g] = rep_matrix(0, T, T); 

    m_t5s1_t6s1_18[g] = rep_matrix(0, T, T); 
    m_t5s1_t7s1_18[g] = rep_matrix(0, T, T); 
    m_t5s1_t7s3_18[g] = rep_matrix(0, T, T); 
    m_t5s1_t7s4_18[g] = rep_matrix(0, T, T); 
    
    m_t6s1_t7s1_18[g] = rep_matrix(0, T, T); 
    m_t6s1_t7s3_18[g] = rep_matrix(0, T, T); 
    m_t6s1_t7s4_18[g] = rep_matrix(0, T, T); 
    m_t6s1_t8s5_18[g] = rep_matrix(0, T, T);
    
    m_t7s1_t8s1_18[g] = rep_matrix(0, T, T); 
    m_t7s1_t8s5_18[g] = rep_matrix(0, T, T); 
    m_t7s1_t9s1_18[g] = rep_matrix(0, T, T); 
    m_t7s1_t10s1_18[g] = rep_matrix(0, T, T); 

    m_t7s3_t9s1_18[g] = rep_matrix(0, T, T);
    m_t7s3_t10s1_18[g] = rep_matrix(0, T, T); 

    m_t7s4_t9s1_18[g] = rep_matrix(0, T, T);
    m_t7s4_t10s1_18[g] = rep_matrix(0, T, T); 

    m_t8s1_t9s1_18[g] = rep_matrix(0, T, T); 
    m_t8s1_t10s1_18[g] = rep_matrix(0, T, T); 

    m_t8s5_t9s5_18[g] = rep_matrix(0, T, T); 

    m_t9s1_t10s1_18[g] = rep_matrix(0, T, T); 
    m_t9s1_t11s1_18[g] = rep_matrix(0, T, T); 

    m_t9s5_t10s1_18[g] = rep_matrix(0, T, T); 

    m_t10s1_t11s1_18[g] = rep_matrix(0, T, T); 
  }
  
  for (i in 1:N_18){
    int r = rls_grp_18[i];
    int first = first_capture(ch_18[2, i, ]);
    int last = last_capture(ch_18[2, i, ]);
    
    if (rls_grp_18[i] == 1)
      N_rls_grp_18[1] += 1;
    else
      N_rls_grp_18[2] += 1;
    
    {
      if (first == 1)
        m_t0s1_t1s1_18[r, ch_18[1, i, first]] += 1;
    }

    if (first > 0 && first != last){
      int now = first;
      int next = next_capture(ch_18[2, i, ], now, last);

      while (next < last){
        if (now == 1){
          if (next == 2)
            m_t1s1_t2s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 3)
            m_t1s1_t3s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 4)
            m_t1s1_t4s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 5)
            m_t1s1_t5s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 6)
            m_t1s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t1s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 2) {
          if (next == 3)
            m_t2s1_t3s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 4)
            m_t2s1_t4s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 5)
            m_t2s1_t5s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 6)
            m_t2s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t2s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 3)
              m_t2s1_t7s3_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 4)
              m_t2s1_t7s4_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_18[2, i, next] == 1)
              m_t2s1_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 3) {
          if (next == 4)
            m_t3s1_t4s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 5)
            m_t3s1_t5s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 6)
            m_t3s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t3s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_18[2, i, next] == 1)
              m_t3s1_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 4) {
          if (next == 5)
            m_t4s1_t5s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 6)
            m_t4s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t4s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 3)
              m_t4s1_t7s3_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 4)
              m_t4s1_t7s4_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 5) {
          if (next == 6)
            m_t5s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t5s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 3)
              m_t5s1_t7s3_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 4)
              m_t5s1_t7s4_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 6) {
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t6s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 3)
              m_t6s1_t7s3_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 4)
              m_t6s1_t7s4_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
          if (next == 8){
            if (ch_18[2, i, next] == 5)
              m_t6s1_t8s5_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 7) {
          if(ch_18[2, i, now] == 1){
            if (next == 8){
              if (ch_18[2, i, next] == 1)
                m_t7s1_t8s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            }
            if (next == 9){
              if (ch_18[2, i, next] == 1)
                m_t7s1_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            }
            if (next == 10)
              m_t7s1_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }

          if(ch_18[2, i, now] == 3){
            if (next == 9)
              m_t7s3_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (next == 10)
              m_t7s3_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }

          if(ch_18[2, i, now] == 4){
            if (next == 9)
              m_t7s4_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (next == 10)
              m_t7s4_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 8) {
          if(ch_18[2, i, now] == 1){
            if (next == 9)
              m_t8s1_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (next == 10)
              m_t8s1_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
          if(ch_18[2, i, now] == 5){
            if (next == 9)
              m_t8s5_t9s5_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 9) {
          if(ch_18[2, i, now] == 1){
            if (next == 10)
              m_t9s1_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
          if(ch_18[2, i, now] == 5){
            if (next == 10)
              m_t9s5_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        now = next;
        next = next_capture(ch_18[2, i, ], now, last);
      }

      // from second-to-last capture to last capture

    {
        if (now == 1){
          if (next == 2)
            m_t1s1_t2s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 3)
            m_t1s1_t3s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (ch_18[2, i, next] == 1)
            m_t1s1_t4s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 5)
            m_t1s1_t5s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 6)
            m_t1s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t1s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 2) {
          if (next == 3)
            m_t2s1_t3s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (ch_18[2, i, next] == 1)
            m_t2s1_t4s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 5)
            m_t2s1_t5s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 6)
            m_t2s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t2s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 3)
              m_t2s1_t7s3_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 4)
              m_t2s1_t7s4_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }

          if (next == 9){
            if (ch_18[2, i, next] == 1)
              m_t2s1_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 3) {
          if (next == 4)
            m_t3s1_t4s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 5)
            m_t3s1_t5s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 6)
            m_t3s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t3s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
          if (next == 9){
            if (ch_18[2, i, next] == 1)
              m_t3s1_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 4) {
          if (next == 5)
            m_t4s1_t5s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 6)
            m_t4s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t4s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 3)
              m_t4s1_t7s3_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 4)
              m_t4s1_t7s4_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 5) {
          if (next == 6)
            m_t5s1_t6s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t5s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 3)
              m_t5s1_t7s3_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 4)
              m_t5s1_t7s4_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 6) {
          if (next == 7){
            if (ch_18[2, i, next] == 1)
              m_t6s1_t7s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 3)
              m_t6s1_t7s3_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (ch_18[2, i, next] == 4)
              m_t6s1_t7s4_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
          if (next == 8){

            if (ch_18[2, i, next] == 5)
              m_t6s1_t8s5_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }

        }

        if (now == 7) {
          if(ch_18[2, i, now] == 1){
            if (next == 8){
              if (ch_18[2, i, next] == 1)
                m_t7s1_t8s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
              if (ch_18[2, i, next] == 5)
                m_t7s1_t8s5_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            }
            if (next == 9){
              if (ch_18[2, i, next] == 1)
                m_t7s1_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;

            }
            if (next == 10)
              m_t7s1_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;

          }

          if(ch_18[2, i, now] == 3){
            if (next == 9)
              m_t7s3_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (next == 10)
              m_t7s3_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }

          if(ch_18[2, i, now] == 4){
            if (next == 9)
              m_t7s4_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (next == 10)
              m_t7s4_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;

          }
        }

        if (now == 8) {
          if(ch_18[2, i, now] == 1){
            if (next == 9)
              m_t8s1_t9s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (next == 10)
              m_t8s1_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;

          }
          if(ch_18[2, i, now] == 5){
            if (next == 9)
              m_t8s5_t9s5_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }

        if (now == 9) {
          if(ch_18[2, i, now] == 1){
            if (next == 10)
              m_t9s1_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
            if (next == 11)
              m_t9s1_t11s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
          if(ch_18[2, i, now] == 5){
            if (next == 10)
              m_t9s5_t10s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;
          }
        }
      if (now == 10)
         if (next == 11)
            m_t10s1_t11s1_18[r, ch_18[1, i, now], ch_18[1, i, next]] += 1;

    }

    }

    {
      if (last == 0)
        l_0_s1_18[r] += 1;
      if (last == 1)
        l_1_s1_18[r, ch_18[1, i, last]] += 1;
      if (last == 2)
        l_2_s1_18[r, ch_18[1, i, last]] += 1;
      if (last == 3)
        l_3_s1_18[r, ch_18[1, i, last]] += 1;
      if (last == 4) {
        if (ch_18[2, i, last] == 1)
          l_4_s1_18[r, ch_18[1, i, last]] += 1;
        if (ch_18[2, i, last] == 2)
          l_4_s2_18[r, ch_18[1, i, last]] += 1;
      }
      if (last == 5)
        l_5_s1_18[r, ch_18[1, i, last]] += 1;
      if (last == 6)
        l_6_s1_18[r, ch_18[1, i, last]] += 1;

      if (last == 7){
        if (ch_18[2, i, last] == 1)
          l_7_s1_18[r, ch_18[1, i, last]] += 1;
        if (ch_18[2, i, last] == 3)
          l_7_s3_18[r, ch_18[1, i, last]] += 1;
        if (ch_18[2, i, last] == 4)
          l_7_s4_18[r, ch_18[1, i, last]] += 1;
      }
      if (last == 8){
        if (ch_18[2, i, last] == 1)
          l_8_s1_18[r, ch_18[1, i, last]] += 1;
        if (ch_18[2, i, last] == 5)
          l_8_s5_18[r, ch_18[1, i, last]] += 1;
      }
      if (last == 9){
        if (ch_18[2, i, last] == 1)
          l_9_s1_18[r, ch_18[1, i, last]] += 1;
        if (ch_18[2, i, last] == 5)
          l_9_s5_18[r, ch_18[1, i, last]] += 1;
      }
      if (last == 10)
        l_10_s1_18[r, ch_18[1, i, last]] +=1;
    }
  }  
  
}

parameters {
  real<lower = 0, upper = 1> chipps2benecia; //survival chipps -> benecia

  real            beta_phi_s1_t0_14;    //survival release to knight's landing
  real            beta_phi_s1_t0_15;    //survival release to knight's landing
  real            beta_phi_s1_t0_16;    //survival release to knight's landing
  real            beta_phi_s1_t0_17;    //survival release to knight's landing
  real            beta_phi_s1_t0_18[2];    //survival release to knight's landing

  real            beta_phi_s1_t1[3]; //survival knight's landing -> abv fremont weir
  //assume survival to blw fremont is 1
  
  real            beta_phi_s1_t3[3]; //survival blw fremont -> feather
  real            beta_phi_s1_t4[3]; //survival feather -> sacramento
  real            beta_phi_s1_t5[3]; //survival sacramento -> freeport
  real            beta_phi_s1_t6[3]; //survival freeport -> sutter/steam complex
  real            beta_phi_s1_t7[3]; //survival sutter/steam complex -> georgiana
  real            beta_phi_s1_t8[3]; //survival georgiana -> rio vista
  real            beta_phi_s1_t9[3]; //survival rio vista -> chipps

  real            beta_phi_s2_t3[3]; //survival fremont weir -> toe drain
  real            beta_phi_s2_t8[3]; //survival yolo -> rio vista
  
  real            beta_phi_s3_t8[3]; //survival sutter -> rio vista
  real            beta_phi_s4_t8[3]; //survival steam -> rio vista
  
  real            beta_phi_s5_t8[3]; //survival georgian -> interior delta; intercept, georgian flow, temp
  real            beta_phi_s5_t9[4]; //survival interior delta -> chipps; intercept, i/e ratio, temp

  real            mu_v_s1_t0_14;    //velocity mean release to knight's landing
  real            mu_v_s1_t0_15;    //velocity mean release to knight's landing
  real            mu_v_s1_t0_16;    //velocity mean release to knight's landing
  real            mu_v_s1_t0_17;    //velocity mean release to knight's landing
  real            mu_v_s1_t0_18[2];    //velocity mean release to knight's landing

  real            beta_mu_s1_t1[2]; //logmean travel knight's landing -> abv fremont weir
  real            beta_mu_s1_t2b[2]; //logmean travel abv fremont weir -> feather (2014, 2015, and 2016)
  real            beta_mu_s1_t2[2]; //logmean travel abv fremont -> blw fremont
  real            beta_mu_s1_t3[2]; //logmean travel blw fremont -> feather
  real            beta_mu_s1_t4[2]; //logmean travel feather -> sacramento
  real            beta_mu_s1_t5[2]; //logmean travel sacramento -> freeport
  real            beta_mu_s1_t6[2]; //logmean travel freeport -> sutter/steam complex
  real            beta_mu_s1_t7[2]; //logmean travel sutter/steam complex -> georgiana
  real            beta_mu_s1_t8[2]; //logmean travel georgiana -> rio vista
  real            beta_mu_s1_t9[2]; //logmean travel rio vista -> chipps
  real            beta_mu_s1_t10; //logmean travel chipps -> benecia

  real            beta_mu_s2_t3[2]; //logmean travel fremont weir -> toe drain
  real            beta_mu_s2_t8[2]; //logmean travel yolo -> rio vista
  
  real            beta_mu_s3_t8[2]; //logmean travel sutter -> rio vista
  real            beta_mu_s4_t8[2]; //logmean travel steam -> rio vista
  
  real            beta_mu_s5_t8[2]; //logmean travel georgian -> interior delta 
  real            beta_mu_s5_t9[2]; //logmean travel interior delta -> chipps

  real<lower = 0> sigma_v_s1_t0_14;   //velocity variance release to knight's landing
  real<lower = 0> sigma_v_s1_t0_15;   //velocity variance release to knight's landing
  real<lower = 0> sigma_v_s1_t0_16;   //velocity variance release to knight's landing
  real<lower = 0> sigma_v_s1_t0_17;   //velocity variance release to knight's landing
  real<lower = 0> sigma_v_s1_t0_18;   //velocity variance release to knight's landing
  
  vector[T]       tt_eps_s1_t0_14;
  vector[T]       tt_eps_s1_t0_15;
  vector[T]       tt_eps_s1_t0_16;
  vector[T]       tt_eps_s1_t0_17;
  vector[T]       tt_eps_s1_t0_18;
  real<lower = 0> sigma_eps_s1_t0;  
  
  real<lower = 0> sigma_tt_s1_t1;   //lodsd travel knight's landing -> abv fremont weir
  real<lower = 0> sigma_tt_s1_t2b;   //lodsd travel knight's landing -> abv fremont weir
  real<lower = 0> sigma_tt_s1_t2;   //lodsd travel abv fremont -> blw fremont
  real<lower = 0> sigma_tt_s1_t3;   //lodsd travel blw fremont -> feather
  real<lower = 0> sigma_tt_s1_t4;   //lodsd travel feather -> sacramento
  real<lower = 0> sigma_tt_s1_t5;   //lodsd travel sacramento -> freeport
  real<lower = 0> sigma_tt_s1_t6;   //lodsd travel freeport -> sutter/steam complex
  real<lower = 0> sigma_tt_s1_t7;   //lodsd travel sutter/steam complex -> georgiana
  real<lower = 0> sigma_tt_s1_t8;   //lodsd travel georgiana -> rio vista
  real<lower = 0> sigma_tt_s1_t9;   //lodsd travel rio vista -> chipps
  real<lower = 0> sigma_tt_s1_t10;   //lodsd travel rio vista -> chipps
  real<lower = 0> sigma_tt_s2_t3;   //lodsd travel fremont weir -> toe drain
  real<lower = 0> sigma_tt_s2_t8;   //lodsd travel yolo -> rio vista
  real<lower = 0> sigma_tt_s3_t8;   //lodsd travel sutter -> rio vista
  real<lower = 0> sigma_tt_s4_t8;   //lodsd travel steam -> rio vista
  real<lower = 0> sigma_tt_s5_t8;   //lodsd travel georgian -> interior delta 
  real<lower = 0> sigma_tt_s5_t9;   //lodsd travel interior delta -> chipps

  real            beta_psi_1to2_t2[2]; // probability of entering yolo from sac
  real            beta_psi_1to3_t6[2]; // probability of entering sutter from sac
  real            beta_psi_1to4_t6[2]; // probability of entering steam from sac
  real            beta_psi_1to5_t7[2]; // probability of entering geo from sac
  
  real            beta_p_s1_t1_14[2];  //detection knight's landing
  real            beta_p_s1_t4_14[2];  //detection sac below feather
  real            beta_p_s1_t5_14[2];  //detection sacramento (tower or I80)
  real            beta_p_s1_t6_14[2];  //detection freeport
  real            beta_p_s1_t7_14[2];  //detection sac blw sutter/steam
  real            beta_p_s1_t8_14[2];  //detection sac blw geogiana
  real            beta_p_s1_t9_14[2];  //detection rio vista

  real            beta_p_s3_t7_14[2]; //detection sutter/steam entrance
  real            beta_p_s4_t7_14[2]; //detection sutter/steam entrance

  real            beta_p_s5_t8_14[2]; //detection georgiana slough entrance
  real            beta_p_s5_t9_14[2]; //detection interior delta

  real            beta_p_s1_t1_15[2];  //detection knight's landing
  real            beta_p_s1_t4_15[2];  //detection sac below feather
  real            beta_p_s1_t5_15[2];  //detection sacramento (tower or I80)
  real            beta_p_s1_t6_15[2];  //detection freeport
  real            beta_p_s1_t7_15[2];  //detection sac blw sutter/steam
  real            beta_p_s1_t8_15[2];  //detection sac blw geogiana
  real            beta_p_s1_t9_15[2];  //detection rio vista

  real            beta_p_s3_t7_15[2]; //detection sutter/steam entrance
  real            beta_p_s4_t7_15[2]; //detection sutter/steam entrance

  real            beta_p_s5_t8_15[2]; //detection georgiana slough entrance
  real            beta_p_s5_t9_15[2]; //detection interior delta

  real            beta_p_s1_t1_16[2];  //detection knight's landing
  real            beta_p_s1_t4_16[2];  //detection sac below feather
  real            beta_p_s1_t5_16[2];  //detection sacramento (tower or I80)
  real            beta_p_s1_t6_16[2];  //detection freeport
  real            beta_p_s1_t7_16[2];  //detection sac blw sutter/steam
  real            beta_p_s1_t9_16[2];  //detection rio vista
  real            beta_p_s1_t10_16[2]; //detection Chipps

  real            beta_p_s2_t4_16[2]; //detection toe drain

  real            beta_p_s3_t7_16[2]; //detection sutter/steam entrance
  real            beta_p_s4_t7_16[2]; //detection sutter/steam entrance

  real            beta_p_s5_t9_16[2]; //detection interior delta

  real            beta_p_s1_t1_17[2];  //detection knight's landing
  real            beta_p_s1_t2_17[2];  //detection abv fremont weir
  real            beta_p_s1_t3_17[2];  //detection blw fremont weir
  real            beta_p_s1_t4_17[2];  //detection sac below feather
  real            beta_p_s1_t5_17[2];  //detection sacramento (tower or I80)
  real            beta_p_s1_t6_17[2];  //detection freeport
  real            beta_p_s1_t7_17[2];  //detection sac blw sutter/steam
  real            beta_p_s1_t8_17[2];  //detection sac blw geogiana
  real            beta_p_s1_t9_17[2];  //detection rio vista
  real            beta_p_s1_t10_17[2]; //detection Chipps

  real            beta_p_s2_t4_17[2]; //detection toe drain

  real            beta_p_s3_t7_17[2]; //detection sutter/steam entrance
  real            beta_p_s4_t7_17[2]; //detection sutter/steam entrance

  real            beta_p_s5_t9_17[2]; //detection interior delta
  
  real            beta_p_s1_t1_18[2];  //detection knight's landing
  real            beta_p_s1_t2_18[2];  //detection abv fremont weir
  real            beta_p_s1_t3_18[2];  //detection blw fremont weir
  real            beta_p_s1_t4_18[2];  //detection sac below feather
  real            beta_p_s1_t5_18[2];  //detection sacramento (tower or I80)
  real            beta_p_s1_t6_18[2];  //detection freeport
  real            beta_p_s1_t7_18[2];  //detection sac blw sutter/steam
  real            beta_p_s1_t8_18[2];  //detection sac blw geogiana
  real            beta_p_s1_t9_18[2];  //detection rio vista
  real            beta_p_s1_t10_18[2]; //detection Chipps
  
  real            beta_p_s2_t4_18[2]; //detection toe drain

  real            beta_p_s3_t7_18[2]; //detection sutter/steam entrance
  real            beta_p_s4_t7_18[2]; //detection sutter/steam entrance

  real            beta_p_s5_t8_18[2]; //detection georgiana slough entrance
  real            beta_p_s5_t9_18[2]; //detection interior delta  
}

transformed parameters{
  //Tran Param 2014-----------------------------------------------------------------\\  
  real                 phi_s1_t0_14; //survival release to knight's landing
  vector[T]            phi_s1_t1_14; //survival knight's landing -> abv fremont weir
  vector[T]            phi_s1_t3_14; //survival blw fremont -> feather
  vector[T]            phi_s1_t4_14; //survival feather -> sacramento
  vector[T]            phi_s1_t5_14; //survival sacramento -> freeport
  vector[T]            phi_s1_t6_14; //survival freeport -> sutter/steam complex
  vector[T]            phi_s1_t7_14; //survival sutter/steam complex -> georgiana
  vector[T]            phi_s1_t8_14; //survival georgiana -> rio vista
  vector[T]            phi_s1_t9_14; //survival rio vista -> chipps

  vector[T]            phi_s3_t8_14; //survival sutter -> rio vista
  vector[T]            phi_s4_t8_14; //survival steam -> rio vista

  vector[T]            phi_s5_t8_14; //survival georgian -> interior delta
  vector[T]            phi_s5_t9_14; //survival interior delta -> chipps

  vector[T]            psi_1to3_t6_14; //prob entering sutter
  vector[T]            psi_1to4_t6_14; //prob entering steam
  vector[T]            psi_1to5_t7_14; //prob entering georgiana

  vector[T]            p_s1_t1_14; //p detect knight's landing
  vector[T]            p_s1_t2_14 = rep_vector(0, T); //no detection abv fremont weir in 2014
  vector[T]            p_s1_t3_14 = rep_vector(0, T); //no detection blw fremont in 2014
  vector[T]            p_s1_t4_14; //p detect feather
  vector[T]            p_s1_t5_14; //p detect sacramento
  vector[T]            p_s1_t6_14; //p detect freeport
  vector[T]            p_s1_t7_14; //p detect sac blw sutter/steam
  vector[T]            p_s1_t8_14; //p detect sac blw georgiana
  vector[T]            p_s1_t9_14; //p detect rio vista
  vector[T]            p_s1_t10_14 = rep_vector(0, T); //no detection chipps in 2014

  vector[T]            p_s3_t7_14; // p detect sutter
  vector[T]            p_s4_t7_14; // p detect steam

  vector[T]            p_s5_t8_14; // p detect georgiana
  vector[T]            p_s5_t9_14; // p detect interior delta

  row_vector[T]        alpha_s1_t0_14 = rep_row_vector(0, T);
  matrix[T, T]         alpha_s1_t1_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t2_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t3_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t4_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t5_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t6_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t7_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t8_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t9_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t10_14 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s3_t8_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s4_t8_14 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s5_t8_14 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s5_t9_14 = rep_matrix(0, T, T);

  row_vector[T]        lambda_t0s1_t1s1_14;
  row_vector[T]        lambda_t0s1_t4s1_14;
  row_vector[T]        lambda_t0s1_t5s1_14;

  matrix[T, T]         lambda_t1s1_t4s1_14;
  matrix[T, T]         lambda_t1s1_t5s1_14;

  matrix[T, T]         lambda_t4s1_t5s1_14;

  matrix[T, T]         lambda_t5s1_t6s1_14;
  matrix[T, T]         lambda_t5s1_t7s1_14;
  matrix[T, T]         lambda_t5s1_t8s1_14;
  matrix[T, T]         lambda_t5s1_t9s1_14;

  matrix[T, T]         lambda_t5s1_t7s3_14;
  matrix[T, T]         lambda_t5s1_t7s4_14;

  matrix[T, T]         lambda_t6s1_t7s1_14;
  matrix[T, T]         lambda_t6s1_t8s1_14;
  matrix[T, T]         lambda_t6s1_t9s1_14;
  matrix[T, T]         lambda_t6s1_t11s1_14;

  matrix[T, T]         lambda_t6s1_t7s3_14;
  matrix[T, T]         lambda_t6s1_t7s4_14;
  matrix[T, T]         lambda_t6s1_t8s5_14;
  matrix[T, T]         lambda_t6s1_t9s5_14;

  matrix[T, T]         lambda_t7s1_t8s1_14;
  matrix[T, T]         lambda_t7s1_t8s5_14;
  matrix[T, T]         lambda_t7s1_t9s1_14;
  
  matrix[T, T]         lambda_t7s3_t9s1_14;
  matrix[T, T]         lambda_t7s4_t9s1_14;

  matrix[T, T]         lambda_t8s1_t9s1_14;
  matrix[T, T]         lambda_t8s1_t11s1_14;


  matrix[T, T]         lambda_t9s1_t11s1_14;

  vector[T]            chi_9_s5_14; //prob last seen interior delta
  vector[T]            chi_9_s1_14; //prob last seen rio vista
  vector[T]            chi_8_s5_14; //prob last seen georgiana
  vector[T]            chi_8_s1_14; //prob last seen sac blw georgiana
  vector[T]            chi_7_s4_14; //prob last seen steam
  vector[T]            chi_7_s3_14; //prob last seen sutter
  vector[T]            chi_7_s1_14; //prob last seen sac blw sutter/steam
  vector[T]            chi_6_s1_14; //prob last seen freeport
  vector[T]            chi_5_s1_14; //prob last seen sacramento
  vector[T]            chi_4_s2_14; //prob last seen toe drain
  vector[T]            chi_4_s1_14; //prob last seen feather
  vector[T]            chi_1_s1_14; //prob last seen knight's landing
  real                 chi_0_s1_14; //prob last seen release
  
  //Tran Param 2015-----------------------------------------------------------------\\  
  vector[2]            phi_s1_t0_15; //survival release to knight's landing
  vector[T]            phi_s1_t1_15; //survival knight's landing -> abv fremont weir
  vector[T]            phi_s1_t3_15; //survival blw fremont -> feather
  vector[T]            phi_s1_t4_15; //survival feather -> sacramento
  vector[T]            phi_s1_t5_15; //survival sacramento -> freeport
  vector[T]            phi_s1_t6_15; //survival freeport -> sutter/steam complex
  vector[T]            phi_s1_t7_15; //survival sutter/steam complex -> georgiana
  vector[T]            phi_s1_t8_15; //survival georgiana -> rio vista
  vector[T]            phi_s1_t9_15; //survival rio vista -> chipps

  vector[T]            phi_s3_t8_15; //survival sutter -> rio vista
  vector[T]            phi_s4_t8_15; //survival steam -> rio vista

  vector[T]            phi_s5_t8_15; //survival georgian -> interior delta
  vector[T]            phi_s5_t9_15; //survival interior delta -> chipps

  vector[T]            psi_1to3_t6_15; //prob entering sutter
  vector[T]            psi_1to4_t6_15; //prob entering steam
  vector[T]            psi_1to5_t7_15; //prob entering georgiana

  vector[T]            p_s1_t1_15; //p detect knight's landing
  vector[T]            p_s1_t2_15 = rep_vector(0, T); //no detection abv fremont weir in 2015
  vector[T]            p_s1_t3_15 = rep_vector(0, T); //no detection blw fremont in 2015
  vector[T]            p_s1_t4_15; //p detect feather
  vector[T]            p_s1_t5_15; //p detect sacramento
  vector[T]            p_s1_t6_15; //p detect freeport
  vector[T]            p_s1_t7_15; //p detect sac blw sutter/steam
  vector[T]            p_s1_t8_15; //p detect sac blw georgiana
  vector[T]            p_s1_t9_15; //p detect rio vista
  vector[T]            p_s1_t10_15 = rep_vector(0, T); //no detection chipps in 2015

  vector[T]            p_s3_t7_15; // p detect sutter
  vector[T]            p_s4_t7_15; // p detect steam

  vector[T]            p_s5_t8_15; // p detect georgiana
  vector[T]            p_s5_t9_15; // p detect interior delta

  matrix[2, T]         alpha_s1_t0_15 = rep_matrix(0, 2, T);
  matrix[T, T]         alpha_s1_t1_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t2_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t3_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t4_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t5_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t6_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t7_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t8_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t9_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t10_15 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s3_t8_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s4_t8_15 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s5_t8_15 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s5_t9_15 = rep_matrix(0, T, T);

  matrix[2, T]         lambda_t0s1_t1s1_15;
  matrix[2, T]         lambda_t0s1_t4s1_15;
  matrix[2, T]         lambda_t0s1_t5s1_15;
  matrix[2, T]         lambda_t0s1_t6s1_15;
  matrix[2, T]         lambda_t0s1_t7s1_15;
  matrix[2, T]         lambda_t0s1_t8s1_15;

  matrix[2, T]         lambda_t0s1_t7s3_15;
  matrix[2, T]         lambda_t0s1_t7s4_15;

  matrix[T, T]         lambda_t1s1_t4s1_15;
  matrix[T, T]         lambda_t1s1_t5s1_15;
  matrix[T, T]         lambda_t1s1_t6s1_15;
  matrix[T, T]         lambda_t1s1_t7s1_15;
  matrix[T, T]         lambda_t1s1_t8s1_15;

  matrix[T, T]         lambda_t1s1_t7s3_15;
  matrix[T, T]         lambda_t1s1_t7s4_15;

  matrix[T, T]         lambda_t4s1_t5s1_15;
  matrix[T, T]         lambda_t4s1_t6s1_15;
  matrix[T, T]         lambda_t4s1_t7s1_15;
  matrix[T, T]         lambda_t4s1_t8s1_15;
  matrix[T, T]         lambda_t4s1_t9s1_15;

  matrix[T, T]         lambda_t4s1_t7s3_15;
  matrix[T, T]         lambda_t4s1_t7s4_15;
  matrix[T, T]         lambda_t4s1_t8s5_15;

  matrix[T, T]         lambda_t5s1_t6s1_15;
  matrix[T, T]         lambda_t5s1_t7s1_15;
  matrix[T, T]         lambda_t5s1_t8s1_15;
  matrix[T, T]         lambda_t5s1_t9s1_15;

  matrix[T, T]         lambda_t5s1_t7s3_15;
  matrix[T, T]         lambda_t5s1_t7s4_15;
  matrix[T, T]         lambda_t5s1_t8s5_15;

  matrix[T, T]         lambda_t6s1_t7s1_15;
  
  matrix[T, T]         lambda_t6s1_t7s3_15;
  matrix[T, T]         lambda_t6s1_t7s4_15;

  matrix[T, T]         lambda_t7s1_t8s1_15;
  matrix[T, T]         lambda_t7s1_t9s1_15;

  matrix[T, T]         lambda_t7s3_t9s1_15;
  matrix[T, T]         lambda_t7s3_t11s1_15;

  matrix[T, T]         lambda_t7s4_t9s1_15;
  matrix[T, T]         lambda_t7s4_t11s1_15;

  matrix[T, T]         lambda_t7s1_t8s5_15;

  matrix[T, T]         lambda_t8s1_t9s1_15;
  matrix[T, T]         lambda_t8s1_t11s1_15;

  matrix[T, T]         lambda_t8s5_t9s5_15;
  matrix[T, T]         lambda_t8s5_t11s1_15;

  matrix[T, T]         lambda_t9s1_t11s1_15;
  matrix[T, T]         lambda_t9s5_t11s1_15;

  vector[T]            chi_9_s5_15; //prob last seen interior delta
  vector[T]            chi_9_s1_15; //prob last seen rio vista
  vector[T]            chi_8_s5_15; //prob last seen georgiana
  vector[T]            chi_8_s1_15; //prob last seen sac blw georgiana
  vector[T]            chi_7_s4_15; //prob last seen steam
  vector[T]            chi_7_s3_15; //prob last seen sutter
  vector[T]            chi_7_s1_15; //prob last seen sac blw sutter/steam
  vector[T]            chi_6_s1_15; //prob last seen freeport
  vector[T]            chi_5_s1_15; //prob last seen sacramento
  vector[T]            chi_4_s2_15; //prob last seen toe drain
  vector[T]            chi_4_s1_15; //prob last seen feather
  vector[T]            chi_1_s1_15; //prob last seen knight's landing
  vector[2]            chi_0_s1_15; //prob last seen release
  
  //Tran Param 2016-----------------------------------------------------------------\\  
  vector[2]            phi_s1_t0_16; //survival release to knight's landing
  vector[T]            phi_s1_t1_16; //survival knight's landing -> abv fremont weir
  vector[T]            phi_s1_t3_16; //survival blw fremont -> feather
  vector[T]            phi_s1_t4_16; //survival feather -> sacramento
  vector[T]            phi_s1_t5_16; //survival sacramento -> freeport
  vector[T]            phi_s1_t6_16; //survival freeport -> sutter/steam complex
  vector[T]            phi_s1_t7_16; //survival sutter/steam complex -> georgiana
  vector[T]            phi_s1_t8_16; //survival georgiana -> rio vista
  vector[T]            phi_s1_t9_16; //survival rio vista -> chipps

  vector[T]            phi_s2_t3_16; //survival fremont weir -> toe drain
  vector[T]            phi_s2_t8_16; //survival yolo -> rio vista

  vector[T]            phi_s3_t8_16; //survival sutter -> rio vista
  vector[T]            phi_s4_t8_16; //survival steam -> rio vista

  vector[T]            phi_s5_t8_16; //survival georgian -> interior delta
  vector[T]            phi_s5_t9_16; //survival interior delta -> chipps

  vector[T]            psi_1to2_t2_16;
  vector[T]            psi_1to3_t6_16;
  vector[T]            psi_1to4_t6_16;
  vector[T]            psi_1to5_t7_16;

  vector[T]            p_s1_t1_16;
  vector[T]            p_s1_t2_16 = rep_vector(0, T); //no detection abv fremont weir in 2016
  vector[T]            p_s1_t3_16 = rep_vector(0, T); //no detection blw fremont weir in 2016
  vector[T]            p_s1_t4_16;
  vector[T]            p_s1_t5_16;
  vector[T]            p_s1_t6_16;
  vector[T]            p_s1_t7_16;
  vector[T]            p_s1_t8_16 = rep_vector(0, T); //no detection sac blw georgina in 2016
  vector[T]            p_s1_t9_16;
  vector[T]            p_s1_t10_16;

  vector[T]            p_s2_t4_16;

  vector[T]            p_s3_t7_16;
  vector[T]            p_s4_t7_16;

  vector[T]            p_s5_t8_16  = rep_vector(0, T); //no detection georgiana entrance in 2016
  vector[T]            p_s5_t9_16;

  matrix[2, T]         alpha_s1_t0_16 = rep_matrix(0, 2, T);
  matrix[T, T]         alpha_s1_t1_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t2_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t3_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t4_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t5_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t6_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t7_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t8_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t9_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t10_16 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s2_t3_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s2_t8_16 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s3_t8_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s4_t8_16 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s5_t8_16 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s5_t9_16 = rep_matrix(0, T, T);

  matrix[2, T]      lambda_t0s1_t1s1_16;
  matrix[2, T]      lambda_t0s1_t4s1_16;
  matrix[2, T]      lambda_t0s1_t5s1_16;
  matrix[2, T]      lambda_t0s1_t6s1_16;
  matrix[2, T]      lambda_t0s1_t7s1_16;
  matrix[2, T]      lambda_t0s1_t9s1_16;
  matrix[2, T]      lambda_t0s1_t10s1_16;

  matrix[2, T]      lambda_t0s1_t4s2_16;
  matrix[2, T]      lambda_t0s1_t7s3_16;
  matrix[2, T]      lambda_t0s1_t7s4_16;
  matrix[2, T]      lambda_t0s1_t9s5_16;

  matrix[T, T]         lambda_t1s1_t4s1_16;
  matrix[T, T]         lambda_t1s1_t5s1_16;
  matrix[T, T]         lambda_t1s1_t6s1_16;

  matrix[T, T]         lambda_t1s1_t4s2_16;
  matrix[T, T]         lambda_t1s1_t7s4_16;

  matrix[T, T]         lambda_t4s1_t5s1_16;
  matrix[T, T]         lambda_t4s1_t6s1_16;
  matrix[T, T]         lambda_t4s1_t7s1_16;

  matrix[T, T]         lambda_t4s2_t9s1_16;

  matrix[T, T]         lambda_t5s1_t6s1_16;
  matrix[T, T]         lambda_t5s1_t7s1_16;
  matrix[T, T]         lambda_t5s1_t9s1_16;

  matrix[T, T]         lambda_t5s1_t7s3_16;
  matrix[T, T]         lambda_t5s1_t7s4_16;

  matrix[T, T]         lambda_t6s1_t7s1_16;
  matrix[T, T]         lambda_t6s1_t9s1_16;

  matrix[T, T]         lambda_t6s1_t7s3_16;
  matrix[T, T]         lambda_t6s1_t7s4_16;

  matrix[T, T]         lambda_t7s1_t9s1_16;
  matrix[T, T]         lambda_t7s1_t9s5_16;
  matrix[T, T]         lambda_t7s1_t10s1_16;

  matrix[T, T]         lambda_t7s3_t9s1_16;
  matrix[T, T]         lambda_t7s3_t10s1_16;

  matrix[T, T]         lambda_t7s4_t9s1_16;
  matrix[T, T]         lambda_t7s4_t10s1_16;

  matrix[T, T]         lambda_t9s1_t10s1_16;
  matrix[T, T]         lambda_t9s1_t11s1_16;

  matrix[T, T]         lambda_t9s5_t10s1_16;
  matrix[T, T]         lambda_t9s5_t11s1_16;

  vector[T]             chi_10_s1_16; //prob last seen Chipps
  vector[T]             chi_9_s1_16; //prob last seen rio vista
  vector[T]             chi_9_s5_16; //prob last seen interior delta
  vector[T]             chi_7_s1_16; //prob last seen sac blw sutter/steam
  vector[T]             chi_7_s3_16; //prob last seen sutter
  vector[T]             chi_7_s4_16; //prob last seen steam
  vector[T]             chi_6_s1_16; //prob last seen freeport
  vector[T]             chi_5_s1_16; //prob last seen sacramento
  vector[T]             chi_4_s2_16; //prob last seen toe drain
  vector[T]             chi_4_s1_16; //prob last seen feather
  vector[T]             chi_1_s1_16; //prob last seen knight's landing
  vector[2]             chi_0_s1_16; //prob last seen release
  
  //Tran Param 2017-----------------------------------------------------------------\\  

  real                 phi_s1_t0_17; //survival release to knight's landing
  vector[T]            phi_s1_t1_17; //survival knight's landing -> abv fremont weir
  vector[T]            phi_s1_t3_17; //survival blw fremont -> feather
  vector[T]            phi_s1_t4_17; //survival feather -> sacramento
  vector[T]            phi_s1_t5_17; //survival sacramento -> freeport
  vector[T]            phi_s1_t6_17; //survival freeport -> sutter/steam complex
  vector[T]            phi_s1_t7_17; //survival sutter/steam complex -> georgiana
  vector[T]            phi_s1_t8_17; //survival georgiana -> rio vista
  vector[T]            phi_s1_t9_17; //survival rio vista -> chipps

  vector[T]            phi_s2_t3_17; //survival fremont weir -> toe drain
  vector[T]            phi_s2_t8_17; //survival yolo -> rio vista

  vector[T]            phi_s3_t8_17; //survival sutter -> rio vista
  vector[T]            phi_s4_t8_17; //survival steam -> rio vista

  vector[T]            phi_s5_t8_17; //survival georgian -> interior delta
  vector[T]            phi_s5_t9_17; //survival interior delta -> chipps

  vector[T]            psi_1to2_t2_17;
  vector[T]            psi_1to3_t6_17;
  vector[T]            psi_1to4_t6_17;
  vector[T]            psi_1to5_t7_17;

  vector[T]            p_s1_t1_17;
  vector[T]            p_s1_t2_17;
  vector[T]            p_s1_t3_17;
  vector[T]            p_s1_t4_17;
  vector[T]            p_s1_t5_17;
  vector[T]            p_s1_t6_17;
  vector[T]            p_s1_t7_17;
  vector[T]            p_s1_t8_17;
  vector[T]            p_s1_t9_17;
  vector[T]            p_s1_t10_17;

  vector[T]            p_s2_t4_17;

  vector[T]            p_s3_t7_17;
  vector[T]            p_s4_t7_17;

  vector[T]            p_s5_t8_17 = rep_vector(0, T);
  vector[T]            p_s5_t9_17;

  row_vector[T]        alpha_s1_t0_17 = rep_row_vector(0, T);
  matrix[T, T]         alpha_s1_t1_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t2_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t3_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t4_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t5_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t6_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t7_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t8_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t9_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t10_17 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s2_t3_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s2_t8_17 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s3_t8_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s4_t8_17 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s5_t8_17 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s5_t9_17 = rep_matrix(0, T, T);

  row_vector[T]        lambda_t0s1_t1s1_17;
  row_vector[T]        lambda_t0s1_t2s1_17;
  row_vector[T]        lambda_t0s1_t3s1_17;
  row_vector[T]        lambda_t0s1_t4s1_17;
  row_vector[T]        lambda_t0s1_t5s1_17;
  row_vector[T]        lambda_t0s1_t6s1_17;
  row_vector[T]        lambda_t0s1_t7s1_17;
  row_vector[T]        lambda_t0s1_t8s1_17;
  row_vector[T]        lambda_t0s1_t9s1_17;
  row_vector[T]        lambda_t0s1_t10s1_17;
  row_vector[T]        lambda_t0s1_t4s2_17;
  row_vector[T]        lambda_t0s1_t7s3_17;
  row_vector[T]        lambda_t0s1_t7s4_17;
  row_vector[T]        lambda_t0s1_t8s5_17;
  row_vector[T]        lambda_t0s1_t9s5_17;

  matrix[T, T]         lambda_t1s1_t2s1_17;
  matrix[T, T]         lambda_t1s1_t3s1_17;
  matrix[T, T]         lambda_t1s1_t4s1_17;
  matrix[T, T]         lambda_t1s1_t5s1_17;
  matrix[T, T]         lambda_t1s1_t6s1_17;
  matrix[T, T]         lambda_t1s1_t7s1_17;
  matrix[T, T]         lambda_t1s1_t8s1_17;
  matrix[T, T]         lambda_t1s1_t9s1_17;
  matrix[T, T]         lambda_t1s1_t4s2_17;
  matrix[T, T]         lambda_t1s1_t7s3_17;
  matrix[T, T]         lambda_t1s1_t7s4_17;


  matrix[T, T]         lambda_t2s1_t3s1_17;
  matrix[T, T]         lambda_t2s1_t4s1_17;
  matrix[T, T]         lambda_t2s1_t5s1_17;
  matrix[T, T]         lambda_t2s1_t6s1_17;
  matrix[T, T]         lambda_t2s1_t7s1_17;
  matrix[T, T]         lambda_t2s1_t8s1_17;
  matrix[T, T]         lambda_t2s1_t9s1_17;
  matrix[T, T]         lambda_t2s1_t10s1_17;
  matrix[T, T]         lambda_t2s1_t4s2_17;
  matrix[T, T]         lambda_t2s1_t7s3_17;
  matrix[T, T]         lambda_t2s1_t7s4_17;
  matrix[T, T]         lambda_t2s1_t8s5_17;
  matrix[T, T]         lambda_t2s1_t9s5_17;

  matrix[T, T]         lambda_t3s1_t4s1_17;
  matrix[T, T]         lambda_t3s1_t5s1_17;
  matrix[T, T]         lambda_t3s1_t6s1_17;
  matrix[T, T]         lambda_t3s1_t7s1_17;
  matrix[T, T]         lambda_t3s1_t8s1_17;
  matrix[T, T]         lambda_t3s1_t9s1_17;
  matrix[T, T]         lambda_t3s1_t10s1_17;
  matrix[T, T]         lambda_t3s1_t7s3_17;
  matrix[T, T]         lambda_t3s1_t7s4_17;
  matrix[T, T]         lambda_t3s1_t8s5_17;
  matrix[T, T]         lambda_t3s1_t9s5_17;

  matrix[T, T]         lambda_t4s1_t5s1_17;
  matrix[T, T]         lambda_t4s1_t6s1_17;
  matrix[T, T]         lambda_t4s1_t7s1_17;
 
  matrix[T, T]         lambda_t4s2_t9s1_17;
  matrix[T, T]         lambda_t4s2_t10s1_17;
  matrix[T, T]         lambda_t4s2_t11s1_17;

  matrix[T, T]         lambda_t5s1_t6s1_17;
  matrix[T, T]         lambda_t5s1_t7s1_17;
  matrix[T, T]         lambda_t5s1_t8s1_17;
  matrix[T, T]         lambda_t5s1_t9s1_17;
  matrix[T, T]         lambda_t5s1_t10s1_17;
  matrix[T, T]         lambda_t5s1_t7s3_17;
  matrix[T, T]         lambda_t5s1_t7s4_17;
  matrix[T, T]         lambda_t5s1_t8s5_17;
  matrix[T, T]         lambda_t5s1_t9s5_17;

  matrix[T, T]         lambda_t6s1_t7s1_17;
  matrix[T, T]         lambda_t6s1_t8s1_17;
  matrix[T, T]         lambda_t6s1_t9s1_17;
  matrix[T, T]         lambda_t6s1_t10s1_17;
  matrix[T, T]         lambda_t6s1_t7s3_17;
  matrix[T, T]         lambda_t6s1_t7s4_17;
  matrix[T, T]         lambda_t6s1_t8s5_17;
  matrix[T, T]         lambda_t6s1_t9s5_17;

  matrix[T, T]         lambda_t7s1_t8s1_17;
  matrix[T, T]         lambda_t7s1_t9s1_17;
  matrix[T, T]         lambda_t7s1_t10s1_17;
  matrix[T, T]         lambda_t7s3_t9s1_17;
  matrix[T, T]         lambda_t7s3_t10s1_17;
  matrix[T, T]         lambda_t7s3_t11s1_17;
  matrix[T, T]         lambda_t7s4_t9s1_17;
  matrix[T, T]         lambda_t7s4_t10s1_17;
  matrix[T, T]         lambda_t7s4_t11s1_17;
  matrix[T, T]         lambda_t7s1_t8s5_17;
  matrix[T, T]         lambda_t7s1_t9s5_17;

  matrix[T, T]         lambda_t8s1_t9s1_17;
  matrix[T, T]         lambda_t8s1_t10s1_17;


  matrix[T, T]         lambda_t9s1_t10s1_17;
  matrix[T, T]         lambda_t9s1_t11s1_17;

  matrix[T, T]         lambda_t9s5_t10s1_17;

  vector[T]             chi_10_s1_17; //prob last seen Chipps
  vector[T]             chi_9_s1_17; //prob last seen rio vista
  vector[T]             chi_9_s5_17; //prob last seen interior delta
  vector[T]             chi_8_s1_17; //prob last seen sac blw georgiana
  vector[T]             chi_8_s5_17; //prob last seen georgiana
  vector[T]             chi_7_s1_17; //prob last seen sac blw sutter/steam
  vector[T]             chi_7_s3_17; //prob last seen sutter
  vector[T]             chi_7_s4_17; //prob last seen steam
  vector[T]             chi_6_s1_17; //prob last seen freeport
  vector[T]             chi_5_s1_17; //prob last seen sacramento
  vector[T]             chi_4_s2_17; //prob last seen toe drain
  vector[T]             chi_4_s1_17; //prob last seen feather
  vector[T]             chi_3_s1_17; //prob last seen blw fremont weir
  vector[T]             chi_2_s1_17; //prob last seen abv fremont weir
  vector[T]             chi_1_s1_17; //prob last seen knight's landing
  real                  chi_0_s1_17; //prob last seen release
  
  //Tran Param 2018-----------------------------------------------------------------\\  
  vector[2]            phi_s1_t0_18; //survival release to knight's landing
  vector[T]            phi_s1_t1_18; //survival knight's landing -> abv fremont weir
  vector[T]            phi_s1_t3_18; //survival blw fremont -> feather
  vector[T]            phi_s1_t4_18; //survival feather -> sacramento
  vector[T]            phi_s1_t5_18; //survival sacramento -> freeport
  vector[T]            phi_s1_t6_18; //survival freeport -> sutter/steam complex
  vector[T]            phi_s1_t7_18; //survival sutter/steam complex -> georgiana
  vector[T]            phi_s1_t8_18; //survival georgiana -> rio vista
  vector[T]            phi_s1_t9_18; //survival rio vista -> chipps

  vector[T]            phi_s2_t3_18; //survival fremont weir -> toe drain
  vector[T]            phi_s2_t8_18; //survival yolo -> rio vista

  vector[T]            phi_s3_t8_18; //survival sutter -> rio vista
  vector[T]            phi_s4_t8_18; //survival steam -> rio vista

  vector[T]            phi_s5_t8_18; //survival georgian -> interior delta
  vector[T]            phi_s5_t9_18; //survival interior delta -> chipps

  vector[T]            psi_1to2_t2_18;
  vector[T]            psi_1to3_t6_18;
  vector[T]            psi_1to4_t6_18;
  vector[T]            psi_1to5_t7_18;

  vector[T]            p_s1_t1_18;
  vector[T]            p_s1_t2_18;
  vector[T]            p_s1_t3_18;
  vector[T]            p_s1_t4_18;
  vector[T]            p_s1_t5_18;
  vector[T]            p_s1_t6_18;
  vector[T]            p_s1_t7_18;
  vector[T]            p_s1_t8_18;
  vector[T]            p_s1_t9_18;
  vector[T]            p_s1_t10_18;

  vector[T]            p_s2_t4_18;

  vector[T]            p_s3_t7_18;
  vector[T]            p_s4_t7_18;

  vector[T]            p_s5_t8_18;
  vector[T]            p_s5_t9_18;

  matrix[2, T]         alpha_s1_t0_18 = rep_matrix(0, 2, T);
  matrix[T, T]         alpha_s1_t1_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t2_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t3_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t4_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t5_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t6_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t7_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t8_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t9_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s1_t10_18 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s2_t3_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s2_t8_18 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s3_t8_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s4_t8_18 = rep_matrix(0, T, T);

  matrix[T, T]         alpha_s5_t8_18 = rep_matrix(0, T, T);
  matrix[T, T]         alpha_s5_t9_18 = rep_matrix(0, T, T);

  matrix[2, T]         lambda_t0s1_t1s1_18;

  matrix[T, T]         lambda_t1s1_t2s1_18;
  matrix[T, T]         lambda_t1s1_t3s1_18;
  matrix[T, T]         lambda_t1s1_t4s1_18;
  matrix[T, T]         lambda_t1s1_t5s1_18;
  matrix[T, T]         lambda_t1s1_t6s1_18;
  matrix[T, T]         lambda_t1s1_t7s1_18;

  matrix[T, T]         lambda_t2s1_t3s1_18;
  matrix[T, T]         lambda_t2s1_t4s1_18;
  matrix[T, T]         lambda_t2s1_t5s1_18;
  matrix[T, T]         lambda_t2s1_t6s1_18;
  matrix[T, T]         lambda_t2s1_t7s1_18;
  matrix[T, T]         lambda_t2s1_t7s3_18;
  matrix[T, T]         lambda_t2s1_t7s4_18;  
  matrix[T, T]         lambda_t2s1_t8s1_18;
  matrix[T, T]         lambda_t2s1_t9s1_18;

  matrix[T, T]         lambda_t3s1_t4s1_18;
  matrix[T, T]         lambda_t3s1_t5s1_18;
  matrix[T, T]         lambda_t3s1_t6s1_18;
  matrix[T, T]         lambda_t3s1_t7s1_18;
  matrix[T, T]         lambda_t3s1_t7s3_18;
  matrix[T, T]         lambda_t3s1_t7s4_18;
  matrix[T, T]         lambda_t3s1_t8s1_18;
  matrix[T, T]         lambda_t3s1_t9s1_18;

  matrix[T, T]         lambda_t4s1_t5s1_18;
  matrix[T, T]         lambda_t4s1_t6s1_18;
  matrix[T, T]         lambda_t4s1_t7s1_18;
  matrix[T, T]         lambda_t4s1_t7s3_18;
  matrix[T, T]         lambda_t4s1_t7s4_18;

  matrix[T, T]         lambda_t5s1_t6s1_18;
  matrix[T, T]         lambda_t5s1_t7s1_18;
  matrix[T, T]         lambda_t5s1_t7s3_18;
  matrix[T, T]         lambda_t5s1_t7s4_18;  

  matrix[T, T]         lambda_t6s1_t7s1_18;
  matrix[T, T]         lambda_t6s1_t7s3_18;
  matrix[T, T]         lambda_t6s1_t7s4_18;
  matrix[T, T]         lambda_t6s1_t8s5_18;  

  matrix[T, T]         lambda_t7s1_t8s1_18;
  matrix[T, T]         lambda_t7s1_t8s5_18;
  matrix[T, T]         lambda_t7s1_t9s1_18;
  matrix[T, T]         lambda_t7s1_t9s5_18;  
  matrix[T, T]         lambda_t7s1_t10s1_18;
  
  matrix[T, T]         lambda_t7s3_t9s1_18;
  matrix[T, T]         lambda_t7s3_t10s1_18;
  
  matrix[T, T]         lambda_t7s4_t9s1_18;
  matrix[T, T]         lambda_t7s4_t10s1_18;

  matrix[T, T]         lambda_t8s1_t9s1_18;
  matrix[T, T]         lambda_t8s1_t10s1_18;

  matrix[T, T]         lambda_t8s5_t9s5_18;

  matrix[T, T]         lambda_t9s1_t10s1_18;
  matrix[T, T]         lambda_t9s1_t11s1_18;

  matrix[T, T]         lambda_t9s5_t10s1_18;

  vector[T]             chi_10_s1_18[2]; //prob last seen Chipps
  vector[T]             chi_9_s1_18[2]; //prob last seen rio vista
  vector[T]             chi_9_s5_18[2]; //prob last seen interior delta
  vector[T]             chi_8_s1_18[2]; //prob last seen sac blw georgiana
  vector[T]             chi_8_s5_18[2]; //prob last seen georgiana
  vector[T]             chi_7_s1_18[2]; //prob last seen sac blw sutter/steam
  vector[T]             chi_7_s3_18[2]; //prob last seen sutter
  vector[T]             chi_7_s4_18[2]; //prob last seen steam
  vector[T]             chi_6_s1_18[2]; //prob last seen freeport
  vector[T]             chi_5_s1_18[2]; //prob last seen sacramento
  vector[T]             chi_4_s2_18[2]; //prob last seen yolo
  vector[T]             chi_4_s1_18[2]; //prob last seen feather
  vector[T]             chi_3_s1_18[2]; //prob last seen blw fremont weir
  vector[T]             chi_2_s1_18[2]; //prob last seen abv fremont weir
  vector[T]             chi_1_s1_18[2]; //prob last seen knight's landing
  real                  chi_0_s1_18[2]; //prob last seen release
  
  //-----------------------------------------------------------------------\\
  //***********************************************************************\\
  //Start defining TP------------------------------------------------------\\

  phi_s1_t0_14      = inv_logit(beta_phi_s1_t0_14);
  
  alpha_s1_t0_14 = arrival_probs_re(T, mu_v_s1_t0_14, sigma_v_s1_t0_14, tt_eps_s1_t0_14 * sigma_eps_s1_t0);

  for (g in 1:2){
    phi_s1_t0_15[g]      = inv_logit(beta_phi_s1_t0_15);
    // alpha_s1_t0_15[g, rls_day_15[g]:T] = arrival_probs(T - rls_day_15[g] + 1, beta_mu_s1_t0_15, sigma_tt_s1_t0_15);
    alpha_s1_t0_15[g, rls_day_15[g]:T] = arrival_probs_re(T - rls_day_15[g] + 1,  mu_v_s1_t0_15, sigma_v_s1_t0_15, tt_eps_s1_t0_15[rls_day_15[g]:T] * sigma_eps_s1_t0);
  }

  for (g in 1:2){
    phi_s1_t0_16[g]      = inv_logit(beta_phi_s1_t0_16);
    alpha_s1_t0_16[g, rls_day_16[g]:T] = arrival_probs_re(T - rls_day_16[g] + 1, mu_v_s1_t0_16, sigma_v_s1_t0_16, tt_eps_s1_t0_16[rls_day_16[g]:T] * sigma_eps_s1_t0);
  }

  phi_s1_t0_17      = inv_logit(beta_phi_s1_t0_17);
  alpha_s1_t0_17   = arrival_probs_re(T, mu_v_s1_t0_17, sigma_v_s1_t0_17, tt_eps_s1_t0_17 * sigma_eps_s1_t0);

  for (g in 1:2){
    phi_s1_t0_18[g]      = inv_logit(beta_phi_s1_t0_18[g]);
    alpha_s1_t0_18[g, rls_day_18[g]:T] = arrival_probs_re(T - rls_day_18[g] + 1, mu_v_s1_t0_18[g], sigma_v_s1_t0_18, tt_eps_s1_t0_18[rls_day_18[g]:T] * sigma_eps_s1_t0);
  }
  
  for (t in 1:T){
    // 2014 Deriv Params ---- \\
    {
      phi_s1_t1_14[t]   = inv_logit(beta_phi_s1_t1[1] + beta_phi_s1_t1[2] * KL[t, 1] + beta_phi_s1_t1[3] * TMP[t, 1]);
      phi_s1_t3_14[t]   = inv_logit(beta_phi_s1_t3[1] + beta_phi_s1_t3[2] * KL[t, 1] + beta_phi_s1_t3[3] * TMP[t, 1]);
      phi_s1_t4_14[t]   = inv_logit(beta_phi_s1_t4[1] + beta_phi_s1_t4[2] * SAC[t, 1] + beta_phi_s1_t4[3] * TMP[t, 1]);
      phi_s1_t5_14[t]   = inv_logit(beta_phi_s1_t5[1] + beta_phi_s1_t5[2] * SAC[t, 1] + beta_phi_s1_t5[3] * TMP[t, 1]);
      phi_s1_t6_14[t]   = inv_logit(beta_phi_s1_t6[1] + beta_phi_s1_t6[2] * SAC[t, 1] + beta_phi_s1_t6[3] * TMP[t, 1]);
      phi_s1_t7_14[t]   = inv_logit(beta_phi_s1_t7[1] + beta_phi_s1_t7[2] * SAC[t, 1] + beta_phi_s1_t7[3] * TMP[t, 1]);
      phi_s1_t8_14[t]   = inv_logit(beta_phi_s1_t8[1] + beta_phi_s1_t8[2] * SAC[t, 1] + beta_phi_s1_t8[3] * TMP[t, 1]);
      phi_s1_t9_14[t]   = inv_logit(beta_phi_s1_t9[1] + beta_phi_s1_t9[2] * RIO[t, 1] + beta_phi_s1_t9[3] * TMP[t, 1]);

      phi_s3_t8_14[t]   = inv_logit(beta_phi_s3_t8[1] + beta_phi_s3_t8[2] * SAC[t, 1] + beta_phi_s3_t8[3] * TMP[t, 1]);
      phi_s4_t8_14[t]   = inv_logit(beta_phi_s4_t8[1] + beta_phi_s4_t8[2] * SAC[t, 1] + beta_phi_s4_t8[3] * TMP[t, 1]);

      phi_s5_t8_14[t]   = inv_logit(beta_phi_s5_t8[1] + beta_phi_s5_t8[2] * SAC[t, 1] + beta_phi_s5_t8[3] * TMP[t, 1]);
      phi_s5_t9_14[t]   = inv_logit(beta_phi_s5_t9[1] + beta_phi_s5_t9[2] * SAC[t, 1] +
                                    beta_phi_s5_t9[3] * EI[t, 1] + beta_phi_s5_t9[4] * TMP[t, 1]);

      p_s1_t1_14[t]     = inv_logit(beta_p_s1_t1_14[1] + beta_p_s1_t1_14[2] * KL[t, 1]);
      p_s1_t4_14[t]     = inv_logit(beta_p_s1_t4_14[1] + beta_p_s1_t4_14[2] * SAC[t, 1]);
      p_s1_t5_14[t]     = inv_logit(beta_p_s1_t5_14[1] + beta_p_s1_t5_14[2] * SAC[t, 1]);
      p_s1_t6_14[t]     = inv_logit(beta_p_s1_t6_14[1] + beta_p_s1_t6_14[2] * SAC[t, 1]);
      p_s1_t7_14[t]     = inv_logit(beta_p_s1_t7_14[1] + beta_p_s1_t7_14[2] * SAC[t, 1]);
      p_s1_t8_14[t]     = inv_logit(beta_p_s1_t8_14[1] + beta_p_s1_t8_14[2] * SAC[t, 1]);
      p_s1_t9_14[t]     = inv_logit(beta_p_s1_t9_14[1] + beta_p_s1_t9_14[2] * RIO[t, 1]);

      p_s3_t7_14[t]     = inv_logit(beta_p_s3_t7_14[1] + beta_p_s3_t7_14[2] * SAC[t, 1]);
      p_s4_t7_14[t]     = inv_logit(beta_p_s4_t7_14[1] + beta_p_s4_t7_14[2] * SAC[t, 1]);

      p_s5_t8_14[t]     = inv_logit(beta_p_s5_t8_14[1] + beta_p_s5_t8_14[2] * SAC[t, 1]);
      p_s5_t9_14[t]     = inv_logit(beta_p_s5_t9_14[1] + beta_p_s5_t9_14[2] * SAC[t, 1]);

      psi_1to3_t6_14[t] = inv_logit(beta_psi_1to3_t6[1] + beta_psi_1to3_t6[2] * SAC[t, 1]);
      psi_1to4_t6_14[t] = inv_logit(beta_psi_1to4_t6[1] + beta_psi_1to4_t6[2] * SAC[t, 1]);
      psi_1to5_t7_14[t] = inv_logit(beta_psi_1to5_t7[1] + beta_psi_1to5_t7[2] * SAC[t, 1]);

      alpha_s1_t1_14[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t1[1] + beta_mu_s1_t1[2] * KL[t, 1],
                              sigma_tt_s1_t1);
                              
      alpha_s1_t2_14[t, t:T] = arrival_probs(T - t + 1,
                                             beta_mu_s1_t2b[1] + beta_mu_s1_t2b[2] * KL[t, 1],
                                             sigma_tt_s1_t2b);
      alpha_s1_t3_14[t, t]   = 1; 

      alpha_s1_t4_14[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t4[1] + beta_mu_s1_t4[2] * SAC[t, 1],
                              sigma_tt_s1_t4);
      alpha_s1_t5_14[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t5[1] + beta_mu_s1_t5[2] * SAC[t, 1],
                              sigma_tt_s1_t5);
      alpha_s1_t6_14[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t6[1] + beta_mu_s1_t6[2] * SAC[t, 1],
                               sigma_tt_s1_t6);
      alpha_s1_t7_14[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t7[1] + beta_mu_s1_t7[2] * SAC[t, 1],
                               sigma_tt_s1_t7);
      alpha_s1_t8_14[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t8[1] + beta_mu_s1_t8[2] * SAC[t, 1],
                               sigma_tt_s1_t8);
      alpha_s1_t9_14[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t9[1] + beta_mu_s1_t9[2] * RIO[t, 1],
                               sigma_tt_s1_t9);
      alpha_s1_t10_14[t, t:T] = arrival_probs(T - t + 1, beta_mu_s1_t10, sigma_tt_s1_t10);


      alpha_s3_t8_14[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s3_t8[1] + beta_mu_s3_t8[2] * SAC[t, 1],
                               sigma_tt_s3_t8);
      alpha_s4_t8_14[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s4_t8[1] + beta_mu_s4_t8[2] * SAC[t, 1],
                               sigma_tt_s4_t8);
      alpha_s5_t8_14[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t8[1] + beta_mu_s5_t8[2] * SAC[t, 1],
                               sigma_tt_s5_t8);
      alpha_s5_t9_14[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t9[1] + beta_mu_s5_t9[2] * SAC[t, 1],
                               sigma_tt_s5_t9);
    }

    // 2015 Deriv Params ---- \\
    {
      phi_s1_t1_15[t]   = inv_logit(beta_phi_s1_t1[1] + beta_phi_s1_t1[2] * KL[t, 2] + beta_phi_s1_t1[3] * TMP[t, 2]);
      phi_s1_t3_15[t]   = inv_logit(beta_phi_s1_t3[1] + beta_phi_s1_t3[2] * KL[t, 2] + beta_phi_s1_t3[3] * TMP[t, 2]);
      phi_s1_t4_15[t]   = inv_logit(beta_phi_s1_t4[1] + beta_phi_s1_t4[2] * SAC[t, 2] + beta_phi_s1_t4[3] * TMP[t, 2]);
      phi_s1_t5_15[t]   = inv_logit(beta_phi_s1_t5[1] + beta_phi_s1_t5[2] * SAC[t, 2] + beta_phi_s1_t5[3] * TMP[t, 2]);
      phi_s1_t6_15[t]   = inv_logit(beta_phi_s1_t6[1] + beta_phi_s1_t6[2] * SAC[t, 2] + beta_phi_s1_t6[3] * TMP[t, 2]);
      phi_s1_t7_15[t]   = inv_logit(beta_phi_s1_t7[1] + beta_phi_s1_t7[2] * SAC[t, 2] + beta_phi_s1_t7[3] * TMP[t, 2]);
      phi_s1_t8_15[t]   = inv_logit(beta_phi_s1_t8[1] + beta_phi_s1_t8[2] * SAC[t, 2] + beta_phi_s1_t8[3] * TMP[t, 2]);
      phi_s1_t9_15[t]   = inv_logit(beta_phi_s1_t9[1] + beta_phi_s1_t9[2] * RIO[t, 2] + beta_phi_s1_t9[3] * TMP[t, 2]);

      phi_s3_t8_15[t]   = inv_logit(beta_phi_s3_t8[1] + beta_phi_s3_t8[2] * SAC[t, 2] + beta_phi_s3_t8[3] * TMP[t, 2]);
      phi_s4_t8_15[t]   = inv_logit(beta_phi_s4_t8[1] + beta_phi_s4_t8[2] * SAC[t, 2] + beta_phi_s4_t8[3] * TMP[t, 2]);

      phi_s5_t8_15[t]   = inv_logit(beta_phi_s5_t8[1] + beta_phi_s5_t8[2] * SAC[t, 2] + beta_phi_s5_t8[3] * TMP[t, 2]);
      phi_s5_t9_15[t]   = inv_logit(beta_phi_s5_t9[1] + beta_phi_s5_t9[2] * SAC[t, 2] +
                                    beta_phi_s5_t9[3] * EI[t, 2] + beta_phi_s5_t9[4] * TMP[t, 2]);
                                    
      p_s1_t1_15[t]     = inv_logit(beta_p_s1_t1_15[1] + beta_p_s1_t1_15[2] * KL[t, 2]);
      p_s1_t4_15[t]     = inv_logit(beta_p_s1_t4_15[1] + beta_p_s1_t4_15[2] * SAC[t, 2]);
      p_s1_t5_15[t]     = inv_logit(beta_p_s1_t5_15[1] + beta_p_s1_t5_15[2] * SAC[t, 2]);
      p_s1_t7_15[t]     = inv_logit(beta_p_s1_t7_15[1] + beta_p_s1_t7_15[2] * SAC[t, 2]);
      p_s1_t8_15[t]     = inv_logit(beta_p_s1_t8_15[1] + beta_p_s1_t8_15[2] * SAC[t, 2]);
      p_s1_t9_15[t]     = inv_logit(beta_p_s1_t9_15[1] + beta_p_s1_t9_15[2] * RIO[t, 2]);

      p_s3_t7_15[t]     = inv_logit(beta_p_s3_t7_15[1] + beta_p_s3_t7_15[2] * SAC[t, 2]);
      p_s4_t7_15[t]     = inv_logit(beta_p_s4_t7_15[1] + beta_p_s4_t7_15[2] * SAC[t, 2]);

      p_s5_t8_15[t]     = inv_logit(beta_p_s5_t8_15[1] + beta_p_s5_t8_15[2] * SAC[t, 2]);
      p_s5_t9_15[t]     = inv_logit(beta_p_s5_t9_15[1] + beta_p_s5_t9_15[2] * SAC[t, 2]);

      psi_1to3_t6_15[t] = inv_logit(beta_psi_1to3_t6[1] + beta_psi_1to3_t6[2] * SAC[t, 2]);
      psi_1to4_t6_15[t] = inv_logit(beta_psi_1to4_t6[1] + beta_psi_1to4_t6[2] * SAC[t, 2]);
      psi_1to5_t7_15[t] = inv_logit(beta_psi_1to5_t7[1] + beta_psi_1to5_t7[2] * SAC[t, 2]);

     alpha_s1_t1_15[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t1[1] + beta_mu_s1_t1[2] * KL[t, 2],
                              sigma_tt_s1_t1);
      alpha_s1_t2_15[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t2b[1] + beta_mu_s1_t2b[2] * KL[t, 2],
                              sigma_tt_s1_t2b);
      alpha_s1_t3_15[t, t]  = 1;               
      alpha_s1_t4_15[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t4[1] + beta_mu_s1_t4[2] * SAC[t, 2],
                              sigma_tt_s1_t4);
      alpha_s1_t5_15[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t5[1] + beta_mu_s1_t5[2] * SAC[t, 2],
                              sigma_tt_s1_t5);
      alpha_s1_t6_15[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t6[1] + beta_mu_s1_t6[2] * SAC[t, 2],
                               sigma_tt_s1_t6);
      alpha_s1_t7_15[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t7[1] + beta_mu_s1_t7[2] * SAC[t, 2],
                               sigma_tt_s1_t7);
      alpha_s1_t8_15[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t8[1] + beta_mu_s1_t8[2] * SAC[t, 2],
                               sigma_tt_s1_t8);
      alpha_s1_t9_15[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t9[1] + beta_mu_s1_t9[2] * RIO[t, 2],
                               sigma_tt_s1_t9);
      alpha_s1_t10_15[t, t:T] = arrival_probs(T - t + 1, beta_mu_s1_t10, sigma_tt_s1_t10);


      alpha_s3_t8_15[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s3_t8[1] + beta_mu_s3_t8[2] * SAC[t, 2],
                               sigma_tt_s3_t8);
      alpha_s4_t8_15[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s4_t8[1] + beta_mu_s4_t8[2] * SAC[t, 2],
                               sigma_tt_s4_t8);
      alpha_s5_t8_15[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t8[1] + beta_mu_s5_t8[2] * SAC[t, 2],
                               sigma_tt_s5_t8);
      alpha_s5_t9_15[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t9[1] + beta_mu_s5_t9[2] * SAC[t, 2],
                               sigma_tt_s5_t9);
                               
    }

    // 2016 Deriv Params ---- \\
    {
      phi_s1_t1_16[t]   = inv_logit(beta_phi_s1_t1[1] + beta_phi_s1_t1[2] * KL[t, 3] + beta_phi_s1_t1[3] * TMP[t, 3]);
      phi_s1_t3_16[t]   = inv_logit(beta_phi_s1_t3[1] + beta_phi_s1_t3[2] * KL[t, 3] + beta_phi_s1_t3[3] * TMP[t, 3]);
      phi_s1_t4_16[t]   = inv_logit(beta_phi_s1_t4[1] + beta_phi_s1_t4[2] * SAC[t, 3] + beta_phi_s1_t4[3] * TMP[t, 3]);
      phi_s1_t5_16[t]   = inv_logit(beta_phi_s1_t5[1] + beta_phi_s1_t5[2] * SAC[t, 3] + beta_phi_s1_t5[3] * TMP[t, 3]);
      phi_s1_t6_16[t]   = inv_logit(beta_phi_s1_t6[1] + beta_phi_s1_t6[2] * SAC[t, 3] + beta_phi_s1_t6[3] * TMP[t, 3]);
      phi_s1_t7_16[t]   = inv_logit(beta_phi_s1_t7[1] + beta_phi_s1_t7[2] * SAC[t, 3] + beta_phi_s1_t7[3] * TMP[t, 3]);
      phi_s1_t8_16[t]   = inv_logit(beta_phi_s1_t8[1] + beta_phi_s1_t8[2] * SAC[t, 3] + beta_phi_s1_t8[3] * TMP[t, 3]);
      phi_s1_t9_16[t]   = inv_logit(beta_phi_s1_t9[1] + beta_phi_s1_t9[2] * RIO[t, 3] + beta_phi_s1_t9[3] * TMP[t, 3]);

      phi_s2_t3_16[t]   = inv_logit(beta_phi_s2_t3[1] + beta_phi_s2_t3[2] * YOL[t, 1] + beta_phi_s2_t3[3] * YOL_T[t, 1]);
      phi_s2_t8_16[t]   = inv_logit(beta_phi_s2_t8[1] + beta_phi_s2_t8[2] * YOL[t, 1] + beta_phi_s2_t8[3] * YOL_T[t, 1]);

      phi_s3_t8_16[t]   = inv_logit(beta_phi_s3_t8[1] + beta_phi_s3_t8[2] * SAC[t, 3] + beta_phi_s3_t8[3] * TMP[t, 3]);
      phi_s4_t8_16[t]   = inv_logit(beta_phi_s4_t8[1] + beta_phi_s4_t8[2] * SAC[t, 3] + beta_phi_s4_t8[3] * TMP[t, 3]);

      phi_s5_t8_16[t]   = inv_logit(beta_phi_s5_t8[1] + beta_phi_s5_t8[2] * SAC[t, 3] + beta_phi_s5_t8[3] * TMP[t, 3]);
      phi_s5_t9_16[t]   = inv_logit(beta_phi_s5_t9[1] + beta_phi_s5_t9[2] * SAC[t, 3] +
                                    beta_phi_s5_t9[3] * EI[t, 3] + beta_phi_s5_t9[4] * TMP[t, 3]);
                                    
      p_s1_t1_16[t]     = inv_logit(beta_p_s1_t1_16[1] + beta_p_s1_t1_16[2] * KL[t, 3]);
      p_s1_t4_16[t]     = inv_logit(beta_p_s1_t4_16[1] + beta_p_s1_t4_16[2] * SAC[t, 3]);
      p_s1_t5_16[t]     = inv_logit(beta_p_s1_t5_16[1] + beta_p_s1_t5_16[2] * SAC[t, 3]);
      p_s1_t6_16[t]     = inv_logit(beta_p_s1_t6_16[1] + beta_p_s1_t6_16[2] * SAC[t, 3]);
      p_s1_t7_16[t]     = inv_logit(beta_p_s1_t7_16[1] + beta_p_s1_t7_16[2] * SAC[t, 3]);
      p_s1_t9_16[t]     = inv_logit(beta_p_s1_t9_16[1] + beta_p_s1_t9_16[2] * RIO[t, 3]);
      p_s1_t10_16[t]    = inv_logit(beta_p_s1_t10_16[1] + beta_p_s1_t10_16[2] * RIO[t, 3]);

      p_s2_t4_16[t]     = inv_logit(beta_p_s2_t4_16[1] + beta_p_s2_t4_16[2] * YOL[t, 3]);

      p_s3_t7_16[t]     = inv_logit(beta_p_s3_t7_16[1] + beta_p_s3_t7_16[2] * SAC[t, 3]);
      p_s4_t7_16[t]     = inv_logit(beta_p_s4_t7_16[1] + beta_p_s4_t7_16[2] * SAC[t, 3]);

      p_s5_t9_16[t]     = inv_logit(beta_p_s5_t9_16[1] + beta_p_s5_t9_16[2] * SAC[t, 3]);

      psi_1to2_t2_16[t] = YOL_O[t, 1] * inv_logit(beta_psi_1to2_t2[1] + beta_psi_1to2_t2[2] * YOL_R[t, 1]);
      psi_1to3_t6_16[t] = inv_logit(beta_psi_1to3_t6[1] + beta_psi_1to3_t6[2] * SAC[t, 3]);
      psi_1to4_t6_16[t] = inv_logit(beta_psi_1to4_t6[1] + beta_psi_1to4_t6[2] * SAC[t, 3]);
      psi_1to5_t7_16[t] = inv_logit(beta_psi_1to5_t7[1] + beta_psi_1to5_t7[2] * SAC[t, 3]);

      alpha_s1_t1_16[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t1[1] + beta_mu_s1_t1[2] * KL[t, 3],
                              sigma_tt_s1_t1);
      alpha_s1_t2_16[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t2b[1] + beta_mu_s1_t2b[2] * KL[t, 3],
                              sigma_tt_s1_t2b);
                                   
      alpha_s1_t3_16[t, t] = 1;

      alpha_s1_t4_16[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t4[1] + beta_mu_s1_t4[2] * SAC[t, 3],
                              sigma_tt_s1_t4);
      alpha_s1_t5_16[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t5[1] + beta_mu_s1_t5[2] * SAC[t, 3],
                              sigma_tt_s1_t5);
      alpha_s1_t6_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t6[1] + beta_mu_s1_t6[2] * SAC[t, 3],
                               sigma_tt_s1_t6);
      alpha_s1_t7_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t7[1] + beta_mu_s1_t7[2] * SAC[t, 3],
                               sigma_tt_s1_t7);
      alpha_s1_t8_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t8[1] + beta_mu_s1_t8[2] * SAC[t, 3],
                               sigma_tt_s1_t8);
      alpha_s1_t9_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t9[1] + beta_mu_s1_t9[2] * RIO[t, 3],
                               sigma_tt_s1_t9);
      alpha_s1_t10_16[t, t:T] = arrival_probs(T - t + 1, beta_mu_s1_t10, sigma_tt_s1_t10);

      alpha_s2_t3_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s2_t3[1] + beta_mu_s2_t3[2] * YOL[t, 3],
                               sigma_tt_s2_t3);
      alpha_s2_t8_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s2_t8[1] + beta_mu_s2_t8[2] * YOL[t, 3],
                               sigma_tt_s2_t8);

      alpha_s3_t8_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s3_t8[1] + beta_mu_s3_t8[2] * SAC[t, 3],
                               sigma_tt_s3_t8);
      alpha_s4_t8_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s4_t8[1] + beta_mu_s4_t8[2] * SAC[t, 3],
                               sigma_tt_s4_t8);
      alpha_s5_t8_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t8[1] + beta_mu_s5_t8[2] * SAC[t, 3],
                               sigma_tt_s5_t8);
      alpha_s5_t9_16[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t9[1] + beta_mu_s5_t9[2] * SAC[t, 3],
                               sigma_tt_s5_t9);
    }

    // 2017 Deriv Params ---- \\
    {
      phi_s1_t1_17[t]   = inv_logit(beta_phi_s1_t1[1] + beta_phi_s1_t1[2] * KL[t, 4] + beta_phi_s1_t1[3] * TMP[t, 4]);
      phi_s1_t3_17[t]   = inv_logit(beta_phi_s1_t3[1] + beta_phi_s1_t3[2] * KL[t, 4] + beta_phi_s1_t3[3] * TMP[t, 4]);
      phi_s1_t4_17[t]   = inv_logit(beta_phi_s1_t4[1] + beta_phi_s1_t4[2] * SAC[t, 4] + beta_phi_s1_t4[3] * TMP[t, 4]);
      phi_s1_t5_17[t]   = inv_logit(beta_phi_s1_t5[1] + beta_phi_s1_t5[2] * SAC[t, 4] + beta_phi_s1_t5[3] * TMP[t, 4]);
      phi_s1_t6_17[t]   = inv_logit(beta_phi_s1_t6[1] + beta_phi_s1_t6[2] * SAC[t, 4] + beta_phi_s1_t6[3] * TMP[t, 4]);
      phi_s1_t7_17[t]   = inv_logit(beta_phi_s1_t7[1] + beta_phi_s1_t7[2] * SAC[t, 4] + beta_phi_s1_t7[3] * TMP[t, 4]);
      phi_s1_t8_17[t]   = inv_logit(beta_phi_s1_t8[1] + beta_phi_s1_t8[2] * SAC[t, 4] + beta_phi_s1_t8[3] * TMP[t, 4]);
      phi_s1_t9_17[t]   = inv_logit(beta_phi_s1_t9[1] + beta_phi_s1_t9[2] * RIO[t, 4] + beta_phi_s1_t9[3] * TMP[t, 4]);
      
      phi_s2_t3_17[t]   = inv_logit(beta_phi_s2_t3[1] + beta_phi_s2_t3[2] * YOL[t, 2] + beta_phi_s2_t3[3] * YOL_T[t, 2]);
      phi_s2_t8_17[t]   = inv_logit(beta_phi_s2_t8[1] + beta_phi_s2_t8[2] * YOL[t, 2] + beta_phi_s2_t8[3] * YOL_T[t, 2]);

      phi_s3_t8_17[t]   = inv_logit(beta_phi_s3_t8[1] + beta_phi_s3_t8[2] * SAC[t, 4] + beta_phi_s3_t8[3] * TMP[t, 4]);
      phi_s4_t8_17[t]   = inv_logit(beta_phi_s4_t8[1] + beta_phi_s4_t8[2] * SAC[t, 4] + beta_phi_s4_t8[3] * TMP[t, 4]);

      phi_s5_t8_17[t]   = inv_logit(beta_phi_s5_t8[1] + beta_phi_s5_t8[2] * SAC[t, 4] + beta_phi_s5_t8[3] * TMP[t, 4]);
      phi_s5_t9_17[t]   =  inv_logit(beta_phi_s5_t9[1] + beta_phi_s5_t9[2] * SAC[t, 4] +
                                    beta_phi_s5_t9[3] * EI[t, 4] + beta_phi_s5_t9[4] * TMP[t, 4]);

      p_s1_t1_17[t]     = inv_logit(beta_p_s1_t1_17[1] + beta_p_s1_t1_17[2] * KL[t, 4]);
      p_s1_t2_17[t]     = inv_logit(beta_p_s1_t2_17[1] + beta_p_s1_t2_17[2] * KL[t, 4]);
      p_s1_t3_17[t]     = inv_logit(beta_p_s1_t3_17[1] + beta_p_s1_t3_17[2] * SAC[t, 4]);
      // p_s1_t4_17[t]     = inv_logit(beta_p_s1_t4_17[1] + beta_p_s1_t4_17[2] * SAC[t, 4]);
      p_s1_t5_17[t]     = inv_logit(beta_p_s1_t5_17[1] + beta_p_s1_t5_17[2] * SAC[t, 4]);
      p_s1_t6_17[t]     = inv_logit(beta_p_s1_t6_17[1] + beta_p_s1_t6_17[2] * SAC[t, 4]);
      p_s1_t7_17[t]     = inv_logit(beta_p_s1_t7_17[1] + beta_p_s1_t7_17[2] * SAC[t, 4]);
      p_s1_t9_17[t]     = inv_logit(beta_p_s1_t9_17[1] + beta_p_s1_t9_17[2] * RIO[t, 4]);
      p_s1_t10_17[t]    = inv_logit(beta_p_s1_t10_17[1] + beta_p_s1_t10_17[2] * RIO[t, 4]);

      p_s2_t4_17[t]     = inv_logit(beta_p_s2_t4_17[1] + beta_p_s2_t4_17[2] * YOL[t, 2]);

      // p_s3_t7_17[t]     = inv_logit(beta_p_s3_t7_17[1] + beta_p_s3_t7_17[2] * SAC[t, 4]);
      p_s4_t7_17[t]     = inv_logit(beta_p_s4_t7_17[1] + beta_p_s4_t7_17[2] * SAC[t, 4]);
      p_s5_t9_17[t]     = inv_logit(beta_p_s5_t9_17[1] + beta_p_s5_t9_17[2] * SAC[t, 4]);

      psi_1to2_t2_17[t] = YOL_O[t, 2] * inv_logit(beta_psi_1to2_t2[1] + beta_psi_1to2_t2[2] * YOL_R[t, 2]);
      psi_1to3_t6_17[t] = inv_logit(beta_psi_1to3_t6[1] + beta_psi_1to3_t6[2] * SAC[t, 4]);
      psi_1to4_t6_17[t] = inv_logit(beta_psi_1to4_t6[1] + beta_psi_1to4_t6[2] * SAC[t, 4]);
      psi_1to5_t7_17[t] = inv_logit(beta_psi_1to5_t7[1] + beta_psi_1to5_t7[2] * SAC[t, 4]);

      alpha_s1_t1_17[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t1[1] + beta_mu_s1_t1[2] * KL[t, 4],
                              sigma_tt_s1_t1);
                              
      alpha_s1_t2_17[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t2[1] + beta_mu_s1_t2[2] * KL[t, 4],
                              sigma_tt_s1_t2);
                              
      alpha_s1_t3_17[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t3[1] + beta_mu_s1_t3[2] * KL[t, 4],
                              sigma_tt_s1_t3);
                              
      // alpha_s1_t2_17[t, t:T] = arrival_probs(T - t + 1,
      //                         beta_mu_s1_t2[1] + beta_mu_s1_t2[2] * KL[t, 4], 
      //                         sigma_tt_s1_t2);
      // alpha_s1_t3_17[t, t:T] = arrival_probs(T - t + 1,
      //                         beta_mu_s1_t3[1] + beta_mu_s1_t3[2] * KL[t, 4],
      //                         sigma_tt_s1_t3);
      alpha_s1_t4_17[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t4[1] + beta_mu_s1_t4[2] * SAC[t, 4],
                              sigma_tt_s1_t4);
      alpha_s1_t5_17[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t5[1] + beta_mu_s1_t5[2] * SAC[t, 4],
                              sigma_tt_s1_t5);
      alpha_s1_t6_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t6[1] + beta_mu_s1_t6[2] * SAC[t, 4],
                               sigma_tt_s1_t6);
      alpha_s1_t7_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t7[1] + beta_mu_s1_t7[2] * SAC[t, 4],
                               sigma_tt_s1_t7);
      alpha_s1_t8_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t8[1] + beta_mu_s1_t8[2] * SAC[t, 4],
                               sigma_tt_s1_t8);
      alpha_s1_t9_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t9[1] + beta_mu_s1_t9[2] * RIO[t, 4],
                               sigma_tt_s1_t9);
      alpha_s1_t10_17[t, t:T] = arrival_probs(T - t + 1, beta_mu_s1_t10, sigma_tt_s1_t10);

      alpha_s2_t3_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s2_t3[1] + beta_mu_s2_t3[2] * YOL[t, 2],
                               sigma_tt_s2_t3);
      alpha_s2_t8_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s2_t8[1] + beta_mu_s2_t8[2] * YOL[t, 2],
                               sigma_tt_s2_t8);

      alpha_s3_t8_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s3_t8[1] + beta_mu_s3_t8[2] * SAC[t, 4],
                               sigma_tt_s3_t8);
      alpha_s4_t8_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s4_t8[1] + beta_mu_s4_t8[2] * SAC[t, 4],
                               sigma_tt_s4_t8);
      alpha_s5_t8_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t8[1] + beta_mu_s5_t8[2] * SAC[t, 4],
                               sigma_tt_s5_t8);
      alpha_s5_t9_17[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t9[1] + beta_mu_s5_t9[2] * SAC[t, 4],
                               sigma_tt_s5_t9);
    }
  
    // 2018 Deriv Params ---- \\
    {
      phi_s1_t1_18[t]   = inv_logit(beta_phi_s1_t1[1] + beta_phi_s1_t1[2] * KL[t, 5] + beta_phi_s1_t1[3] * TMP[t, 5]);
      phi_s1_t3_18[t]   = inv_logit(beta_phi_s1_t3[1] + beta_phi_s1_t3[2] * KL[t, 5] + beta_phi_s1_t3[3] * TMP[t, 5]);
      phi_s1_t4_18[t]   = inv_logit(beta_phi_s1_t4[1] + beta_phi_s1_t4[2] * SAC[t, 5] + beta_phi_s1_t4[3] * TMP[t, 5]);
      phi_s1_t5_18[t]   = inv_logit(beta_phi_s1_t5[1] + beta_phi_s1_t5[2] * SAC[t, 5] + beta_phi_s1_t5[3] * TMP[t, 5]);
      phi_s1_t6_18[t]   = inv_logit(beta_phi_s1_t6[1] + beta_phi_s1_t6[2] * SAC[t, 5] + beta_phi_s1_t6[3] * TMP[t, 5]);
      phi_s1_t7_18[t]   = inv_logit(beta_phi_s1_t7[1] + beta_phi_s1_t7[2] * SAC[t, 5] + beta_phi_s1_t7[3] * TMP[t, 5]);
      phi_s1_t8_18[t]   = inv_logit(beta_phi_s1_t8[1] + beta_phi_s1_t8[2] * SAC[t, 5] + beta_phi_s1_t8[3] * TMP[t, 5]);
      phi_s1_t9_18[t]   = inv_logit(beta_phi_s1_t9[1] + beta_phi_s1_t9[2] * RIO[t, 5] + beta_phi_s1_t9[3] * TMP[t, 5]);
      
      phi_s2_t3_18[t]   = inv_logit(beta_phi_s2_t3[1] + beta_phi_s2_t3[2] * YOL[t, 3] + beta_phi_s2_t3[2] * YOL_T[t, 3]);
      phi_s2_t8_18[t]   = inv_logit(beta_phi_s2_t8[1] + beta_phi_s2_t8[2] * YOL[t, 3] + beta_phi_s2_t8[2] * YOL_T[t, 3]);

      phi_s3_t8_18[t]   = inv_logit(beta_phi_s3_t8[1] + beta_phi_s3_t8[2] * SAC[t, 5] + beta_phi_s3_t8[3] * TMP[t, 5]);
      phi_s4_t8_18[t]   = inv_logit(beta_phi_s4_t8[1] + beta_phi_s4_t8[2] * SAC[t, 5] + beta_phi_s4_t8[3] * TMP[t, 5]);

      phi_s5_t8_18[t]   = inv_logit(beta_phi_s5_t8[1] + beta_phi_s5_t8[2] * SAC[t, 5] + beta_phi_s5_t8[3] * TMP[t, 5]);
      phi_s5_t9_18[t]   = inv_logit(beta_phi_s5_t9[1] + beta_phi_s5_t9[2] * SAC[t, 5] +
                                    beta_phi_s5_t9[3] * EI[t, 5] + beta_phi_s5_t9[4] * TMP[t, 5]);
                                    
      p_s1_t1_18[t]     = inv_logit(beta_p_s1_t1_18[1] + beta_p_s1_t1_18[2] * KL[t, 5]);
      p_s1_t2_18[t]     = inv_logit(beta_p_s1_t2_18[1] + beta_p_s1_t2_18[2] * KL[t, 5]);
      p_s1_t3_18[t]     = inv_logit(beta_p_s1_t3_18[1] + beta_p_s1_t3_18[2] * SAC[t, 5]);
      p_s1_t4_18[t]     = inv_logit(beta_p_s1_t4_18[1] + beta_p_s1_t4_18[2] * SAC[t, 5]);
      p_s1_t5_18[t]     = inv_logit(beta_p_s1_t5_18[1] + beta_p_s1_t5_18[2] * SAC[t, 5]);
      p_s1_t6_18[t]     = inv_logit(beta_p_s1_t6_18[1] + beta_p_s1_t6_18[2] * SAC[t, 5]);
      p_s1_t7_18[t]     = inv_logit(beta_p_s1_t7_18[1] + beta_p_s1_t7_18[2] * SAC[t, 5]);
      p_s1_t8_18[t]     = inv_logit(beta_p_s1_t8_18[1] + beta_p_s1_t8_18[2] * SAC[t, 5]);
      p_s1_t9_18[t]     = inv_logit(beta_p_s1_t9_18[1] + beta_p_s1_t9_18[2] * RIO[t, 5]);
      p_s1_t10_18[t]    = inv_logit(beta_p_s1_t10_18[1] + beta_p_s1_t10_18[2] * RIO[t, 5]);

      p_s2_t4_18[t]     = inv_logit(beta_p_s2_t4_18[1] + beta_p_s2_t4_18[2] * YOL[t, 3]);

      p_s3_t7_18[t]     = inv_logit(beta_p_s3_t7_18[1] + beta_p_s3_t7_18[2] * SAC[t, 5]);
      p_s4_t7_18[t]     = inv_logit(beta_p_s4_t7_18[1] + beta_p_s4_t7_18[2] * SAC[t, 5]);

      p_s5_t8_18[t]     = inv_logit(beta_p_s5_t8_18[1] + beta_p_s5_t8_18[2] * SAC[t, 5]);
      p_s5_t9_18[t]     = inv_logit(beta_p_s5_t9_18[1] + beta_p_s5_t9_18[2] * SAC[t, 5]);

      psi_1to2_t2_18[t] = YOL_O[t, 3] * inv_logit(beta_psi_1to2_t2[1] + beta_psi_1to2_t2[2] * YOL_R[t, 3]);
      psi_1to3_t6_18[t] = inv_logit(beta_psi_1to3_t6[1] + beta_psi_1to3_t6[2] * SAC[t, 5]);
      psi_1to4_t6_18[t] = inv_logit(beta_psi_1to4_t6[1] + beta_psi_1to4_t6[2] * SAC[t, 5]);
      psi_1to5_t7_18[t] = inv_logit(beta_psi_1to5_t7[1] + beta_psi_1to5_t7[2] * SAC[t, 5]);

      alpha_s1_t1_18[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t1[1] + beta_mu_s1_t1[2] * KL[t, 5],
                              sigma_tt_s1_t1);
                              
      alpha_s1_t2_18[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t2[1] + beta_mu_s1_t2[2] * KL[t, 5],
                              sigma_tt_s1_t2);
                              
      alpha_s1_t3_18[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t3[1] + beta_mu_s1_t3[2] * KL[t, 5],
                              sigma_tt_s1_t3);                        
      // alpha_s1_t2_18[t, t:T] = arrival_probs(T - t + 1,
      //                         beta_mu_s1_t2[1] + beta_mu_s1_t2[2] * KL[t, 5],
      //                         sigma_tt_s1_t2);
      // alpha_s1_t3_18[t, t:T] = arrival_probs(T - t + 1,
      //                         beta_mu_s1_t3[1] + beta_mu_s1_t3[2] * KL[t, 5],
      //                         sigma_tt_s1_t3);
      alpha_s1_t4_18[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t4[1] + beta_mu_s1_t4[2] * SAC[t, 5],
                              sigma_tt_s1_t4);
      alpha_s1_t5_18[t, t:T] = arrival_probs(T - t + 1,
                              beta_mu_s1_t5[1] + beta_mu_s1_t5[2] * SAC[t, 5],
                              sigma_tt_s1_t5);
      alpha_s1_t6_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t6[1] + beta_mu_s1_t6[2] * SAC[t, 5],
                               sigma_tt_s1_t6);
      alpha_s1_t7_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t7[1] + beta_mu_s1_t7[2] * SAC[t, 5],
                               sigma_tt_s1_t7);
      alpha_s1_t8_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t8[1] + beta_mu_s1_t8[2] * SAC[t, 5],
                               sigma_tt_s1_t8);
      alpha_s1_t9_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s1_t9[1] + beta_mu_s1_t9[2] * RIO[t, 5],
                               sigma_tt_s1_t9);
      alpha_s1_t10_18[t, t:T] = arrival_probs(T - t + 1, beta_mu_s1_t10, sigma_tt_s1_t10);

      alpha_s2_t3_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s2_t3[1] + beta_mu_s2_t3[2] * YOL[t, 3],
                               sigma_tt_s2_t3);
      alpha_s2_t8_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s2_t8[1] + beta_mu_s2_t8[2] * YOL[t, 3],
                               sigma_tt_s2_t8);

      alpha_s3_t8_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s3_t8[1] + beta_mu_s3_t8[2] * SAC[t, 5],
                               sigma_tt_s3_t8);
      alpha_s4_t8_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s4_t8[1] + beta_mu_s4_t8[2] * SAC[t, 5],
                               sigma_tt_s4_t8);
      alpha_s5_t8_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t8[1] + beta_mu_s5_t8[2] * SAC[t, 5],
                               sigma_tt_s5_t8);
      alpha_s5_t9_18[t, t:T] = arrival_probs(T - t + 1,
                               beta_mu_s5_t9[1] + beta_mu_s5_t9[2] * SAC[t, 5],
                               sigma_tt_s5_t9);
    }    
  }


  // Define exclusion to account for receiver outages
  for (t in 1:20){
    p_s1_t6_15[t]     = 0.0;
  }

  for (t in 21:T){
    p_s1_t6_15[t]     = inv_logit(beta_p_s1_t6_15[1] + beta_p_s1_t6_15[2] * SAC[t, 2]);
  }
  
  for (t in 1:53){
    p_s1_t4_17[t]     = 0.0;
  }

  for (t in 54:T){
    p_s1_t4_17[t]     = inv_logit(beta_p_s1_t4_17[1] + beta_p_s1_t4_17[2] * SAC[t, 4]);
  }
  
  for (t in 1:34){
    p_s3_t7_17[t]     = 0.0;
  }

  for (t in 35:T){
    p_s3_t7_17[t]     = inv_logit(beta_p_s3_t7_17[1] + beta_p_s3_t7_17[2] * SAC[t, 4]);
  }

  for (t in 1:40){
    p_s1_t8_17[t]     = inv_logit(beta_p_s1_t8_17[1] + beta_p_s1_t8_17[2] * SAC[t, 4]);
  }

  for (t in 41:T){
    p_s1_t8_17[t]     = 0.0;
  }
  
  //Below Fremont Weir out from 3/23-3/29
  for (t in 23:28){
    p_s1_t3_18[t]    = 0.0;
  }
  
  //Above and Below Fremont Weir out from 4/4 - 4/6
  for (t in 35:37){
    p_s1_t2_18[t]    = 0.0;
    p_s1_t3_18[t]    = 0.0;
  }
  
  // Below Goergiana out from 3/22 to 4/16
  for (t in 22:47){
    p_s1_t8_18[t]     = 0.0;
  }
  
 
  //lambdas 2014------------------------------------------------------------------------------------------\\
  {
    row_vector[T]      lambda_t0s1_t2s1_14;
    row_vector[T]      lambda_t0s1_t3s1_14;

    lambda_t0s1_t1s1_14 = phi_s1_t0_14 * alpha_s1_t0_14;
    lambda_t0s1_t2s1_14 = (lambda_t0s1_t1s1_14 .* ((1 - p_s1_t1_14) .* phi_s1_t1_14)') * alpha_s1_t1_14;
    lambda_t0s1_t3s1_14 = lambda_t0s1_t2s1_14 * alpha_s1_t2_14;
    lambda_t0s1_t4s1_14 = (lambda_t0s1_t3s1_14 .* phi_s1_t3_14') * alpha_s1_t3_14;
    lambda_t0s1_t5s1_14 = (lambda_t0s1_t4s1_14 .* ((1 - p_s1_t4_14) .* phi_s1_t4_14)') * alpha_s1_t4_14;
  }

  {
    lambda_t1s1_t4s1_14  = diag_post_multiply(diag_pre_multiply(phi_s1_t1_14, alpha_s1_t1_14) * alpha_s1_t2_14, phi_s1_t3_14) * alpha_s1_t3_14;
    lambda_t1s1_t5s1_14  = diag_post_multiply(lambda_t1s1_t4s1_14, (1 - p_s1_t4_14) .* phi_s1_t4_14) * alpha_s1_t4_14;
  }

  {
    lambda_t4s1_t5s1_14  = diag_pre_multiply(phi_s1_t4_14, alpha_s1_t4_14);
  }

  {
    lambda_t5s1_t6s1_14  = diag_pre_multiply(phi_s1_t5_14, alpha_s1_t5_14);
    lambda_t5s1_t7s1_14  = diag_post_multiply(lambda_t5s1_t6s1_14, (1 - p_s1_t6_14) .* phi_s1_t6_14) *
                           diag_post_multiply(alpha_s1_t6_14, (1 - psi_1to3_t6_14) .* (1 - psi_1to4_t6_14));
    lambda_t5s1_t7s3_14  = diag_post_multiply(lambda_t5s1_t6s1_14, (1 - p_s1_t6_14) .* phi_s1_t6_14) *
                           diag_post_multiply(alpha_s1_t6_14, psi_1to3_t6_14);
    lambda_t5s1_t7s4_14  = diag_post_multiply(lambda_t5s1_t6s1_14, (1 - p_s1_t6_14) .* phi_s1_t6_14) *
                           diag_post_multiply(alpha_s1_t6_14, (1 - psi_1to3_t6_14) .* psi_1to4_t6_14);
    lambda_t5s1_t8s1_14  = diag_post_multiply(lambda_t5s1_t7s1_14, (1 - p_s1_t7_14) .* phi_s1_t7_14) *
                           diag_post_multiply(alpha_s1_t7_14, (1 - psi_1to5_t7_14));
    lambda_t5s1_t9s1_14  = (diag_post_multiply(lambda_t5s1_t8s1_14, (1 - p_s1_t8_14) .* phi_s1_t8_14) * alpha_s1_t8_14) +
                           (diag_post_multiply(lambda_t5s1_t7s3_14, (1 - p_s3_t7_14) .* phi_s3_t8_14) * alpha_s3_t8_14) +
                           (diag_post_multiply(lambda_t5s1_t7s4_14, (1 - p_s4_t7_14) .* phi_s4_t8_14) * alpha_s4_t8_14);
  }

  {
    lambda_t6s1_t7s1_14  = diag_pre_multiply(phi_s1_t6_14,
                           diag_post_multiply(alpha_s1_t6_14, ((1 - psi_1to3_t6_14) .* (1 - psi_1to4_t6_14))));
    lambda_t6s1_t7s3_14  = diag_pre_multiply(phi_s1_t6_14,
                           diag_post_multiply(alpha_s1_t6_14, psi_1to3_t6_14));
    lambda_t6s1_t7s4_14  = diag_pre_multiply(phi_s1_t6_14,
                           diag_post_multiply(alpha_s1_t6_14, (1 - psi_1to3_t6_14) .* psi_1to4_t6_14));
    lambda_t6s1_t8s1_14  = diag_post_multiply(lambda_t6s1_t7s1_14, (1 - p_s1_t7_14) .* phi_s1_t7_14) *
                           diag_post_multiply(alpha_s1_t7_14, (1 - psi_1to5_t7_14));

    lambda_t6s1_t9s1_14  = (diag_post_multiply(lambda_t6s1_t8s1_14, (1 - p_s1_t8_14) .* phi_s1_t8_14) * alpha_s1_t8_14) +
                           (diag_post_multiply(lambda_t6s1_t7s3_14, (1 - p_s3_t7_14) .* phi_s3_t8_14) * alpha_s3_t8_14) +
                           (diag_post_multiply(lambda_t6s1_t7s4_14, (1 - p_s4_t7_14) .* phi_s4_t8_14) * alpha_s4_t8_14);
  }

  {

    lambda_t7s1_t8s1_14  = diag_pre_multiply(phi_s1_t7_14,
                           diag_post_multiply(alpha_s1_t7_14, (1 - psi_1to5_t7_14)));
    lambda_t7s1_t8s5_14  = diag_pre_multiply(phi_s1_t7_14,
                           diag_post_multiply(alpha_s1_t7_14,  psi_1to5_t7_14));
    lambda_t7s1_t9s1_14  = diag_post_multiply(lambda_t7s1_t8s1_14, (1 - p_s1_t8_14) .* phi_s1_t8_14) * alpha_s1_t8_14;
  
    lambda_t7s3_t9s1_14  = diag_pre_multiply(phi_s3_t8_14, alpha_s3_t8_14);
    lambda_t7s4_t9s1_14  = diag_pre_multiply(phi_s4_t8_14, alpha_s4_t8_14);
  }

  {
    matrix[T, T]      lambda_t8s1_t10s1_14;
    
    lambda_t8s1_t9s1_14  = diag_pre_multiply(phi_s1_t8_14, alpha_s1_t8_14);
    lambda_t8s1_t10s1_14 = diag_post_multiply(lambda_t8s1_t9s1_14, (1 - p_s1_t9_14) .* phi_s1_t9_14) * alpha_s1_t9_14;
    lambda_t8s1_t11s1_14 = diag_post_multiply(lambda_t8s1_t10s1_14, (1 - p_s1_t10_14) * chipps2benecia) * alpha_s1_t10_14;
  }

  {
    lambda_t9s1_t11s1_14 = diag_post_multiply(diag_pre_multiply(phi_s1_t9_14, alpha_s1_t9_14), (1 - p_s1_t10_14) * chipps2benecia) * alpha_s1_t10_14;
  }
  
  //lambdas 2015------------------------------------------------------------------------------------------\\
  {
    matrix[2, T]      lambda_t0s1_t2s1_15;
    matrix[2, T]      lambda_t0s1_t3s1_15;

    lambda_t0s1_t1s1_15 = diag_pre_multiply(phi_s1_t0_15, alpha_s1_t0_15);
    lambda_t0s1_t2s1_15 = diag_post_multiply(lambda_t0s1_t1s1_15, (1 - p_s1_t1_15) .* phi_s1_t1_15) * alpha_s1_t1_15;
    lambda_t0s1_t3s1_15 = lambda_t0s1_t2s1_15 * alpha_s1_t2_15;
    lambda_t0s1_t4s1_15 = diag_post_multiply(lambda_t0s1_t3s1_15, phi_s1_t3_15) * alpha_s1_t3_15;

    lambda_t0s1_t5s1_15 = diag_post_multiply(lambda_t0s1_t4s1_15, (1 - p_s1_t4_15) .* phi_s1_t4_15) * alpha_s1_t4_15;
    lambda_t0s1_t6s1_15 = diag_post_multiply(lambda_t0s1_t5s1_15, (1 - p_s1_t5_15) .* phi_s1_t5_15) * alpha_s1_t5_15;

    lambda_t0s1_t7s1_15 = diag_post_multiply(lambda_t0s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                          diag_post_multiply(alpha_s1_t6_15, ((1 - psi_1to3_t6_15) .* (1 - psi_1to4_t6_15)));

    lambda_t0s1_t7s3_15 = diag_post_multiply(lambda_t0s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                          diag_post_multiply(alpha_s1_t6_15, psi_1to3_t6_15);
    lambda_t0s1_t7s4_15 = diag_post_multiply(lambda_t0s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                          diag_post_multiply(alpha_s1_t6_15, (1 - psi_1to3_t6_15) .* psi_1to4_t6_15);

    lambda_t0s1_t8s1_15 = diag_post_multiply(lambda_t0s1_t7s1_15, (1 - p_s1_t7_15) .* phi_s1_t7_15) *
                          diag_post_multiply(alpha_s1_t7_15, (1 - psi_1to5_t7_15));
  }

  {
    lambda_t1s1_t4s1_15  = diag_post_multiply(diag_pre_multiply(phi_s1_t1_15, alpha_s1_t1_15) * alpha_s1_t2_15, phi_s1_t3_15) * alpha_s1_t3_15;
    lambda_t1s1_t5s1_15  = diag_post_multiply(lambda_t1s1_t4s1_15, (1 - p_s1_t4_15) .* phi_s1_t4_15) * alpha_s1_t4_15;
    lambda_t1s1_t6s1_15  = diag_post_multiply(lambda_t1s1_t5s1_15, (1 - p_s1_t5_15) .* phi_s1_t5_15) * alpha_s1_t5_15;
    lambda_t1s1_t7s1_15  = diag_post_multiply(lambda_t1s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, ((1 - psi_1to3_t6_15) .* (1 - psi_1to4_t6_15)));
    lambda_t1s1_t7s3_15  = diag_post_multiply(lambda_t1s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, psi_1to3_t6_15);
    lambda_t1s1_t7s4_15  = diag_post_multiply(lambda_t1s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, (1 - psi_1to3_t6_15) .* psi_1to4_t6_15);
    lambda_t1s1_t8s1_15  = diag_post_multiply(lambda_t1s1_t7s1_15, (1 - p_s1_t7_15) .* phi_s1_t7_15) *
                           diag_post_multiply(alpha_s1_t7_15, (1 - psi_1to5_t7_15));
  }

  {
    lambda_t4s1_t5s1_15  = diag_pre_multiply(phi_s1_t4_15, alpha_s1_t4_15);
    lambda_t4s1_t6s1_15  = diag_post_multiply(lambda_t4s1_t5s1_15, (1 - p_s1_t5_15) .* phi_s1_t5_15) * alpha_s1_t5_15;
    lambda_t4s1_t7s1_15  = diag_post_multiply(lambda_t4s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, ((1 - psi_1to3_t6_15) .* (1 - psi_1to4_t6_15)));
    lambda_t4s1_t7s3_15  = diag_post_multiply(lambda_t4s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, psi_1to3_t6_15);
    lambda_t4s1_t7s4_15  = diag_post_multiply(lambda_t4s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, (1 - psi_1to3_t6_15) .* psi_1to4_t6_15);
    lambda_t4s1_t8s1_15  = diag_post_multiply(lambda_t4s1_t7s1_15, (1 - p_s1_t7_15) .* phi_s1_t7_15) *
                           diag_post_multiply(alpha_s1_t7_15, (1 - psi_1to5_t7_15));
    lambda_t4s1_t8s5_15  = diag_post_multiply(lambda_t4s1_t7s1_15, (1 - p_s1_t7_15) .* phi_s1_t7_15) *
                           diag_post_multiply(alpha_s1_t7_15,  psi_1to5_t7_15);
    lambda_t4s1_t9s1_15  = (diag_post_multiply(lambda_t4s1_t8s1_15, (1 - p_s1_t8_15) .* phi_s1_t8_15) * alpha_s1_t8_15) +
                           (diag_post_multiply(lambda_t4s1_t7s3_15, (1 - p_s3_t7_15) .* phi_s3_t8_15) * alpha_s3_t8_15) +
                           (diag_post_multiply(lambda_t4s1_t7s4_15, (1 - p_s4_t7_15) .* phi_s4_t8_15) * alpha_s4_t8_15);
  }

  {
    lambda_t5s1_t6s1_15  = diag_pre_multiply(phi_s1_t5_15, alpha_s1_t5_15);
    lambda_t5s1_t7s1_15  = diag_post_multiply(lambda_t5s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, ((1 - psi_1to3_t6_15) .* (1 - psi_1to4_t6_15)));
    lambda_t5s1_t7s3_15  = diag_post_multiply(lambda_t5s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, psi_1to3_t6_15);
    lambda_t5s1_t7s4_15  = diag_post_multiply(lambda_t5s1_t6s1_15, (1 - p_s1_t6_15) .* phi_s1_t6_15) *
                           diag_post_multiply(alpha_s1_t6_15, (1 - psi_1to3_t6_15) .* psi_1to4_t6_15);
    lambda_t5s1_t8s1_15  = diag_post_multiply(lambda_t5s1_t7s1_15, (1 - p_s1_t7_15) .* phi_s1_t7_15) *
                           diag_post_multiply(alpha_s1_t7_15, (1 - psi_1to5_t7_15));
    lambda_t5s1_t8s5_15  = diag_post_multiply(lambda_t5s1_t7s1_15, (1 - p_s1_t7_15) .* phi_s1_t7_15) *
                           diag_post_multiply(alpha_s1_t7_15,  psi_1to5_t7_15);
    lambda_t5s1_t9s1_15  = (diag_post_multiply(lambda_t5s1_t8s1_15, (1 - p_s1_t8_15) .* phi_s1_t8_15) * alpha_s1_t8_15) +
                           (diag_post_multiply(lambda_t5s1_t7s3_15, (1 - p_s3_t7_15) .* phi_s3_t8_15) * alpha_s3_t8_15) +
                           (diag_post_multiply(lambda_t5s1_t7s4_15, (1 - p_s4_t7_15) .* phi_s4_t8_15) * alpha_s4_t8_15);
  }

  {
    lambda_t6s1_t7s1_15  = diag_pre_multiply(phi_s1_t6_15,
                           diag_post_multiply(alpha_s1_t6_15, ((1 - psi_1to3_t6_15) .* (1 - psi_1to4_t6_15))));
    lambda_t6s1_t7s3_15  = diag_pre_multiply(phi_s1_t6_15,
                           diag_post_multiply(alpha_s1_t6_15, psi_1to3_t6_15));
    lambda_t6s1_t7s4_15  = diag_pre_multiply(phi_s1_t6_15,
                           diag_post_multiply(alpha_s1_t6_15, (1 - psi_1to3_t6_15) .* psi_1to4_t6_15));
  }

  {
    matrix[T, T]      lambda_t7s3_t10s1_15;
    matrix[T, T]      lambda_t7s4_t10s1_15;

    lambda_t7s1_t8s1_15  = diag_pre_multiply(phi_s1_t7_15,
                           diag_post_multiply(alpha_s1_t7_15, (1 - psi_1to5_t7_15)));
    lambda_t7s1_t8s5_15  = diag_pre_multiply(phi_s1_t7_15,
                           diag_post_multiply(alpha_s1_t7_15,  psi_1to5_t7_15));
    lambda_t7s1_t9s1_15  = diag_post_multiply(lambda_t7s1_t8s1_15, (1 - p_s1_t8_15) .* phi_s1_t8_15) * alpha_s1_t8_15;

    lambda_t7s3_t9s1_15  = diag_pre_multiply(phi_s3_t8_15, alpha_s3_t8_15);
    lambda_t7s3_t10s1_15 = diag_post_multiply(lambda_t7s3_t9s1_15, (1 - p_s1_t9_15) .* phi_s1_t9_15) * alpha_s1_t9_15;
    lambda_t7s3_t11s1_15 = (lambda_t7s3_t10s1_15 * chipps2benecia) * alpha_s1_t10_15;

    lambda_t7s4_t9s1_15  = diag_pre_multiply(phi_s4_t8_15, alpha_s4_t8_15);
    lambda_t7s4_t10s1_15 = diag_post_multiply(lambda_t7s4_t9s1_15, (1 - p_s1_t9_15) .* phi_s1_t9_15) * alpha_s1_t9_15;
    lambda_t7s4_t11s1_15 = (lambda_t7s4_t10s1_15 * chipps2benecia) * alpha_s1_t10_15;
  }

  {
    matrix[T, T]      lambda_t8s1_t10s1_15;
    matrix[T, T]      lambda_t8s5_t10s1_15;

    lambda_t8s1_t9s1_15  = diag_pre_multiply(phi_s1_t8_15, alpha_s1_t8_15);
    lambda_t8s1_t10s1_15 = diag_post_multiply(lambda_t8s1_t9s1_15, (1 - p_s1_t9_15) .* phi_s1_t9_15) * alpha_s1_t9_15;
    lambda_t8s1_t11s1_15 = diag_post_multiply(lambda_t8s1_t10s1_15, (1 - p_s1_t10_15) * chipps2benecia) * alpha_s1_t10_15;

    lambda_t8s5_t9s5_15  = diag_pre_multiply(phi_s5_t8_15, alpha_s5_t8_15);
    lambda_t8s5_t10s1_15 = diag_post_multiply(lambda_t8s5_t9s5_15, (1 - p_s5_t9_15) .* phi_s5_t9_15) * alpha_s5_t9_15;
    lambda_t8s5_t11s1_15 = diag_post_multiply(lambda_t8s5_t10s1_15, (1 - p_s1_t10_15) * chipps2benecia) * alpha_s1_t10_15;
  }

  {
    lambda_t9s1_t11s1_15 = diag_post_multiply(diag_pre_multiply(phi_s1_t9_15, alpha_s1_t9_15), (1 - p_s1_t10_15) * chipps2benecia) * alpha_s1_t10_15;
    lambda_t9s5_t11s1_15 = diag_post_multiply(diag_pre_multiply(phi_s5_t9_15, alpha_s5_t9_15), (1 - p_s1_t10_15) * chipps2benecia) * alpha_s1_t10_15;
  }
  
  //2016 lambdas------------------------------------------------------------------------------------------\\
  {
    matrix[2, T]      lambda_t0s1_t2s1_16;
    matrix[2, T]      lambda_t0s1_t3s1_16;
    matrix[2, T]      lambda_t0s1_t8s1_16;
    matrix[2, T]      lambda_t0s1_t8s5_16;

    lambda_t0s1_t1s1_16 = diag_pre_multiply(phi_s1_t0_16, alpha_s1_t0_16);
    lambda_t0s1_t2s1_16 = diag_post_multiply(lambda_t0s1_t1s1_16, (1 - p_s1_t1_16) .* phi_s1_t1_16) * alpha_s1_t1_16;
    lambda_t0s1_t3s1_16 = diag_post_multiply(lambda_t0s1_t2s1_16, (1 - psi_1to2_t2_16)) * alpha_s1_t2_16;

    lambda_t0s1_t4s1_16 = diag_post_multiply(lambda_t0s1_t3s1_16, phi_s1_t3_16) * alpha_s1_t3_16;
    lambda_t0s1_t4s2_16 = diag_post_multiply(lambda_t0s1_t2s1_16, psi_1to2_t2_16 .* phi_s2_t3_16) *  alpha_s2_t3_16;

    lambda_t0s1_t5s1_16 = diag_post_multiply(lambda_t0s1_t4s1_16, (1 - p_s1_t4_16) .* phi_s1_t4_16) * alpha_s1_t4_16;
    lambda_t0s1_t6s1_16 = diag_post_multiply(lambda_t0s1_t5s1_16, (1 - p_s1_t5_16) .* phi_s1_t5_16) * alpha_s1_t5_16;


    lambda_t0s1_t7s1_16 = diag_post_multiply(lambda_t0s1_t6s1_16, (1 - p_s1_t6_16) .* phi_s1_t6_16) *
                          diag_post_multiply(alpha_s1_t6_16, ((1 - psi_1to3_t6_16) .* (1 - psi_1to4_t6_16)));

    lambda_t0s1_t7s3_16 = diag_post_multiply(lambda_t0s1_t6s1_16, (1 - p_s1_t6_16) .* phi_s1_t6_16) *
                          diag_post_multiply(alpha_s1_t6_16, psi_1to3_t6_16);
    lambda_t0s1_t7s4_16 = diag_post_multiply(lambda_t0s1_t6s1_16, (1 - p_s1_t6_16) .* phi_s1_t6_16) *
                          diag_post_multiply(alpha_s1_t6_16, (1 - psi_1to3_t6_16) .*psi_1to4_t6_16);

    lambda_t0s1_t8s1_16 = diag_post_multiply(lambda_t0s1_t7s1_16, (1 - p_s1_t7_16) .* phi_s1_t7_16) *
                       diag_post_multiply(alpha_s1_t7_16, (1 - psi_1to5_t7_16));
    lambda_t0s1_t8s5_16 = diag_post_multiply(lambda_t0s1_t7s1_16, (1 - p_s1_t7_16) .* phi_s1_t7_16) *
                           diag_post_multiply(alpha_s1_t7_16,  psi_1to5_t7_16);

    lambda_t0s1_t9s1_16 = (diag_post_multiply(lambda_t0s1_t8s1_16, phi_s1_t8_16) * alpha_s1_t8_16) +
                          (diag_post_multiply(lambda_t0s1_t4s2_16, (1 - p_s2_t4_16) .* phi_s2_t8_16) * alpha_s2_t8_16) +
                          (diag_post_multiply(lambda_t0s1_t7s3_16, (1 - p_s3_t7_16) .* phi_s3_t8_16) * alpha_s3_t8_16) +
                          (diag_post_multiply(lambda_t0s1_t7s4_16, (1 - p_s4_t7_16) .* phi_s4_t8_16) * alpha_s4_t8_16);

    lambda_t0s1_t9s5_16 = diag_post_multiply(lambda_t0s1_t8s5_16, phi_s5_t8_16) * alpha_s5_t8_16;

    lambda_t0s1_t10s1_16 = (diag_post_multiply(lambda_t0s1_t9s1_16, (1 - p_s1_t9_16) .* phi_s1_t9_16) * alpha_s1_t9_16) +
                           (diag_post_multiply(lambda_t0s1_t9s5_16, (1 - p_s5_t9_16) .* phi_s5_t9_16) * alpha_s5_t9_16);

  }

  {
    matrix[T, T]      lambda_t1s1_t3s1_16;

    lambda_t1s1_t3s1_16  = diag_post_multiply(diag_pre_multiply(phi_s1_t1_16, alpha_s1_t1_16), (1 - psi_1to2_t2_16)) * alpha_s1_t2_16;
    lambda_t1s1_t4s1_16  = diag_post_multiply(lambda_t1s1_t3s1_16, phi_s1_t3_16) * alpha_s1_t3_16;
    lambda_t1s1_t4s2_16  = diag_post_multiply(diag_pre_multiply(phi_s1_t1_16, alpha_s1_t1_16), psi_1to2_t2_16 .* phi_s2_t3_16) *  alpha_s2_t3_16;
    lambda_t1s1_t5s1_16  = diag_post_multiply(lambda_t1s1_t4s1_16, (1 - p_s1_t4_16) .* phi_s1_t4_16) * alpha_s1_t4_16;
    lambda_t1s1_t6s1_16  = diag_post_multiply(lambda_t1s1_t5s1_16, (1 - p_s1_t5_16) .* phi_s1_t5_16) * alpha_s1_t5_16;

    lambda_t1s1_t7s4_16  = diag_post_multiply(lambda_t1s1_t6s1_16, (1 - p_s1_t6_16) .* phi_s1_t6_16) *
                           diag_post_multiply(alpha_s1_t6_16, (1 - psi_1to3_t6_16) .*psi_1to4_t6_16);
  }

  {
    lambda_t4s1_t5s1_16  = diag_pre_multiply(phi_s1_t4_16, alpha_s1_t4_16);
    lambda_t4s1_t6s1_16  = diag_post_multiply(lambda_t4s1_t5s1_16, (1 - p_s1_t5_16) .* phi_s1_t5_16) * alpha_s1_t5_16;
    lambda_t4s1_t7s1_16  = diag_post_multiply(lambda_t4s1_t6s1_16, (1 - p_s1_t6_16) .* phi_s1_t6_16) *
                           diag_post_multiply(alpha_s1_t6_16, ((1 - psi_1to3_t6_16) .* (1 - psi_1to4_t6_16)));

    lambda_t4s2_t9s1_16  = diag_pre_multiply(phi_s2_t8_16,  alpha_s2_t8_16);
  }

  {
    matrix[T, T]      lambda_t5s1_t8s1_16;

    lambda_t5s1_t6s1_16  = diag_pre_multiply(phi_s1_t5_16, alpha_s1_t5_16);
    lambda_t5s1_t7s1_16  = diag_post_multiply(lambda_t5s1_t6s1_16, (1 - p_s1_t6_16) .* phi_s1_t6_16) *
                           diag_post_multiply(alpha_s1_t6_16, ((1 - psi_1to3_t6_16) .* (1 - psi_1to4_t6_16)));
    lambda_t5s1_t7s3_16  = diag_post_multiply(lambda_t5s1_t6s1_16, (1 - p_s1_t6_16) .* phi_s1_t6_16) *
                           diag_post_multiply(alpha_s1_t6_16, psi_1to3_t6_16);
    lambda_t5s1_t7s4_16  = diag_post_multiply(lambda_t5s1_t6s1_16, (1 - p_s1_t6_16) .* phi_s1_t6_16) *
                           diag_post_multiply(alpha_s1_t6_16, (1 - psi_1to3_t6_16) .*psi_1to4_t6_16);
    lambda_t5s1_t8s1_16  = diag_post_multiply(lambda_t5s1_t7s1_16, (1 - p_s1_t7_16) .* phi_s1_t7_16) *
                           diag_post_multiply(alpha_s1_t7_16, (1 - psi_1to5_t7_16));
    lambda_t5s1_t9s1_16  = (diag_post_multiply(lambda_t5s1_t8s1_16, phi_s1_t8_16) * alpha_s1_t8_16) +
                           (diag_post_multiply(lambda_t5s1_t7s3_16, (1 - p_s3_t7_16) .* phi_s3_t8_16) * alpha_s3_t8_16) +
                           (diag_post_multiply(lambda_t5s1_t7s4_16, (1 - p_s4_t7_16) .* phi_s4_t8_16) * alpha_s4_t8_16);
  }

  {
    matrix[T, T]      lambda_t6s1_t8s1_16;

    lambda_t6s1_t7s1_16  = diag_pre_multiply(phi_s1_t6_16,
                           diag_post_multiply(alpha_s1_t6_16, ((1 - psi_1to3_t6_16) .* (1 - psi_1to4_t6_16))));
    lambda_t6s1_t7s3_16  = diag_pre_multiply(phi_s1_t6_16,
                           diag_post_multiply(alpha_s1_t6_16, psi_1to3_t6_16));
    lambda_t6s1_t7s4_16  = diag_pre_multiply(phi_s1_t6_16,
                           diag_post_multiply(alpha_s1_t6_16, (1 - psi_1to3_t6_16) .* psi_1to4_t6_16));
    lambda_t6s1_t8s1_16  = diag_post_multiply(lambda_t6s1_t7s1_16, (1 - p_s1_t7_16) .* phi_s1_t7_16) *
                           diag_post_multiply(alpha_s1_t7_16, (1 - psi_1to5_t7_16));
    lambda_t6s1_t9s1_16  = (diag_post_multiply(lambda_t6s1_t8s1_16, phi_s1_t8_16) * alpha_s1_t8_16) +
                           (diag_post_multiply(lambda_t6s1_t7s3_16, (1 - p_s3_t7_16) .* phi_s3_t8_16) * alpha_s3_t8_16) +
                           (diag_post_multiply(lambda_t6s1_t7s4_16, (1 - p_s4_t7_16) .* phi_s4_t8_16) * alpha_s4_t8_16);
  }

  {
    matrix[T, T]      lambda_t7s1_t8s1_16;
    matrix[T, T]      lambda_t7s1_t8s5_16;

    lambda_t7s1_t8s1_16  = diag_pre_multiply(phi_s1_t7_16, diag_post_multiply(alpha_s1_t7_16, (1 - psi_1to5_t7_16)));
    lambda_t7s1_t8s5_16  = diag_pre_multiply(phi_s1_t7_16, diag_post_multiply(alpha_s1_t7_16,  psi_1to5_t7_16));
    lambda_t7s1_t9s1_16  = diag_post_multiply(lambda_t7s1_t8s1_16, phi_s1_t8_16) * alpha_s1_t8_16;
    lambda_t7s1_t9s5_16  = diag_post_multiply(lambda_t7s1_t8s5_16, phi_s5_t8_16) * alpha_s5_t8_16;
    lambda_t7s1_t10s1_16 = (diag_post_multiply(lambda_t7s1_t9s1_16, (1 - p_s1_t9_16) .* phi_s1_t9_16) * alpha_s1_t9_16) +
                           (diag_post_multiply(lambda_t7s1_t9s5_16, (1 - p_s5_t9_16) .* phi_s5_t9_16) * alpha_s5_t9_16);

    lambda_t7s3_t9s1_16  = diag_pre_multiply(phi_s3_t8_16, alpha_s3_t8_16);
    lambda_t7s3_t10s1_16 = diag_post_multiply(lambda_t7s3_t9s1_16, (1 - p_s1_t9_16) .* phi_s1_t9_16) * alpha_s1_t9_16;

    lambda_t7s4_t9s1_16  = diag_pre_multiply(phi_s4_t8_16, alpha_s4_t8_16);
    lambda_t7s4_t10s1_16 = diag_post_multiply(lambda_t7s4_t9s1_16, (1 - p_s1_t9_16) .* phi_s1_t9_16) * alpha_s1_t9_16;
  }

  {
    lambda_t9s1_t10s1_16 = diag_pre_multiply(phi_s1_t9_16, alpha_s1_t9_16);
    lambda_t9s1_t11s1_16 = diag_post_multiply(lambda_t9s1_t10s1_16, (1 - p_s1_t10_16) * chipps2benecia) * alpha_s1_t10_16;

    lambda_t9s5_t10s1_16 = diag_pre_multiply(phi_s5_t9_16, alpha_s5_t9_16);
    lambda_t9s5_t11s1_16 = diag_post_multiply(lambda_t9s5_t10s1_16, (1 - p_s1_t10_16) * chipps2benecia) * alpha_s1_t10_16;
  }
  
  //2017 lambdas------------------------------------------------------------------------------------------\\
  {
    lambda_t0s1_t1s1_17 = phi_s1_t0_17 * alpha_s1_t0_17;
    lambda_t0s1_t2s1_17 = (lambda_t0s1_t1s1_17 .* ((1 - p_s1_t1_17) .* phi_s1_t1_17)') * alpha_s1_t1_17;
    lambda_t0s1_t3s1_17 = (lambda_t0s1_t2s1_17 .* ((1 - p_s1_t2_17) .* (1 - psi_1to2_t2_17))') * alpha_s1_t2_17;

    lambda_t0s1_t4s1_17 = (lambda_t0s1_t3s1_17 .* ((1 - p_s1_t3_17) .* phi_s1_t3_17)') * alpha_s1_t3_17;
    lambda_t0s1_t4s2_17 = (lambda_t0s1_t2s1_17 .* ((1 - p_s1_t2_17) .* psi_1to2_t2_17 .* phi_s2_t3_17)') *  alpha_s2_t3_17;

    lambda_t0s1_t5s1_17 = (lambda_t0s1_t4s1_17 .* ((1 - p_s1_t4_17) .* phi_s1_t4_17)') * alpha_s1_t4_17;
    lambda_t0s1_t6s1_17 = (lambda_t0s1_t5s1_17 .* ((1 - p_s1_t5_17) .* phi_s1_t5_17)') * alpha_s1_t5_17;


    lambda_t0s1_t7s1_17 = (lambda_t0s1_t6s1_17 .* ((1 - p_s1_t6_17) .* phi_s1_t6_17)') *
                          diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17)));

    lambda_t0s1_t7s3_17 = (lambda_t0s1_t6s1_17 .* ((1 - p_s1_t6_17) .* phi_s1_t6_17)') *
                          diag_post_multiply(alpha_s1_t6_17, psi_1to3_t6_17);
    lambda_t0s1_t7s4_17 = (lambda_t0s1_t6s1_17 .* ((1 - p_s1_t6_17) .* phi_s1_t6_17)') *
                          diag_post_multiply(alpha_s1_t6_17, (1 - psi_1to3_t6_17) .* psi_1to4_t6_17);

    lambda_t0s1_t8s1_17 = (lambda_t0s1_t7s1_17 .* ((1 - p_s1_t7_17) .* phi_s1_t7_17)') *
                       diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17));
    lambda_t0s1_t8s5_17 = (lambda_t0s1_t7s1_17 .* ((1 - p_s1_t7_17) .* phi_s1_t7_17)') *
                       diag_post_multiply(alpha_s1_t7_17,  psi_1to5_t7_17);

    lambda_t0s1_t9s1_17 = ((lambda_t0s1_t8s1_17 .* ((1 - p_s1_t8_17) .* phi_s1_t8_17)') * alpha_s1_t8_17) +
                          ((lambda_t0s1_t4s2_17 .* ((1 - p_s2_t4_17) .* phi_s2_t8_17)') * alpha_s2_t8_17) +
                          ((lambda_t0s1_t7s3_17 .* ((1 - p_s3_t7_17) .* phi_s3_t8_17)') * alpha_s3_t8_17) +
                          ((lambda_t0s1_t7s4_17 .* ((1 - p_s4_t7_17) .* phi_s4_t8_17)') * alpha_s4_t8_17);

    lambda_t0s1_t9s5_17 = (lambda_t0s1_t8s5_17  .* ((1 - p_s5_t8_17) .* phi_s5_t8_17)') * alpha_s5_t8_17;

    lambda_t0s1_t10s1_17 = ((lambda_t0s1_t9s1_17 .* ((1 - p_s1_t9_17) .* phi_s1_t9_17)') * alpha_s1_t9_17) +
                           ((lambda_t0s1_t9s5_17 .* ((1 - p_s5_t9_17) .* phi_s5_t9_17)') * alpha_s5_t9_17);
  }

  {
    lambda_t1s1_t2s1_17  = diag_pre_multiply(phi_s1_t1_17, alpha_s1_t1_17);
    lambda_t1s1_t3s1_17  = diag_post_multiply(lambda_t1s1_t2s1_17, (1 - p_s1_t2_17) .* (1 - psi_1to2_t2_17)) * alpha_s1_t2_17;
    lambda_t1s1_t4s1_17  = diag_post_multiply(lambda_t1s1_t3s1_17, (1 - p_s1_t3_17) .* phi_s1_t3_17) * alpha_s1_t3_17;
    lambda_t1s1_t4s2_17  = diag_post_multiply(lambda_t1s1_t2s1_17, (1 - p_s1_t2_17) .* psi_1to2_t2_17 .* phi_s2_t3_17) *  alpha_s2_t3_17;
    lambda_t1s1_t5s1_17  = diag_post_multiply(lambda_t1s1_t4s1_17, (1 - p_s1_t4_17) .* phi_s1_t4_17) * alpha_s1_t4_17;
    lambda_t1s1_t6s1_17  = diag_post_multiply(lambda_t1s1_t5s1_17, (1 - p_s1_t5_17) .* phi_s1_t5_17) * alpha_s1_t5_17;
    lambda_t1s1_t7s1_17  = diag_post_multiply(lambda_t1s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17)));
    lambda_t1s1_t7s3_17  = diag_post_multiply(lambda_t1s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, psi_1to3_t6_17);
    lambda_t1s1_t7s4_17  = diag_post_multiply(lambda_t1s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, (1 - psi_1to3_t6_17) .* psi_1to4_t6_17);
    lambda_t1s1_t8s1_17  = diag_post_multiply(lambda_t1s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                           diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17));
    lambda_t1s1_t9s1_17  = (diag_post_multiply(lambda_t1s1_t8s1_17, (1 - p_s1_t8_17) .* phi_s1_t8_17) * alpha_s1_t8_17) +
                           (diag_post_multiply(lambda_t1s1_t4s2_17, (1 - p_s2_t4_17) .* phi_s2_t8_17) * alpha_s2_t8_17) +
                           (diag_post_multiply(lambda_t1s1_t7s3_17, (1 - p_s3_t7_17) .* phi_s3_t8_17) * alpha_s3_t8_17) +
                           (diag_post_multiply(lambda_t1s1_t7s4_17, (1 - p_s4_t7_17) .* phi_s4_t8_17) * alpha_s4_t8_17);
  }

  {
    lambda_t2s1_t3s1_17  = diag_pre_multiply((1 - psi_1to2_t2_17), alpha_s1_t2_17);
    lambda_t2s1_t4s2_17  = diag_pre_multiply(psi_1to2_t2_17 .* phi_s2_t3_17,  alpha_s2_t3_17);

    lambda_t2s1_t4s1_17  = diag_post_multiply(lambda_t2s1_t3s1_17, (1 - p_s1_t3_17) .* phi_s1_t3_17) * alpha_s1_t3_17;
    lambda_t2s1_t5s1_17  = diag_post_multiply(lambda_t2s1_t4s1_17, (1 - p_s1_t4_17) .* phi_s1_t4_17) * alpha_s1_t4_17;
    lambda_t2s1_t6s1_17  = diag_post_multiply(lambda_t2s1_t5s1_17, (1 - p_s1_t5_17) .* phi_s1_t5_17) * alpha_s1_t5_17;
    lambda_t2s1_t7s1_17  = diag_post_multiply(lambda_t2s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                          diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17)));
    lambda_t2s1_t7s3_17  = diag_post_multiply(lambda_t2s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                          diag_post_multiply(alpha_s1_t6_17, psi_1to3_t6_17);
    lambda_t2s1_t7s4_17  = diag_post_multiply(lambda_t2s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                          diag_post_multiply(alpha_s1_t6_17, (1 - psi_1to3_t6_17) .* psi_1to4_t6_17);
    lambda_t2s1_t8s1_17  = diag_post_multiply(lambda_t2s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                          diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17));
    lambda_t2s1_t8s5_17  = diag_post_multiply(lambda_t2s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                          diag_post_multiply(alpha_s1_t7_17,  psi_1to5_t7_17);
    lambda_t2s1_t9s1_17  = (diag_post_multiply(lambda_t2s1_t8s1_17, (1 - p_s1_t8_17) .* phi_s1_t8_17) * alpha_s1_t8_17) +
                          (diag_post_multiply(lambda_t2s1_t4s2_17, (1 - p_s2_t4_17) .* phi_s2_t8_17) * alpha_s2_t8_17) +
                          (diag_post_multiply(lambda_t2s1_t7s3_17, (1 - p_s3_t7_17) .* phi_s3_t8_17) * alpha_s3_t8_17) +
                          (diag_post_multiply(lambda_t2s1_t7s4_17, (1 - p_s4_t7_17) .* phi_s4_t8_17) * alpha_s4_t8_17);
    lambda_t2s1_t9s5_17  = diag_post_multiply(lambda_t2s1_t8s5_17, (1 - p_s5_t8_17) .* phi_s5_t8_17) * alpha_s5_t8_17;
    lambda_t2s1_t10s1_17 = (diag_post_multiply(lambda_t2s1_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17) +
                           (diag_post_multiply(lambda_t2s1_t9s5_17, (1 - p_s5_t9_17) .* phi_s5_t9_17) * alpha_s5_t9_17);
  }

  {
    lambda_t3s1_t4s1_17  = diag_pre_multiply(phi_s1_t3_17,  alpha_s1_t3_17);
    lambda_t3s1_t5s1_17  = diag_post_multiply(lambda_t3s1_t4s1_17, (1 - p_s1_t4_17) .* phi_s1_t4_17) * alpha_s1_t4_17;
    lambda_t3s1_t6s1_17  = diag_post_multiply(lambda_t3s1_t5s1_17, (1 - p_s1_t5_17) .* phi_s1_t5_17) * alpha_s1_t5_17;
    lambda_t3s1_t7s1_17  = diag_post_multiply(lambda_t3s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17)));
    lambda_t3s1_t7s3_17  = diag_post_multiply(lambda_t3s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, psi_1to3_t6_17);
    lambda_t3s1_t7s4_17  = diag_post_multiply(lambda_t3s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, (1 - psi_1to3_t6_17) .* psi_1to4_t6_17);
    lambda_t3s1_t8s1_17  = diag_post_multiply(lambda_t3s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                           diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17));
    lambda_t3s1_t8s5_17  = diag_post_multiply(lambda_t3s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                           diag_post_multiply(alpha_s1_t7_17,  psi_1to5_t7_17);
    lambda_t3s1_t9s1_17  = (diag_post_multiply(lambda_t3s1_t8s1_17, (1 - p_s1_t8_17) .* phi_s1_t8_17) * alpha_s1_t8_17) +
                           (diag_post_multiply(lambda_t3s1_t7s3_17, (1 - p_s3_t7_17) .* phi_s3_t8_17) * alpha_s3_t8_17) +
                           (diag_post_multiply(lambda_t3s1_t7s4_17, (1 - p_s4_t7_17) .* phi_s4_t8_17) * alpha_s4_t8_17);
    lambda_t3s1_t9s5_17  = diag_post_multiply(lambda_t3s1_t8s5_17, (1 - p_s5_t8_17) .* phi_s5_t8_17) * alpha_s5_t8_17;
    lambda_t3s1_t10s1_17 = (diag_post_multiply(lambda_t3s1_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17) +
                           (diag_post_multiply(lambda_t3s1_t9s5_17, (1 - p_s5_t9_17) .* phi_s5_t9_17) * alpha_s5_t9_17);
  }

  {
    lambda_t4s1_t5s1_17  = diag_pre_multiply(phi_s1_t4_17, alpha_s1_t4_17);
    lambda_t4s1_t6s1_17  = diag_post_multiply(lambda_t4s1_t5s1_17, (1 - p_s1_t5_17) .* phi_s1_t5_17) * alpha_s1_t5_17;
    lambda_t4s1_t7s1_17  = diag_post_multiply(lambda_t4s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17)));

    lambda_t4s2_t9s1_17  = diag_pre_multiply(phi_s2_t8_17,  alpha_s2_t8_17);
    lambda_t4s2_t10s1_17 = diag_post_multiply(lambda_t4s2_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17;
    lambda_t4s2_t11s1_17 = diag_post_multiply(lambda_t4s2_t10s1_17, (1 - p_s1_t10_17) * chipps2benecia) * alpha_s1_t10_17;
  }

  {
    lambda_t5s1_t6s1_17  = diag_pre_multiply(phi_s1_t5_17, alpha_s1_t5_17);
    lambda_t5s1_t7s1_17  = diag_post_multiply(lambda_t5s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17)));
    lambda_t5s1_t7s3_17  = diag_post_multiply(lambda_t5s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, psi_1to3_t6_17);
    lambda_t5s1_t7s4_17  = diag_post_multiply(lambda_t5s1_t6s1_17, (1 - p_s1_t6_17) .* phi_s1_t6_17) *
                           diag_post_multiply(alpha_s1_t6_17, (1 - psi_1to3_t6_17) .* psi_1to4_t6_17);
    lambda_t5s1_t8s1_17  = diag_post_multiply(lambda_t5s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                           diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17));
    lambda_t5s1_t8s5_17  = diag_post_multiply(lambda_t5s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                           diag_post_multiply(alpha_s1_t7_17,  psi_1to5_t7_17);
    lambda_t5s1_t9s1_17  = (diag_post_multiply(lambda_t5s1_t8s1_17, (1 - p_s1_t8_17) .* phi_s1_t8_17) * alpha_s1_t8_17) +
                           (diag_post_multiply(lambda_t5s1_t7s3_17, (1 - p_s3_t7_17) .* phi_s3_t8_17) * alpha_s3_t8_17) +
                           (diag_post_multiply(lambda_t5s1_t7s4_17, (1 - p_s4_t7_17) .* phi_s4_t8_17) * alpha_s4_t8_17);
    lambda_t5s1_t9s5_17  = diag_post_multiply(lambda_t5s1_t8s5_17, (1 - p_s5_t8_17) .* phi_s5_t8_17) * alpha_s5_t8_17;
    lambda_t5s1_t10s1_17 = (diag_post_multiply(lambda_t5s1_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17) +
                           (diag_post_multiply(lambda_t5s1_t9s5_17, (1 - p_s5_t9_17) .* phi_s5_t9_17) * alpha_s5_t9_17);
  }

  {
    lambda_t6s1_t7s1_17  = diag_pre_multiply(phi_s1_t6_17,
                           diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17))));
    lambda_t6s1_t7s3_17  = diag_pre_multiply(phi_s1_t6_17,
                           diag_post_multiply(alpha_s1_t6_17, psi_1to3_t6_17));
    lambda_t6s1_t7s4_17  = diag_pre_multiply(phi_s1_t6_17,
                           diag_post_multiply(alpha_s1_t6_17, (1 - psi_1to3_t6_17) .* psi_1to4_t6_17));
    lambda_t6s1_t8s1_17  = diag_post_multiply(lambda_t6s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                           diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17));
    lambda_t6s1_t8s5_17  = diag_post_multiply(lambda_t6s1_t7s1_17, (1 - p_s1_t7_17) .* phi_s1_t7_17) *
                           diag_post_multiply(alpha_s1_t7_17,  psi_1to5_t7_17);
    lambda_t6s1_t9s1_17  = (diag_post_multiply(lambda_t6s1_t8s1_17, (1 - p_s1_t8_17) .* phi_s1_t8_17) * alpha_s1_t8_17) +
                           (diag_post_multiply(lambda_t6s1_t7s3_17, (1 - p_s3_t7_17) .* phi_s3_t8_17) * alpha_s3_t8_17) +
                           (diag_post_multiply(lambda_t6s1_t7s4_17, (1 - p_s4_t7_17) .* phi_s4_t8_17) * alpha_s4_t8_17);
    lambda_t6s1_t9s5_17  = diag_post_multiply(lambda_t6s1_t8s5_17, (1 - p_s5_t8_17) .* phi_s5_t8_17) * alpha_s5_t8_17;
    lambda_t6s1_t10s1_17 = (diag_post_multiply(lambda_t6s1_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17) +
                           (diag_post_multiply(lambda_t6s1_t9s5_17, (1 - p_s5_t9_17) .* phi_s5_t9_17) * alpha_s5_t9_17);
  }

  {
    lambda_t7s1_t8s1_17  = diag_pre_multiply(phi_s1_t7_17,
                           diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17)));
    lambda_t7s1_t8s5_17  = diag_pre_multiply(phi_s1_t7_17,
                           diag_post_multiply(alpha_s1_t7_17,  psi_1to5_t7_17));
    lambda_t7s1_t9s1_17  = diag_post_multiply(lambda_t7s1_t8s1_17, (1 - p_s1_t8_17) .* phi_s1_t8_17) * alpha_s1_t8_17;
    lambda_t7s1_t9s5_17  = diag_post_multiply(lambda_t7s1_t8s5_17, (1 - p_s5_t8_17) .* phi_s5_t8_17) * alpha_s5_t8_17;
    lambda_t7s1_t10s1_17 = (diag_post_multiply(lambda_t7s1_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17) +
                           (diag_post_multiply(lambda_t7s1_t9s5_17, (1 - p_s5_t9_17) .* phi_s5_t9_17) * alpha_s5_t9_17);

    lambda_t7s3_t9s1_17  = diag_pre_multiply(phi_s3_t8_17, alpha_s3_t8_17);
    lambda_t7s3_t10s1_17 = diag_post_multiply(lambda_t7s3_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17;
    lambda_t7s3_t11s1_17 = diag_post_multiply(lambda_t7s3_t10s1_17, (1 - p_s1_t10_17) * chipps2benecia) * alpha_s1_t10_17;

    lambda_t7s4_t9s1_17  = diag_pre_multiply(phi_s4_t8_17, alpha_s4_t8_17);
    lambda_t7s4_t10s1_17 = diag_post_multiply(lambda_t7s4_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17;
    lambda_t7s4_t11s1_17 = diag_post_multiply(lambda_t7s4_t10s1_17, (1 - p_s1_t10_17) * chipps2benecia) * alpha_s1_t10_17;
  }

  {
    lambda_t8s1_t9s1_17  = diag_pre_multiply(phi_s1_t8_17, alpha_s1_t8_17);
    lambda_t8s1_t10s1_17 = diag_post_multiply(lambda_t8s1_t9s1_17, (1 - p_s1_t9_17) .* phi_s1_t9_17) * alpha_s1_t9_17;
  }

  {
    lambda_t9s1_t10s1_17 = diag_pre_multiply(phi_s1_t9_17, alpha_s1_t9_17);
    lambda_t9s1_t11s1_17 = diag_post_multiply(lambda_t9s1_t10s1_17, (1 - p_s1_t10_17) * chipps2benecia) * alpha_s1_t10_17;

    lambda_t9s5_t10s1_17 = diag_pre_multiply(phi_s5_t9_17, alpha_s5_t9_17);
  }
  
  //2018 lambdas ---------------------------------------------------------------------------------\\
  {
    lambda_t0s1_t1s1_18 = diag_pre_multiply(phi_s1_t0_18, alpha_s1_t0_18);
  }

  {
    lambda_t1s1_t2s1_18  = diag_pre_multiply(phi_s1_t1_18, alpha_s1_t1_18);
    lambda_t1s1_t3s1_18  = diag_post_multiply(lambda_t1s1_t2s1_18, (1 - p_s1_t2_18)) * alpha_s1_t2_18;
    lambda_t1s1_t4s1_18  = diag_post_multiply(lambda_t1s1_t3s1_18, (1 - p_s1_t3_18) .* phi_s1_t3_18) * alpha_s1_t3_18;
    lambda_t1s1_t5s1_18  = diag_post_multiply(lambda_t1s1_t4s1_18, (1 - p_s1_t4_18) .* phi_s1_t4_18) * alpha_s1_t4_18;
    lambda_t1s1_t6s1_18  = diag_post_multiply(lambda_t1s1_t5s1_18, (1 - p_s1_t5_18) .* phi_s1_t5_18) * alpha_s1_t5_18;
    lambda_t1s1_t7s1_18  = diag_post_multiply(lambda_t1s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, ((1 - psi_1to3_t6_18) .* (1 - psi_1to4_t6_18)));
  }

  {
    matrix[T, T]      lambda_t2s1_t4s2_18;

    lambda_t2s1_t3s1_18  = diag_pre_multiply((1 - psi_1to2_t2_18), alpha_s1_t2_18);
    lambda_t2s1_t4s2_18  = diag_pre_multiply(psi_1to2_t2_18 .* phi_s2_t3_18,  alpha_s2_t3_18);
    
    lambda_t2s1_t4s1_18  = diag_post_multiply(lambda_t2s1_t3s1_18, (1 - p_s1_t3_18) .* phi_s1_t3_18) * alpha_s1_t3_18;
    lambda_t2s1_t5s1_18  = diag_post_multiply(lambda_t2s1_t4s1_18, (1 - p_s1_t4_18) .* phi_s1_t4_18) * alpha_s1_t4_18;
    lambda_t2s1_t6s1_18  = diag_post_multiply(lambda_t2s1_t5s1_18, (1 - p_s1_t5_18) .* phi_s1_t5_18) * alpha_s1_t5_18;
    lambda_t2s1_t7s1_18  = diag_post_multiply(lambda_t2s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                          diag_post_multiply(alpha_s1_t6_18, ((1 - psi_1to3_t6_18) .* (1 - psi_1to4_t6_18)));
    lambda_t2s1_t7s3_18  = diag_post_multiply(lambda_t2s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                          diag_post_multiply(alpha_s1_t6_18, psi_1to3_t6_18);
    lambda_t2s1_t7s4_18  = diag_post_multiply(lambda_t2s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                          diag_post_multiply(alpha_s1_t6_18, (1 - psi_1to3_t6_18) .* psi_1to4_t6_18);
    lambda_t2s1_t8s1_18  = diag_post_multiply(lambda_t2s1_t7s1_18, (1 - p_s1_t7_18) .* phi_s1_t7_18) *
                          diag_post_multiply(alpha_s1_t7_18, (1 - psi_1to5_t7_18));
                  
    lambda_t2s1_t9s1_18  = (diag_post_multiply(lambda_t2s1_t8s1_18, (1 - p_s1_t8_18) .* phi_s1_t8_18) * alpha_s1_t8_18) +
                           (diag_post_multiply(lambda_t2s1_t4s2_18, (1 - p_s2_t4_18) .* phi_s2_t8_18) * alpha_s2_t8_18) +
                           (diag_post_multiply(lambda_t2s1_t7s3_18, (1 - p_s3_t7_18) .* phi_s3_t8_18) * alpha_s3_t8_18) +
                           (diag_post_multiply(lambda_t2s1_t7s4_18, (1 - p_s4_t7_18) .* phi_s4_t8_18) * alpha_s4_t8_18);
  }

  {
    lambda_t3s1_t4s1_18  = diag_pre_multiply(phi_s1_t3_18,  alpha_s1_t3_18);
    lambda_t3s1_t5s1_18  = diag_post_multiply(lambda_t3s1_t4s1_18, (1 - p_s1_t4_18) .* phi_s1_t4_18) * alpha_s1_t4_18;
    lambda_t3s1_t6s1_18  = diag_post_multiply(lambda_t3s1_t5s1_18, (1 - p_s1_t5_18) .* phi_s1_t5_18) * alpha_s1_t5_18;
    lambda_t3s1_t7s1_18  = diag_post_multiply(lambda_t3s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, ((1 - psi_1to3_t6_18) .* (1 - psi_1to4_t6_18)));
    lambda_t3s1_t7s3_18  = diag_post_multiply(lambda_t3s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, psi_1to3_t6_18);
    lambda_t3s1_t7s4_18  = diag_post_multiply(lambda_t3s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, (1 - psi_1to3_t6_18) .* psi_1to4_t6_18);
    lambda_t3s1_t8s1_18  = diag_post_multiply(lambda_t3s1_t7s1_18, (1 - p_s1_t7_18) .* phi_s1_t7_18) *
                           diag_post_multiply(alpha_s1_t7_18, (1 - psi_1to5_t7_18));
    lambda_t3s1_t9s1_18  = (diag_post_multiply(lambda_t3s1_t8s1_18, (1 - p_s1_t8_18) .* phi_s1_t8_18) * alpha_s1_t8_18) +
                           (diag_post_multiply(lambda_t3s1_t7s3_18, (1 - p_s3_t7_18) .* phi_s3_t8_18) * alpha_s3_t8_18) +
                           (diag_post_multiply(lambda_t3s1_t7s4_18, (1 - p_s4_t7_18) .* phi_s4_t8_18) * alpha_s4_t8_18);
  }

  {
    lambda_t4s1_t5s1_18  = diag_pre_multiply(phi_s1_t4_18, alpha_s1_t4_18);
    lambda_t4s1_t6s1_18  = diag_post_multiply(lambda_t4s1_t5s1_18, (1 - p_s1_t5_18) .* phi_s1_t5_18) * alpha_s1_t5_18;
    lambda_t4s1_t7s1_18  = diag_post_multiply(lambda_t4s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, ((1 - psi_1to3_t6_18) .* (1 - psi_1to4_t6_18)));
    lambda_t4s1_t7s3_18  = diag_post_multiply(lambda_t4s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, psi_1to3_t6_18);
    lambda_t4s1_t7s4_18  = diag_post_multiply(lambda_t4s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, (1 - psi_1to3_t6_18) .* psi_1to4_t6_18);
  }

  {
    lambda_t5s1_t6s1_18  = diag_pre_multiply(phi_s1_t5_18, alpha_s1_t5_18);
    lambda_t5s1_t7s1_18  = diag_post_multiply(lambda_t5s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, ((1 - psi_1to3_t6_18) .* (1 - psi_1to4_t6_18)));
    lambda_t5s1_t7s3_18  = diag_post_multiply(lambda_t5s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, psi_1to3_t6_18);
    lambda_t5s1_t7s4_18  = diag_post_multiply(lambda_t5s1_t6s1_18, (1 - p_s1_t6_18) .* phi_s1_t6_18) *
                           diag_post_multiply(alpha_s1_t6_18, (1 - psi_1to3_t6_18) .* psi_1to4_t6_18);
  }

  {
    lambda_t6s1_t7s1_18  = diag_pre_multiply(phi_s1_t6_18,
                           diag_post_multiply(alpha_s1_t6_18, ((1 - psi_1to3_t6_18) .* (1 - psi_1to4_t6_18))));
    lambda_t6s1_t7s3_18  = diag_pre_multiply(phi_s1_t6_18,
                           diag_post_multiply(alpha_s1_t6_18, psi_1to3_t6_18));
    lambda_t6s1_t7s4_18  = diag_pre_multiply(phi_s1_t6_18,
                           diag_post_multiply(alpha_s1_t6_18, (1 - psi_1to3_t6_18) .* psi_1to4_t6_18));
    lambda_t6s1_t8s5_18  = diag_post_multiply(lambda_t6s1_t7s1_18, (1 - p_s1_t7_18) .* phi_s1_t7_18) *
                           diag_post_multiply(alpha_s1_t7_18,  psi_1to5_t7_18);
  }

  {
    lambda_t7s1_t8s1_18  = diag_pre_multiply(phi_s1_t7_18,
                           diag_post_multiply(alpha_s1_t7_18, (1 - psi_1to5_t7_18)));
    lambda_t7s1_t8s5_18  = diag_pre_multiply(phi_s1_t7_18,
                           diag_post_multiply(alpha_s1_t7_18,  psi_1to5_t7_18));
    lambda_t7s1_t9s1_18  = diag_post_multiply(lambda_t7s1_t8s1_18, (1 - p_s1_t8_18) .* phi_s1_t8_18) * alpha_s1_t8_18;
    lambda_t7s1_t9s5_18  = diag_post_multiply(lambda_t7s1_t8s5_18, (1 - p_s5_t8_18) .* phi_s5_t8_18) * alpha_s5_t8_18;
    lambda_t7s1_t10s1_18 = (diag_post_multiply(lambda_t7s1_t9s1_18, (1 - p_s1_t9_18) .* phi_s1_t9_18) * alpha_s1_t9_18) +
                           (diag_post_multiply(lambda_t7s1_t9s5_18, (1 - p_s5_t9_18) .* phi_s5_t9_18) * alpha_s5_t9_18);
    lambda_t7s3_t9s1_18  = diag_pre_multiply(phi_s3_t8_18, alpha_s3_t8_18);
    lambda_t7s3_t10s1_18 = diag_post_multiply(lambda_t7s3_t9s1_18, (1 - p_s1_t9_18) .* phi_s1_t9_18) * alpha_s1_t9_18;
    lambda_t7s4_t9s1_18  = diag_pre_multiply(phi_s4_t8_18, alpha_s4_t8_18);
    lambda_t7s4_t10s1_18 = diag_post_multiply(lambda_t7s4_t9s1_18, (1 - p_s1_t9_18) .* phi_s1_t9_18) * alpha_s1_t9_18;
  }

  {
    lambda_t8s1_t9s1_18  = diag_pre_multiply(phi_s1_t8_18, alpha_s1_t8_18);
    lambda_t8s1_t10s1_18 = diag_post_multiply(lambda_t8s1_t9s1_18, (1 - p_s1_t9_18) .* phi_s1_t9_18) * alpha_s1_t9_18;

    lambda_t8s5_t9s5_18  = diag_pre_multiply(phi_s5_t8_18, alpha_s5_t8_18);
  }

  {
    lambda_t9s1_t10s1_18 = diag_pre_multiply(phi_s1_t9_18, alpha_s1_t9_18);
    lambda_t9s1_t11s1_18 = diag_post_multiply(lambda_t9s1_t10s1_18, (1 - p_s1_t10_18) * chipps2benecia) * alpha_s1_t10_18;

    lambda_t9s5_t10s1_18 = diag_pre_multiply(phi_s5_t9_18, alpha_s5_t9_18);
  }
  //2014 chis----\\
  {
    vector[T]             chi_10_s1_14; //prob last seen chipps
    vector[T]             chi_3_s1_14; //prob last seen blw fremont weir
    vector[T]             chi_2_s1_14; //prob last seen abv fremont weir

    for (t in 1:T)
      chi_10_s1_14[t] =  1 - chipps2benecia * (alpha_s1_t10_14[t,t:T] * (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]));

    for (t in 1:T){
      chi_9_s1_14[t]  = (1 - phi_s1_t9_14[t]) +
                        phi_s1_t9_14[t] * (alpha_s1_t9_14[t,t:T] *
                        ( (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                        chi_10_s1_14[t:T] +
                          (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));
      chi_9_s5_14[t]  = (1 - phi_s5_t9_14[t]) +
                        phi_s5_t9_14[t] * (alpha_s5_t9_14[t, t:T] *
                        ( (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                        chi_10_s1_14[t:T] +
                        (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));
    }

    for (t in 1:T){
      chi_8_s1_14[t]  = (1 - phi_s1_t8_14[t]) +
                        phi_s1_t8_14[t] * (alpha_s1_t8_14[t,t:T] *
                        ( (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                        (1 - p_s1_t9_14[t:T]) .* chi_9_s1_14[t:T] +
                          (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));
      chi_8_s5_14[t]  = (1 - phi_s5_t8_14[t]) +
                        phi_s5_t8_14[t] * (alpha_s5_t8_14[t,t:T] *
                        ( (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                        (1 - p_s5_t9_14[t:T]) .* chi_9_s5_14[t:T] +
                        (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));
    }

    for (t in 1:T){
      chi_7_s1_14[t]  = (1 - phi_s1_t7_14[t]) + phi_s1_t7_14[t] * (alpha_s1_t7_14[t, t:T] *
                          ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                          ((1 - psi_1to5_t7_14[t:T]) .* (1 - p_s1_t8_14[t:T]) .* chi_8_s1_14[t:T] +
                           psi_1to5_t7_14[t:T] .* (1 - p_s5_t8_14[t:T]) .* chi_8_s5_14[t:T]) +
                          (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));

      chi_7_s3_14[t]  = (1 - phi_s3_t8_14[t]) + phi_s3_t8_14[t] * (alpha_s3_t8_14[t, t:T] *
                        ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                        (1 - p_s1_t9_14[t:T]) .* chi_9_s1_14[t:T] +
                        (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));

      chi_7_s4_14[t]  = (1 - phi_s4_t8_14[t]) + phi_s4_t8_14[t] * (alpha_s4_t8_14[t, t:T] *
                        ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                        (1 - p_s1_t9_14[t:T]) .* chi_9_s1_14[t:T] +
                        (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));
    }

    for (t in 1:T)
      chi_6_s1_14[t]  = (1 - phi_s1_t6_14[t]) + phi_s1_t6_14[t] * (alpha_s1_t6_14[t, t:T] *
                    ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                    (((1 - psi_1to3_t6_14[t:T]) .* (1 - psi_1to4_t6_14[t:T])) .* (1 - p_s1_t7_14[t:T]) .* chi_7_s1_14[t:T] +
                     psi_1to3_t6_14[t:T] .* (1 - p_s3_t7_14[t:T]) .* chi_7_s3_14[t:T] +
                    ((1 - psi_1to3_t6_14[t:T]) .* psi_1to4_t6_14[t:T]) .* (1 - p_s4_t7_14[t:T]) .* chi_7_s4_14[t:T]) +
                    (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))) );

    for (t in 1:T)
      chi_5_s1_14[t]  = (1 - phi_s1_t5_14[t]) + phi_s1_t5_14[t] * (alpha_s1_t5_14[t, t:T] *
                        ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                        (1 - p_s1_t6_14[t:T]) .* chi_6_s1_14[t:T] +
                        (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));
    for (t in 1:T)
      chi_4_s1_14[t]  = (1 - phi_s1_t4_14[t]) + phi_s1_t4_14[t] * (alpha_s1_t4_14[t, t:T] *
                          ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                          (1 - p_s1_t5_14[t:T]) .* chi_5_s1_14[t:T] +
                          (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));

    for (t in 1:T)
      chi_3_s1_14[t]  = (1 - phi_s1_t3_14[t]) + phi_s1_t3_14[t] * (alpha_s1_t3_14[t, t:T] *
                          ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                          (1 - p_s1_t4_14[t:T]) .* chi_4_s1_14[t:T] +
                          (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));

    for (t in 1:T)
      chi_2_s1_14[t]  = alpha_s1_t2_14[t, t:T] *
                        ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                         chi_3_s1_14[t:T] +
                        (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t])));
    for (t in 1:T)
      chi_1_s1_14[t]  = (1 - phi_s1_t1_14[t]) + phi_s1_t1_14[t] .* (alpha_s1_t1_14[t, t:T] *
                          ((tag_surv_prob_14[t:T]/tag_surv_prob_14[t]) .*
                          chi_2_s1_14[t:T] +
                          (1 - (tag_surv_prob_14[t:T]/tag_surv_prob_14[t]))));

    chi_0_s1_14 = (1 - phi_s1_t0_14) + phi_s1_t0_14 * (alpha_s1_t0_14 *
                       (tag_surv_prob_14 .*
                       (1 - p_s1_t1_14) .* chi_1_s1_14 +
                       (1 - tag_surv_prob_14)));
  }
  //2015 chis----\\
  {
    vector[T]             chi_10_s1_15; //prob last seen chipps
    vector[T]             chi_3_s1_15; //prob last seen blw fremont weir
    vector[T]             chi_2_s1_15; //prob last seen abv fremont weir

    for (t in 1:T)
      chi_10_s1_15[t] =  1 - chipps2benecia * (alpha_s1_t10_15[t,t:T] * (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]));


    for (t in 1:T){
      chi_9_s1_15[t]  = (1 - phi_s1_t9_15[t]) +
                        phi_s1_t9_15[t] * (alpha_s1_t9_15[t,t:T] *
                        ( (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                        chi_10_s1_15[t:T] +
                          (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));
      chi_9_s5_15[t]  = (1 - phi_s5_t9_15[t]) +
                        phi_s5_t9_15[t] * (alpha_s5_t9_15[t, t:T] *
                        ( (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                        chi_10_s1_15[t:T] +
                        (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));
    }

    for (t in 1:T){
      chi_8_s1_15[t]  = (1 - phi_s1_t8_15[t]) +
                        phi_s1_t8_15[t] * (alpha_s1_t8_15[t,t:T] *
                        ( (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                        (1 - p_s1_t9_15[t:T]) .* chi_9_s1_15[t:T] +
                          (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));
      chi_8_s5_15[t]  = (1 - phi_s5_t8_15[t]) +
                        phi_s5_t8_15[t] * (alpha_s5_t8_15[t,t:T] *
                        ( (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                        (1 - p_s5_t9_15[t:T]) .* chi_9_s5_15[t:T] +
                        (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));

    }
    
    for (t in 1:T){
      chi_7_s1_15[t]  = (1 - phi_s1_t7_15[t]) + phi_s1_t7_15[t] * (alpha_s1_t7_15[t, t:T] *
                          ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                          ((1 - psi_1to5_t7_15[t:T]) .* (1 - p_s1_t8_15[t:T]) .* chi_8_s1_15[t:T] +
                           psi_1to5_t7_15[t:T] .* (1 - p_s5_t8_15[t:T]) .* chi_8_s5_15[t:T]) +
                          (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));

      chi_7_s3_15[t]  = (1 - phi_s3_t8_15[t]) + phi_s3_t8_15[t] * (alpha_s3_t8_15[t, t:T] *
                        ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                        (1 - p_s1_t9_15[t:T]) .* chi_9_s1_15[t:T] +
                        (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));

      chi_7_s4_15[t]  = (1 - phi_s4_t8_15[t]) + phi_s4_t8_15[t] * (alpha_s4_t8_15[t, t:T] *
                        ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                        (1 - p_s1_t9_15[t:T]) .* chi_9_s1_15[t:T] +
                        (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));
    }

    for (t in 1:T)
      chi_6_s1_15[t]  = (1 - phi_s1_t6_15[t]) + phi_s1_t6_15[t] * (alpha_s1_t6_15[t, t:T] *
                    ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                    (((1 - psi_1to3_t6_15[t:T]) .* (1 - psi_1to4_t6_15[t:T])) .* (1 - p_s1_t7_15[t:T]) .* chi_7_s1_15[t:T] +
                     psi_1to3_t6_15[t:T] .* (1 - p_s3_t7_15[t:T]) .* chi_7_s3_15[t:T] +
                    ((1 - psi_1to3_t6_15[t:T]) .* psi_1to4_t6_15[t:T]) .* (1 - p_s4_t7_15[t:T]) .* chi_7_s4_15[t:T]) +
                    (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))) );

    for (t in 1:T)
      chi_5_s1_15[t]  = (1 - phi_s1_t5_15[t]) + phi_s1_t5_15[t] * (alpha_s1_t5_15[t, t:T] *
                        ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                        (1 - p_s1_t6_15[t:T]) .* chi_6_s1_15[t:T] +
                        (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));
                        

    for (t in 1:T)
      chi_4_s1_15[t]  = (1 - phi_s1_t4_15[t]) + phi_s1_t4_15[t] * (alpha_s1_t4_15[t, t:T] *
                          ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                          (1 - p_s1_t5_15[t:T]) .* chi_5_s1_15[t:T] +
                          (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));
                          

    for (t in 1:T)
      chi_3_s1_15[t]  = (1 - phi_s1_t3_15[t]) + phi_s1_t3_15[t] * (alpha_s1_t3_15[t, t:T] *
                          ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                          (1 - p_s1_t4_15[t:T]) .* chi_4_s1_15[t:T] +
                          (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));

    for (t in 1:T)
      chi_2_s1_15[t]  = alpha_s1_t2_15[t, t:T] *
                        ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                         chi_3_s1_15[t:T] +
                        (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t])));
                        

    for (t in 1:T)
      chi_1_s1_15[t]  = (1 - phi_s1_t1_15[t]) + phi_s1_t1_15[t] .* (alpha_s1_t1_15[t, t:T] *
                          ((tag_surv_prob_15[t:T]/tag_surv_prob_15[t]) .*
                          chi_2_s1_15[t:T] +
                          (1 - (tag_surv_prob_15[t:T]/tag_surv_prob_15[t]))));
                          

    for (g in 1:2)
      chi_0_s1_15[g] = (1 - phi_s1_t0_15[g]) + phi_s1_t0_15[g] * (alpha_s1_t0_15[g,] *
                       (tag_surv_prob_15 .*
                       (1 - p_s1_t1_15) .* chi_1_s1_15 +
                       (1 - tag_surv_prob_15)));
                       

  }
  //2016 chis----\\
  {
    vector[T]             chi_8_s1_16; //prob last seen sac blw georgiana
    vector[T]             chi_8_s5_16; //prob last seen abv fremont weir
    vector[T]             chi_3_s1_16; //prob last seen blw fremont weir
    vector[T]             chi_2_s1_16; //prob last seen abv fremont weir

    for (t in 1:T)
      chi_10_s1_16[t] =  1 - chipps2benecia * (alpha_s1_t10_16[t,t:T] * (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]));

    for (t in 1:T){
      chi_9_s1_16[t]  = (1 - phi_s1_t9_16[t]) +
                        phi_s1_t9_16[t] * (alpha_s1_t9_16[t,t:T] *
                        ( (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                        (1 - p_s1_t10_16[t:T]) .* chi_10_s1_16[t:T] +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));
      chi_9_s5_16[t]  = (1 - phi_s5_t9_16[t]) +
                        phi_s5_t9_16[t] * (alpha_s5_t9_16[t, t:T] *
                        ( (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                        (1 - p_s1_t10_16[t:T]) .* chi_10_s1_16[t:T] +
                        (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));
    }

    for (t in 1:T){
      chi_8_s1_16[t]  = (1 - phi_s1_t8_16[t]) +
                        phi_s1_t8_16[t] * (alpha_s1_t8_16[t,t:T] *
                        ( (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                        (1 - p_s1_t9_16[t:T]) .* chi_9_s1_16[t:T] +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));
      chi_8_s5_16[t]  = (1 - phi_s5_t8_16[t]) +
                        phi_s5_t8_16[t] * (alpha_s5_t8_16[t,t:T] *
                        ( (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                        (1 - p_s5_t9_16[t:T]) .* chi_9_s5_16[t:T] +
                        (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));

    }

    for (t in 1:T){
      chi_7_s1_16[t]  = (1 - phi_s1_t7_16[t]) + phi_s1_t7_16[t] * (alpha_s1_t7_16[t, t:T] *
                          ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                          ((1 - psi_1to5_t7_16[t:T]) .* chi_8_s1_16[t:T] +
                           psi_1to5_t7_16[t:T] .* chi_8_s5_16[t:T]) +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));

      chi_7_s3_16[t]  = (1 - phi_s3_t8_16[t]) + phi_s3_t8_16[t] * (alpha_s3_t8_16[t, t:T] *
                        ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                        (1 - p_s1_t9_16[t:T]) .* chi_9_s1_16[t:T] +
                        (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));

      chi_7_s4_16[t]  = (1 - phi_s4_t8_16[t]) + phi_s4_t8_16[t] * (alpha_s4_t8_16[t, t:T] *
                        ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                        (1 - p_s1_t9_16[t:T]) .* chi_9_s1_16[t:T] +
                        (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));
    }

    for (t in 1:T)
      chi_6_s1_16[t]  = (1 - phi_s1_t6_16[t]) + phi_s1_t6_16[t] * (alpha_s1_t6_16[t, t:T] *
                    ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                    (((1 - psi_1to3_t6_16[t:T]) .* (1 - psi_1to4_t6_16[t:T])) .* (1 - p_s1_t7_16[t:T]) .* chi_7_s1_16[t:T] +
                     psi_1to3_t6_16[t:T] .* (1 - p_s3_t7_16[t:T]) .* chi_7_s3_16[t:T] +
                    ((1 - psi_1to3_t6_16[t:T]) .* psi_1to4_t6_16[t:T]) .* (1 - p_s4_t7_16[t:T]) .* chi_7_s4_16[t:T]) +
                    (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))) );

    for (t in 1:T)
      chi_5_s1_16[t]  = (1 - phi_s1_t5_16[t]) + phi_s1_t5_16[t] * (alpha_s1_t5_16[t, t:T] *
                        ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                        (1 - p_s1_t6_16[t:T]) .* chi_6_s1_16[t:T] +
                        (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));
    for (t in 1:T){
      chi_4_s1_16[t]  = (1 - phi_s1_t4_16[t]) + phi_s1_t4_16[t] * (alpha_s1_t4_16[t, t:T] *
                          ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                          (1 - p_s1_t5_16[t:T]) .* chi_5_s1_16[t:T] +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));

      chi_4_s2_16[t]  = (1 - phi_s2_t8_16[t]) + phi_s2_t8_16[t] * (alpha_s2_t8_16[t, t:T] *
                          ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                          (1 - p_s1_t9_16[t:T]) .* chi_9_s1_16[t:T] +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));
    }

    for (t in 1:T)
      chi_3_s1_16[t]  = (1 - phi_s1_t3_16[t]) + phi_s1_t3_16[t] * (alpha_s1_t3_16[t, t:T] *
                          ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                          (1 - p_s1_t4_16[t:T]) .* chi_4_s1_16[t:T] +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));

    //Fish enter Yolo right after Abv Fremont Weir phi_s1_t2 = 1 for those that remain in Sacramento R.
    for (t in 1:T)
      chi_2_s1_16[t]  = (1 - psi_1to2_t2_16[t])  * (alpha_s1_t2_16[t, t:T] *
                          ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                          chi_3_s1_16[t:T] +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t])))) +
                        psi_1to2_t2_16[t] * ((1 - phi_s2_t3_16[t]) + phi_s2_t3_16[t] * (alpha_s2_t3_16[t, t:T] *
                          ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                          (1 - p_s2_t4_16[t:T]) .* chi_4_s2_16[t:T] +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t])))));
    for (t in 1:T)
      chi_1_s1_16[t]  = (1 - phi_s1_t1_16[t]) + phi_s1_t1_16[t] .* (alpha_s1_t1_16[t, t:T] *
                          ((tag_surv_prob_16[t:T]/tag_surv_prob_16[t]) .*
                          chi_2_s1_16[t:T] +
                          (1 - (tag_surv_prob_16[t:T]/tag_surv_prob_16[t]))));
    for (g in 1:2)
      chi_0_s1_16[g] = (1 - phi_s1_t0_16[g]) + phi_s1_t0_16[g] * (alpha_s1_t0_16[g,] *
                    (tag_surv_prob_16 .*
                     (1 - p_s1_t1_16) .* chi_1_s1_16 +
                     (1 - tag_surv_prob_16)));
  }
  //2017 chis----\\
  {
    for (t in 1:T)
      chi_10_s1_17[t] =  1 - chipps2benecia * (alpha_s1_t10_17[t,t:T] * (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]));

    for (t in 1:T){
      chi_9_s1_17[t]  = (1 - phi_s1_t9_17[t]) +
                        phi_s1_t9_17[t] * (alpha_s1_t9_17[t,t:T] *
                        ( (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                        (1 - p_s1_t10_17[t:T]) .* chi_10_s1_17[t:T] +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));
      chi_9_s5_17[t]  = (1 - phi_s5_t9_17[t]) +
                        phi_s5_t9_17[t] * (alpha_s5_t9_17[t, t:T] *
                        ( (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                        (1 - p_s1_t10_17[t:T]) .* chi_10_s1_17[t:T] +
                        (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));
    }

    for (t in 1:T){
      chi_8_s1_17[t]  = (1 - phi_s1_t8_17[t]) +
                        phi_s1_t8_17[t] * (alpha_s1_t8_17[t,t:T] *
                        ( (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                        (1 - p_s1_t9_17[t:T]) .* chi_9_s1_17[t:T] +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));
      chi_8_s5_17[t]  = (1 - phi_s5_t8_17[t]) +
                        phi_s5_t8_17[t] * (alpha_s5_t8_17[t,t:T] *
                        ( (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                        (1 - p_s5_t9_17[t:T]) .* chi_9_s5_17[t:T] +
                        (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));

    }

    for (t in 1:T){
      chi_7_s1_17[t]  = (1 - phi_s1_t7_17[t]) + phi_s1_t7_17[t] * (alpha_s1_t7_17[t, t:T] *
                          ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                          ((1 - psi_1to5_t7_17[t:T]) .* (1 - p_s1_t8_17[t:T]) .* chi_8_s1_17[t:T] +
                           psi_1to5_t7_17[t:T] .* (1 - p_s5_t8_17[t:T]) .* chi_8_s5_17[t:T]) +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));

      chi_7_s3_17[t]  = (1 - phi_s3_t8_17[t]) + phi_s3_t8_17[t] * (alpha_s3_t8_17[t, t:T] *
                        ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                        (1 - p_s1_t9_17[t:T]) .* chi_9_s1_17[t:T] +
                        (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));

      chi_7_s4_17[t]  = (1 - phi_s4_t8_17[t]) + phi_s4_t8_17[t] * (alpha_s4_t8_17[t, t:T] *
                        ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                        (1 - p_s1_t9_17[t:T]) .* chi_9_s1_17[t:T] +
                        (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));
    }

    for (t in 1:T)
      chi_6_s1_17[t]  = (1 - phi_s1_t6_17[t]) + phi_s1_t6_17[t] * (alpha_s1_t6_17[t, t:T] *
                    ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                    (((1 - psi_1to3_t6_17[t:T]) .* (1 - psi_1to4_t6_17[t:T])) .* (1 - p_s1_t7_17[t:T]) .* chi_7_s1_17[t:T] +
                     psi_1to3_t6_17[t:T] .* (1 - p_s3_t7_17[t:T]) .* chi_7_s3_17[t:T] +
                    ((1 - psi_1to3_t6_17[t:T]) .* psi_1to4_t6_17[t:T]) .* (1 - p_s4_t7_17[t:T]) .* chi_7_s4_17[t:T]) +
                    (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))) );

    for (t in 1:T)
    chi_5_s1_17[t]  = (1 - phi_s1_t5_17[t]) + phi_s1_t5_17[t] * (alpha_s1_t5_17[t, t:T] *
                      ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                      (1 - p_s1_t6_17[t:T]) .* chi_6_s1_17[t:T] +
                      (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));
    for (t in 1:T){
      chi_4_s1_17[t]  = (1 - phi_s1_t4_17[t]) + phi_s1_t4_17[t] * (alpha_s1_t4_17[t, t:T] *
                          ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                          (1 - p_s1_t5_17[t:T]) .* chi_5_s1_17[t:T] +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));

      chi_4_s2_17[t]  = (1 - phi_s2_t8_17[t]) + phi_s2_t8_17[t] * (alpha_s2_t8_17[t, t:T] *
                          ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                          (1 - p_s1_t9_17[t:T]) .* chi_9_s1_17[t:T] +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));
    }

    for (t in 1:T)
      chi_3_s1_17[t]  = (1 - phi_s1_t3_17[t]) + phi_s1_t3_17[t] * (alpha_s1_t3_17[t, t:T] *
                          ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                          (1 - p_s1_t4_17[t:T]) .* chi_4_s1_17[t:T] +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));

    //Fish enter Yolo right after Abv Fremont Weir phi_s1_t2 = 1 for those that remain in Sacramento R.
    for (t in 1:T)
      chi_2_s1_17[t]  = (1 - psi_1to2_t2_17[t])  * (alpha_s1_t2_17[t, t:T] *
                          ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                          (1 - p_s1_t3_17[t:T]) .* chi_3_s1_17[t:T] +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t])))) +
                        psi_1to2_t2_17[t] * ((1 - phi_s2_t3_17[t]) + phi_s2_t3_17[t] * (alpha_s2_t3_17[t, t:T] *
                          ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                          (1 - p_s2_t4_17[t:T]) .* chi_4_s2_17[t:T] +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t])))));
    for (t in 1:T)
      chi_1_s1_17[t]  = (1 - phi_s1_t1_17[t]) + phi_s1_t1_17[t] .* (alpha_s1_t1_17[t, t:T] *
                          ((tag_surv_prob_17[t:T]/tag_surv_prob_17[t]) .*
                          (1 - p_s1_t2_17[t:T]) .* chi_2_s1_17[t:T] +
                          (1 - (tag_surv_prob_17[t:T]/tag_surv_prob_17[t]))));

    chi_0_s1_17 = (1 - phi_s1_t0_17) + phi_s1_t0_17 * (alpha_s1_t0_17 *
                        (tag_surv_prob_17 .*
                        (1 - p_s1_t1_17) .* chi_1_s1_17 +
                        (1 - tag_surv_prob_17)));
  }
  //2018 chis----\\
  {
    for (t in 1:T){
      chi_10_s1_18[1, t] =  1 - chipps2benecia * (alpha_s1_t10_18[t,t:T] * (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]));
      chi_10_s1_18[2, t] =  1 - chipps2benecia * (alpha_s1_t10_18[t,t:T] * (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]));
    }

    for (t in 1:T){
      chi_9_s1_18[1, t]  = (1 - phi_s1_t9_18[t]) +
                        phi_s1_t9_18[t] * (alpha_s1_t9_18[t,t:T] *
                        ( (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                        (1 - p_s1_t10_18[t:T]) .* chi_10_s1_18[1, t:T] +
                          (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
                          
                          
      chi_9_s5_18[1, t]  = (1 - phi_s5_t9_18[t]) +
                        phi_s5_t9_18[t] * (alpha_s5_t9_18[t, t:T] *
                        ( (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                        (1 - p_s1_t10_18[t:T]) .* chi_10_s1_18[1, t:T] +
                        (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
                        
      chi_9_s1_18[2, t]  = (1 - phi_s1_t9_18[t]) +
                        phi_s1_t9_18[t] * (alpha_s1_t9_18[t,t:T] *
                        ( (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                        (1 - p_s1_t10_18[t:T]) .* chi_10_s1_18[2, t:T] +
                          (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));
      chi_9_s5_18[2, t]  = (1 - phi_s5_t9_18[t]) +
                        phi_s5_t9_18[t] * (alpha_s5_t9_18[t, t:T] *
                        ( (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                        (1 - p_s1_t10_18[t:T]) .* chi_10_s1_18[2, t:T] +
                        (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));                  
    }

    
    for (t in 1:T){
      chi_8_s1_18[1, t]  = (1 - phi_s1_t8_18[t]) +
                        phi_s1_t8_18[t] * (alpha_s1_t8_18[t,t:T] *
                        ( (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                        (1 - p_s1_t9_18[t:T]) .* chi_9_s1_18[1, t:T] +
                          (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1,t]))));
      chi_8_s5_18[1, t]  = (1 - phi_s5_t8_18[t]) +
                        phi_s5_t8_18[t] * (alpha_s5_t8_18[t,t:T] *
                        ( (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                        (1 - p_s5_t9_18[t:T]) .* chi_9_s5_18[1, t:T] +
                        (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
      chi_8_s1_18[2, t]  = (1 - phi_s1_t8_18[t]) +
                        phi_s1_t8_18[t] * (alpha_s1_t8_18[t,t:T] *
                        ( (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                        (1 - p_s1_t9_18[t:T]) .* chi_9_s1_18[2, t:T] +
                          (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2,t]))));
      chi_8_s5_18[2, t]  = (1 - phi_s5_t8_18[t]) +
                        phi_s5_t8_18[t] * (alpha_s5_t8_18[t,t:T] *
                        ( (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                        (1 - p_s5_t9_18[t:T]) .* chi_9_s5_18[2, t:T] +
                        (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));    

    }
    
    for (t in 1:T){
      chi_7_s1_18[1, t]  = (1 - phi_s1_t7_18[t]) + phi_s1_t7_18[t] * (alpha_s1_t7_18[t, t:T] *
                          ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                          ((1 - psi_1to5_t7_18[t:T]) .* (1 - p_s1_t8_18[t:T]) .* chi_8_s1_18[1, t:T] +
                           psi_1to5_t7_18[t:T] .* (1 - p_s5_t8_18[t:T]) .* chi_8_s5_18[1, t:T]) +
                          (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));

      chi_7_s3_18[1, t]  = (1 - phi_s3_t8_18[t]) + phi_s3_t8_18[t] * (alpha_s3_t8_18[t, t:T] *
                        ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                        (1 - p_s1_t9_18[t:T]) .* chi_9_s1_18[1, t:T] +
                        (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));

      chi_7_s4_18[1, t]  = (1 - phi_s4_t8_18[t]) + phi_s4_t8_18[t] * (alpha_s4_t8_18[t, t:T] *
                        ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                        (1 - p_s1_t9_18[t:T]) .* chi_9_s1_18[1, t:T] +
                        (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
                        
      chi_7_s1_18[2, t]  = (1 - phi_s1_t7_18[t]) + phi_s1_t7_18[t] * (alpha_s1_t7_18[t, t:T] *
                          ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                          ((1 - psi_1to5_t7_18[t:T]) .* (1 - p_s1_t8_18[t:T]) .* chi_8_s1_18[2, t:T] +
                           psi_1to5_t7_18[t:T] .* (1 - p_s5_t8_18[t:T]) .* chi_8_s5_18[2, t:T]) +
                          (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));

      chi_7_s3_18[2, t]  = (1 - phi_s3_t8_18[t]) + phi_s3_t8_18[t] * (alpha_s3_t8_18[t, t:T] *
                        ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                        (1 - p_s1_t9_18[t:T]) .* chi_9_s1_18[2, t:T] +
                        (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));

      chi_7_s4_18[2, t]  = (1 - phi_s4_t8_18[t]) + phi_s4_t8_18[t] * (alpha_s4_t8_18[t, t:T] *
                        ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                        (1 - p_s1_t9_18[t:T]) .* chi_9_s1_18[2, t:T] +
                        (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));                        
    }
    for (t in 1:T){
      chi_6_s1_18[1, t]  = (1 - phi_s1_t6_18[t]) + phi_s1_t6_18[t] * (alpha_s1_t6_18[t, t:T] *
                              ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                              (((1 - psi_1to3_t6_18[t:T]) .* (1 - psi_1to4_t6_18[t:T])) .* (1 - p_s1_t7_18[t:T]) .* chi_7_s1_18[1, t:T] +
                               psi_1to3_t6_18[t:T] .* (1 - p_s3_t7_18[t:T]) .* chi_7_s3_18[1, t:T] +
                              ((1 - psi_1to3_t6_18[t:T]) .* psi_1to4_t6_18[t:T]) .* (1 - p_s4_t7_18[t:T]) .* chi_7_s4_18[1, t:T]) +
                              (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))) );
              
      chi_6_s1_18[2, t]  = (1 - phi_s1_t6_18[t]) + phi_s1_t6_18[t] * (alpha_s1_t6_18[t, t:T] *
                            ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                            (((1 - psi_1to3_t6_18[t:T]) .* (1 - psi_1to4_t6_18[t:T])) .* (1 - p_s1_t7_18[t:T]) .* chi_7_s1_18[2, t:T] +
                             psi_1to3_t6_18[t:T] .* (1 - p_s3_t7_18[t:T]) .* chi_7_s3_18[2, t:T] +
                            ((1 - psi_1to3_t6_18[t:T]) .* psi_1to4_t6_18[t:T]) .* (1 - p_s4_t7_18[t:T]) .* chi_7_s4_18[2, t:T]) +
                            (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))) );              
    }

   for (t in 1:T){
      chi_5_s1_18[1, t]  = (1 - phi_s1_t5_18[t]) + phi_s1_t5_18[t] * (alpha_s1_t5_18[t, t:T] *
                      ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                      (1 - p_s1_t6_18[t:T]) .* chi_6_s1_18[1, t:T] +
                      (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
                      
      chi_5_s1_18[2, t]  = (1 - phi_s1_t5_18[t]) + phi_s1_t5_18[t] * (alpha_s1_t5_18[t, t:T] *
                      ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                      (1 - p_s1_t6_18[t:T]) .* chi_6_s1_18[2, t:T] +
                      (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));
   }

    for (t in 1:T){
      chi_4_s1_18[1, t]  = (1 - phi_s1_t4_18[t]) + phi_s1_t4_18[t] * (alpha_s1_t4_18[t, t:T] *
                          ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                          (1 - p_s1_t5_18[t:T]) .* chi_5_s1_18[1, t:T] +
                          (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
                          
      chi_4_s1_18[2, t]  = (1 - phi_s1_t4_18[t]) + phi_s1_t4_18[t] * (alpha_s1_t4_18[t, t:T] *
                          ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                          (1 - p_s1_t5_18[t:T]) .* chi_5_s1_18[2, t:T] +
                          (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));
                          
      chi_4_s2_18[1, t]  = (1 - phi_s2_t8_18[t]) + phi_s2_t8_18[t] * (alpha_s2_t8_18[t, t:T] *
                          ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                          (1 - p_s1_t9_18[t:T]) .* chi_9_s1_18[1, t:T] +
                          (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
                          
      chi_4_s2_18[2, t]  = (1 - phi_s2_t8_18[t]) + phi_s2_t8_18[t] * (alpha_s2_t8_18[t, t:T] *
                          ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                          (1 - p_s1_t9_18[t:T]) .* chi_9_s1_18[1, t:T] +
                          (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));                                 
    }
    

    for (t in 1:T){
      chi_3_s1_18[1, t]  = (1 - phi_s1_t3_18[t]) + phi_s1_t3_18[t] * (alpha_s1_t3_18[t, t:T] *
                    ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                    (1 - p_s1_t4_18[t:T]) .* chi_4_s1_18[1, t:T] +
                    (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
                    
      chi_3_s1_18[2, t]  = (1 - phi_s1_t3_18[t]) + phi_s1_t3_18[t] * (alpha_s1_t3_18[t, t:T] *
                    ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                    (1 - p_s1_t4_18[t:T]) .* chi_4_s1_18[2, t:T] +
                    (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]))));
    }

    //Fish enter Yolo right after Abv Fremont Weir phi_s1_t2 = 1 for those that remain in Sacramento R.

    for (t in 1:T){
      chi_2_s1_18[1, t]  = ((1 - psi_1to2_t2_18[t]) * (alpha_s1_t2_18[t, t:T] *
                          ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                          (1 - p_s1_t3_18[t:T]) .* chi_3_s1_18[1, t:T] +
                          (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))))) +
                          psi_1to2_t2_18[t] * ((1 - phi_s2_t3_18[t]) + phi_s2_t3_18[t] * (alpha_s2_t3_18[t, t:T] *
                          ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                          (1 - p_s2_t4_18[t:T]) .* chi_4_s2_18[1, t:T] +
                          (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t])))));
                          
      chi_2_s1_18[2, t]  = ((1 - psi_1to2_t2_18[t]) * alpha_s1_t2_18[t, t:T] *
                          ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                          (1 - p_s1_t3_18[t:T]) .* chi_3_s1_18[2, t:T] +
                          (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t])))) +
                          psi_1to2_t2_18[t] * ((1 - phi_s2_t3_18[t]) + phi_s2_t3_18[t] * (alpha_s2_t3_18[t, t:T] *
                          ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                          (1 - p_s2_t4_18[t:T]) .* chi_4_s2_18[2, t:T] +
                          (1 - (tag_surv_prob_18[2,t:T]/tag_surv_prob_18[2, t])))))
                          ;  
    }


    for (t in 1:T){
      chi_1_s1_18[1, t]  = (1 - phi_s1_t1_18[t]) + phi_s1_t1_18[t] .* (alpha_s1_t1_18[t, t:T] *
                      ((tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]) .*
                      (1 - p_s1_t2_18[t:T]) .* chi_2_s1_18[1, t:T] +
                      (1 - (tag_surv_prob_18[1, t:T]/tag_surv_prob_18[1, t]))));
                      
      chi_1_s1_18[2, t]  = (1 - phi_s1_t1_18[t]) + phi_s1_t1_18[t] .* (alpha_s1_t1_18[t, t:T] *
                           ((tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t]) .*
                           (1 - p_s1_t2_18[t:T]) .* chi_2_s1_18[2, t:T] +
                           (1 - (tag_surv_prob_18[2, t:T]/tag_surv_prob_18[2, t])))); 
    }


    chi_0_s1_18[1] = (1 - phi_s1_t0_18[1]) + phi_s1_t0_18[1] * (alpha_s1_t0_18[1,] *
                        (tag_surv_prob_18[1] .*
                        (1 - p_s1_t1_18) .* chi_1_s1_18[1,] +
                        (1 - tag_surv_prob_18[1])));
                        
    chi_0_s1_18[2] = (1 - phi_s1_t0_18[2]) + phi_s1_t0_18[2] * (alpha_s1_t0_18[2,] *
                        (tag_surv_prob_18[2,] .*
                        (1 - p_s1_t1_18) .* chi_1_s1_18[2,] +
                        (1 - tag_surv_prob_18[2,])));
                        

  }

}

model {
    // Define priors ----- \\
    //Survival Priors ----\\

    sigma_eps_s1_t0 ~ normal(0, 0.5); 

    beta_phi_s1_t0_14    ~ student_t(7, 0, 2);
    beta_phi_s1_t0_15    ~ student_t(7, 0, 2);
    beta_phi_s1_t0_16    ~ student_t(7, 0, 2);
    beta_phi_s1_t0_17    ~ student_t(7, 0, 2);
    beta_phi_s1_t0_18    ~ student_t(7, 0, 2);

    beta_phi_s1_t1[1] ~ student_t(7, 0, 2);
    beta_phi_s1_t3[1] ~ student_t(7, 0, 2);
    beta_phi_s1_t4[1] ~ student_t(7, 0, 2);
    beta_phi_s1_t5[1] ~ student_t(7, 0, 2);
    beta_phi_s1_t6[1] ~ student_t(7, 0, 2);
    beta_phi_s1_t7[1] ~ student_t(7, 0, 2);
    beta_phi_s1_t8[1] ~ student_t(7, 0, 2);
    beta_phi_s1_t9[1] ~ student_t(7, 0, 2);

    beta_phi_s2_t3[1] ~ student_t(7, 0, 2);
    beta_phi_s2_t8[1] ~ student_t(7, 0, 2);

    beta_phi_s3_t8[1] ~ student_t(7, 0, 2);
    beta_phi_s4_t8[1] ~ student_t(7, 0, 2);
    beta_phi_s5_t8[1] ~ student_t(7, 0, 2);
    beta_phi_s5_t9[1] ~ student_t(7, 0, 2);

    beta_phi_s1_t1[2:3] ~ student_t(7, 0, 1);
    beta_phi_s1_t3[2:3] ~ student_t(7, 0, 1);
    beta_phi_s1_t4[2:3] ~ student_t(7, 0, 1);
    beta_phi_s1_t5[2:3] ~ student_t(7, 0, 1);
    beta_phi_s1_t6[2:3] ~ student_t(7, 0, 1);
    beta_phi_s1_t7[2:3] ~ student_t(7, 0, 1);
    beta_phi_s1_t8[2:3] ~ student_t(7, 0, 1);
    beta_phi_s1_t9[2:3] ~ student_t(7, 0, 1);

    beta_phi_s2_t3[2:3] ~ student_t(7, 0, 1);
    beta_phi_s2_t8[2:3] ~ student_t(7, 0, 1);

    beta_phi_s3_t8[2:3] ~ student_t(7, 0, 1);
    beta_phi_s4_t8[2:3] ~ student_t(7, 0, 1);
    beta_phi_s5_t8[2:3] ~ student_t(7, 0, 1);
    beta_phi_s5_t9[2:4] ~ student_t(7, 0, 1);

    //Travel Time Priors ----\\
    mu_v_s1_t0_14    ~ normal(4, 2);
    mu_v_s1_t0_15    ~ normal(4, 2);
    mu_v_s1_t0_16    ~ normal(4, 2);
    mu_v_s1_t0_17    ~ normal(4, 2);
    mu_v_s1_t0_18    ~ normal(4, 2);

    beta_mu_s1_t1[1] ~ normal(0, 0.5);
    beta_mu_s1_t2[1] ~ normal(0, 0.5);
    beta_mu_s1_t2b[1] ~ normal(0, 1);
    beta_mu_s1_t3[1] ~ normal(0, 0.5);
    beta_mu_s1_t4[1] ~ normal(1, 1);
    beta_mu_s1_t5[1] ~ normal(1, 1);
    beta_mu_s1_t6[1] ~ normal(1, 1);
    beta_mu_s1_t7[1] ~ normal(1, 1);
    beta_mu_s1_t8[1] ~ normal(1, 1);
    beta_mu_s1_t9[1] ~ normal(1, 1);
    beta_mu_s1_t10   ~ normal(1, 1);

    beta_mu_s2_t3[1] ~ normal(1, 1);
    beta_mu_s2_t8[1] ~ normal(1, 1);

    beta_mu_s3_t8[1] ~ normal(1, 1);
    beta_mu_s4_t8[1] ~ normal(1, 1);
    beta_mu_s5_t8[1] ~ normal(1, 1);
    beta_mu_s5_t9[1] ~ normal(2, 1);

    beta_mu_s1_t1[2] ~ normal(0, 1);
    beta_mu_s1_t2b[2] ~ normal(-1, 1);
    beta_mu_s1_t2[2] ~ normal(0, 1);
    beta_mu_s1_t3[2] ~ normal(0, 1);
    beta_mu_s1_t4[2] ~ normal(0, 1);
    beta_mu_s1_t5[2] ~ normal(0, 1);
    beta_mu_s1_t6[2] ~ normal(0, 1);
    beta_mu_s1_t7[2] ~ normal(0, 1);
    beta_mu_s1_t8[2] ~ normal(0, 1);
    beta_mu_s1_t9[2] ~ normal(0, 1);

    beta_mu_s2_t3[2] ~ normal(0, 1);
    beta_mu_s2_t8[2] ~ normal(0, 1);

    beta_mu_s3_t8[2] ~ normal(0, 1);
    beta_mu_s4_t8[2] ~ normal(0, 1);
    beta_mu_s5_t8[2] ~ normal(0, 1);
    beta_mu_s5_t9[2] ~ normal(0, 1);

    tt_eps_s1_t0_14 ~ normal(0, 1);
    tt_eps_s1_t0_15 ~ normal(0, 1);
    tt_eps_s1_t0_16 ~ normal(0, 1);
    tt_eps_s1_t0_17 ~ normal(0, 1);
    tt_eps_s1_t0_18 ~ normal(0, 1);    
    
    sigma_v_s1_t0_14 ~ normal(0, 1);
    sigma_v_s1_t0_15 ~ normal(0, 1);
    sigma_v_s1_t0_16 ~ normal(0, 1);
    sigma_v_s1_t0_17 ~ normal(0, 1);
    sigma_v_s1_t0_18 ~ normal(0, 1);
    
    sigma_tt_s1_t1 ~ normal(0, 1);
    sigma_tt_s1_t2b ~ normal(0, 1);
    sigma_tt_s1_t2 ~ normal(0, 1);
    sigma_tt_s1_t3 ~ normal(0, 1);
    sigma_tt_s1_t4 ~ normal(0, 1);
    sigma_tt_s1_t5 ~ normal(0, 1);
    sigma_tt_s1_t6 ~ normal(0, 1);
    sigma_tt_s1_t7 ~ normal(0, 1);
    sigma_tt_s1_t8 ~ normal(0, 1);
    sigma_tt_s1_t9 ~ normal(0, 1);    
    
    sigma_tt_s2_t3 ~ normal(0, 1);
    sigma_tt_s2_t8 ~ normal(0, 1);
    sigma_tt_s3_t8 ~ normal(0, 1);
    sigma_tt_s4_t8 ~ normal(0, 1);
    sigma_tt_s5_t8 ~ normal(0, 1);
    sigma_tt_s5_t9 ~ normal(0, 1);
    //Routing Priors ----\\
    beta_psi_1to2_t2[1]  ~ student_t(7, 0, 2);
    beta_psi_1to3_t6[1]  ~ student_t(7, 0, 2);
    beta_psi_1to4_t6[1]  ~ student_t(7, 0, 2);
    beta_psi_1to5_t7[1]  ~ student_t(7, 0, 2);
    beta_psi_1to2_t2[2]  ~ student_t(7, 0, 1);
    beta_psi_1to3_t6[2]  ~ student_t(7, 0, 1);
    beta_psi_1to4_t6[2]  ~ student_t(7, 0, 1);
    beta_psi_1to5_t7[2]  ~ student_t(7, 0, 1);
    //Detection 2014 Priors ----\\
    beta_p_s1_t1_14[1] ~ student_t(7, 0, 2);
    beta_p_s1_t4_14[1] ~ student_t(7, 0, 2);
    beta_p_s1_t5_14[1] ~ student_t(7, 0, 2);
    beta_p_s1_t6_14[1] ~ student_t(7, 0, 2);
    beta_p_s1_t7_14[1] ~ student_t(7, 0, 2);
    beta_p_s1_t8_14[1] ~ student_t(7, 0, 2);
    beta_p_s1_t9_14[1] ~ student_t(7, 0, 2);

    beta_p_s3_t7_14[1] ~ student_t(7, 0, 2);
    beta_p_s4_t7_14[1] ~ student_t(7, 0, 2);
    beta_p_s5_t8_14[1] ~ student_t(7, 0, 2);
    beta_p_s5_t9_14[1] ~ student_t(7, 0, 2);

    beta_p_s1_t1_14[2] ~ student_t(7, 0, 1);
    beta_p_s1_t4_14[2] ~ student_t(7, 0, 1);
    beta_p_s1_t5_14[2] ~ student_t(7, 0, 1);
    beta_p_s1_t6_14[2] ~ student_t(7, 0, 1);
    beta_p_s1_t7_14[2] ~ student_t(7, 0, 1);
    beta_p_s1_t8_14[2] ~ student_t(7, 0, 1);
    beta_p_s1_t9_14[2] ~ student_t(7, 0, 1);

    beta_p_s3_t7_14[2] ~ student_t(7, 0, 1);
    beta_p_s4_t7_14[2] ~ student_t(7, 0, 1);
    beta_p_s5_t8_14[2] ~ student_t(7, 0, 1);
    beta_p_s5_t9_14[2] ~ student_t(7, 0, 1);
    
    //Detection 2015 Priors ----\\

    beta_p_s1_t1_15[1] ~ student_t(7, 0, 2);
    beta_p_s1_t4_15[1] ~ student_t(7, 0, 2);
    beta_p_s1_t5_15[1] ~ student_t(7, 0, 2);
    beta_p_s1_t6_15[1] ~ student_t(7, 0, 2);
    beta_p_s1_t7_15[1] ~ student_t(7, 0, 2);
    beta_p_s1_t8_15[1] ~ student_t(7, 0, 2);
    beta_p_s1_t9_15[1] ~ student_t(7, 0, 2);

    beta_p_s3_t7_15[1] ~ student_t(7, 0, 2);
    beta_p_s4_t7_15[1] ~ student_t(7, 0, 2);
    beta_p_s5_t8_15[1] ~ student_t(7, 0, 2);
    beta_p_s5_t9_15[1] ~ student_t(7, 0, 2);

    beta_p_s1_t1_15[2] ~ student_t(7, 0, 1);
    beta_p_s1_t4_15[2] ~ student_t(7, 0, 1);
    beta_p_s1_t5_15[2] ~ student_t(7, 0, 1);
    beta_p_s1_t6_15[2] ~ student_t(7, 0, 1);
    beta_p_s1_t7_15[2] ~ student_t(7, 0, 1);
    beta_p_s1_t8_15[2] ~ student_t(7, 0, 1);
    beta_p_s1_t9_15[2] ~ student_t(7, 0, 1);

    beta_p_s3_t7_15[2] ~ student_t(7, 0, 1);
    beta_p_s4_t7_15[2] ~ student_t(7, 0, 1);
    beta_p_s5_t8_15[2] ~ student_t(7, 0, 1);
    beta_p_s5_t9_15[2] ~ student_t(7, 0, 1);
    
    //Detection 2016 Priors ----\\
    beta_p_s1_t1_16[1] ~ student_t(7, 0, 2);
    beta_p_s1_t4_16[1] ~ student_t(7, 0, 2);
    beta_p_s1_t5_16[1] ~ student_t(7, 0, 2);
    beta_p_s1_t6_16[1] ~ student_t(7, 0, 2);
    beta_p_s1_t7_16[1] ~ student_t(7, 0, 2);
    beta_p_s1_t9_16[1] ~ student_t(7, 0, 2);
    beta_p_s1_t10_16[1] ~ student_t(7, 0, 2);

    beta_p_s2_t4_16[1] ~ student_t(7, 0, 2);

    beta_p_s3_t7_16[1] ~ student_t(7, 0, 2);
    beta_p_s4_t7_16[1] ~ student_t(7, 0, 2);
    beta_p_s5_t9_16[1] ~ student_t(7, 0, 2);

    beta_p_s1_t1_16[2] ~ student_t(7, 0, 1);
    beta_p_s1_t4_16[2] ~ student_t(7, 0, 1);
    beta_p_s1_t5_16[2] ~ student_t(7, 0, 1);
    beta_p_s1_t6_16[2] ~ student_t(7, 0, 1);
    beta_p_s1_t7_16[2] ~ student_t(7, 0, 1);
    beta_p_s1_t9_16[2] ~ student_t(7, 0, 1);
    beta_p_s1_t10_16[2] ~ student_t(7, 0, 1);

    beta_p_s2_t4_16[2] ~ student_t(7, 0, 1);

    beta_p_s3_t7_16[2] ~ student_t(7, 0, 1);
    beta_p_s4_t7_16[2] ~ student_t(7, 0, 1);
    beta_p_s5_t9_16[2] ~ student_t(7, 0, 1);
    
    //Detection 2017 Priors ----\\

    beta_p_s1_t1_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t2_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t3_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t4_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t5_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t6_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t7_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t8_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t9_17[1] ~ student_t(7, 0, 2);
    beta_p_s1_t10_17[1] ~ student_t(7, 0, 2);

    beta_p_s2_t4_17[1] ~ student_t(7, 0, 2);

    beta_p_s3_t7_17[1] ~ student_t(7, 0, 2);
    beta_p_s4_t7_17[1] ~ student_t(7, 0, 2);
    beta_p_s5_t9_17[1] ~ student_t(7, 0, 2);

    beta_p_s1_t1_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t2_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t3_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t4_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t5_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t6_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t7_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t8_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t9_17[2] ~ student_t(7, 0, 1);
    beta_p_s1_t10_17[2] ~ student_t(7, 0, 1);

    beta_p_s2_t4_17[2] ~ student_t(7, 0, 1);

    beta_p_s3_t7_17[2] ~ student_t(7, 0, 1);
    beta_p_s4_t7_17[2] ~ student_t(7, 0, 1);
    beta_p_s5_t9_17[2] ~ student_t(7, 0, 1);
    
    //Detection 2018 Priors ----\\

    beta_p_s1_t1_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t2_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t3_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t4_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t5_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t6_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t7_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t8_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t9_18[1] ~ student_t(7, 0, 2);
    beta_p_s1_t10_18[1] ~ student_t(7, 0, 2);

    beta_p_s2_t4_18[1] ~ student_t(7, 0, 2);
    
    beta_p_s3_t7_18[1] ~ student_t(7, 0, 2);
    beta_p_s4_t7_18[1] ~ student_t(7, 0, 2);
    beta_p_s5_t8_18[1] ~ student_t(7, 0, 2);
    beta_p_s5_t9_18[1] ~ student_t(7, 0, 2);

    beta_p_s1_t1_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t2_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t3_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t4_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t5_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t6_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t7_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t8_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t9_18[2] ~ student_t(7, 0, 1);
    beta_p_s1_t10_18[2] ~ student_t(7, 0, 1);
    beta_p_s2_t4_18[2] ~ student_t(7, 0, 1);

    beta_p_s3_t7_18[2] ~ student_t(7, 0, 1);
    beta_p_s4_t7_18[2] ~ student_t(7, 0, 1);
    beta_p_s5_t8_18[2] ~ student_t(7, 0, 1);
    beta_p_s5_t9_18[2] ~ student_t(7, 0, 1);
  
  
    // Likelihood ---- \\
    target += log(chi_0_s1_14) * l_0_s1_14;


    for (t in 1:max_T){
      target += (log(lambda_t0s1_t1s1_14[t]) + log(p_s1_t1_14[t]) + log(tag_surv_prob_14[t])) * m_t0s1_t1s1_14[t];
      target += (log(lambda_t0s1_t4s1_14[t]) + log(p_s1_t4_14[t]) + log(tag_surv_prob_14[t])) * m_t0s1_t4s1_14[t];
      target += (log(lambda_t0s1_t5s1_14[t]) + log(p_s1_t5_14[t]) + log(tag_surv_prob_14[t])) * m_t0s1_t5s1_14[t];
    }

    for (g in 1:2){
      target += log(chi_0_s1_15[g]) * l_0_s1_15[g];

      for (t in 1:max_T){
        target += (log(lambda_t0s1_t1s1_15[g, t]) + log(p_s1_t1_15[t]) + log(tag_surv_prob_15[t])) * m_t0s1_t1s1_15[g, t];

        target += (log(lambda_t0s1_t4s1_15[g, t]) + log(p_s1_t4_15[t]) + log(tag_surv_prob_15[t])) * m_t0s1_t4s1_15[g, t];

        target += (log(lambda_t0s1_t5s1_15[g, t]) + log(p_s1_t5_15[t]) + log(tag_surv_prob_15[t])) * m_t0s1_t5s1_15[g, t];

        target += (log(lambda_t0s1_t7s1_15[g, t]) + log(p_s1_t7_15[t]) + log(tag_surv_prob_15[t])) * m_t0s1_t7s1_15[g, t];

        target += (log(lambda_t0s1_t7s3_15[g, t]) + log(p_s3_t7_15[t]) + log(tag_surv_prob_15[t])) * m_t0s1_t7s3_15[g, t];

        target += (log(lambda_t0s1_t7s4_15[g, t]) + log(p_s4_t7_15[t]) + log(tag_surv_prob_15[t])) * m_t0s1_t7s4_15[g, t];

      }
    }

    for (g in 1:2){
      target += log(chi_0_s1_16[g]) * l_0_s1_16[g];

      for (t in 1:max_T){
        target += (log(lambda_t0s1_t1s1_16[g, t]) + log(p_s1_t1_16[t]) + log(tag_surv_prob_16[t])) * m_t0s1_t1s1_16[g, t];
        if (m_t0s1_t4s2_16[g, t] > 0)
         target += (log(lambda_t0s1_t4s2_16[g, t]) + log(p_s2_t4_16[t]) + log(tag_surv_prob_16[t])) * m_t0s1_t4s2_16[g, t];
        target += (log(lambda_t0s1_t4s1_16[g, t]) + log(p_s1_t4_16[t]) + log(tag_surv_prob_16[t])) * m_t0s1_t4s1_16[g, t];
        target += (log(lambda_t0s1_t5s1_16[g, t]) + log(p_s1_t5_16[t]) + log(tag_surv_prob_16[t])) * m_t0s1_t5s1_16[g, t];
        target += (log(lambda_t0s1_t7s1_16[g, t]) + log(p_s1_t7_16[t]) + log(tag_surv_prob_16[t])) * m_t0s1_t7s1_16[g, t];
        target += (log(lambda_t0s1_t10s1_16[g, t]) + log(p_s1_t10_16[t]) + log(tag_surv_prob_16[t])) * m_t0s1_t10s1_16[g, t];
      }
    }

    target += log(chi_0_s1_17) * l_0_s1_17;

    for (t in 1:max_T){
      target += (log(lambda_t0s1_t1s1_17[t]) + log(p_s1_t1_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t1s1_17[t];
      target += (log(lambda_t0s1_t3s1_17[t]) + log(p_s1_t3_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t3s1_17[t];
      if (t >= 54){
        target += (log(lambda_t0s1_t4s1_17[t]) + log(p_s1_t4_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t4s1_17[t];
      }
      target += (log(lambda_t0s1_t5s1_17[t]) + log(p_s1_t5_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t5s1_17[t];
      target += (log(lambda_t0s1_t6s1_17[t]) + log(p_s1_t6_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t6s1_17[t];
      target += (log(lambda_t0s1_t7s1_17[t]) + log(p_s1_t7_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t7s1_17[t];
      if (t <= 40){
        target += (log(lambda_t0s1_t8s1_17[t]) + log(p_s1_t8_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t8s1_17[t];
      }
      target += (log(lambda_t0s1_t9s1_17[t]) + log(p_s1_t9_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t9s1_17[t];
      target += (log(lambda_t0s1_t10s1_17[t]) + log(p_s1_t10_17[t]) + log(tag_surv_prob_17[t])) * m_t0s1_t10s1_17[t];
    }
  

  for (t in 1:max_T){
    target += log(chi_1_s1_14[t]) * l_1_s1_14[t];
    target += log(chi_4_s1_14[t]) * l_4_s1_14[t];
    target += log(chi_5_s1_14[t]) * l_5_s1_14[t];
    target += log(chi_6_s1_14[t]) * l_6_s1_14[t];
    target += log(chi_7_s1_14[t]) * l_7_s1_14[t];
    target += log(chi_7_s3_14[t]) * l_7_s3_14[t];
    target += log(chi_7_s4_14[t]) * l_7_s4_14[t];
    target += log(chi_8_s1_14[t]) * l_8_s1_14[t];
    target += log(chi_8_s5_14[t]) * l_8_s5_14[t];
    target += log(chi_9_s1_14[t]) * l_9_s1_14[t];
    target += log(chi_9_s5_14[t]) * l_9_s5_14[t];

    target += log(chi_1_s1_15[t]) * l_1_s1_15[t];
    target += log(chi_4_s1_15[t]) * l_4_s1_15[t];
    target += log(chi_5_s1_15[t]) * l_5_s1_15[t];
    target += log(chi_6_s1_15[t]) * l_6_s1_15[t];
    target += log(chi_7_s1_15[t]) * l_7_s1_15[t];
    target += log(chi_7_s3_15[t]) * l_7_s3_15[t];
    target += log(chi_7_s4_15[t]) * l_7_s4_15[t];
    target += log(chi_8_s1_15[t]) * l_8_s1_15[t];
    target += log(chi_8_s5_15[t]) * l_8_s5_15[t];
    target += log(chi_9_s1_15[t]) * l_9_s1_15[t];
    target += log(chi_9_s5_15[t]) * l_9_s5_15[t];

    target += log(chi_1_s1_16[t]) * l_1_s1_16[t];
    target += log(chi_4_s2_16[t]) * l_4_s2_16[t];
    target += log(chi_4_s1_16[t]) * l_4_s1_16[t];
    target += log(chi_5_s1_16[t]) * l_5_s1_16[t];
    target += log(chi_6_s1_16[t]) * l_6_s1_16[t];
    target += log(chi_7_s1_16[t]) * l_7_s1_16[t];
    target += log(chi_7_s3_16[t]) * l_7_s3_16[t];
    target += log(chi_7_s4_16[t]) * l_7_s4_16[t];
    target += log(chi_9_s1_16[t]) * l_9_s1_16[t];
    target += log(chi_9_s5_16[t]) * l_9_s5_16[t];
    target += log(chi_10_s1_16[t]) * l_10_s1_16[t];

    target += log(chi_1_s1_17[t]) * l_1_s1_17[t];
    target += log(chi_2_s1_17[t]) * l_2_s1_17[t];
    target += log(chi_3_s1_17[t]) * l_3_s1_17[t];
    target += log(chi_4_s2_17[t]) * l_4_s2_17[t];
    target += log(chi_4_s1_17[t]) * l_4_s1_17[t];
    target += log(chi_5_s1_17[t]) * l_5_s1_17[t];
    target += log(chi_6_s1_17[t]) * l_6_s1_17[t];
    target += log(chi_7_s1_17[t]) * l_7_s1_17[t];
    target += log(chi_7_s3_17[t]) * l_7_s3_17[t];
    target += log(chi_7_s4_17[t]) * l_7_s4_17[t];
    target += log(chi_8_s1_17[t]) * l_8_s1_17[t];
    target += log(chi_8_s5_17[t]) * l_8_s5_17[t];
    target += log(chi_9_s1_17[t]) * l_9_s1_17[t];
    target += log(chi_9_s5_17[t]) * l_9_s5_17[t];
    target += log(chi_10_s1_17[t]) * l_10_s1_17[t];
  }


  for (g in 1:2){
    target += log(chi_0_s1_18[g]) * l_0_s1_18[g];

    for (t in 1:max_T){
      target += (log(lambda_t0s1_t1s1_18[g, t]) + log(p_s1_t1_18[t]) + log(tag_surv_prob_18[g, t])) * m_t0s1_t1s1_18[g, t];
      
      target += log(chi_1_s1_18[g, t]) * l_1_s1_18[g, t];
      target += log(chi_2_s1_18[g, t]) * l_2_s1_18[g, t];
      target += log(chi_3_s1_18[g, t]) * l_3_s1_18[g, t];
      target += log(chi_4_s1_18[g, t]) * l_4_s1_18[g, t];
      target += log(chi_5_s1_18[g, t]) * l_5_s1_18[g, t];
      target += log(chi_6_s1_18[g, t]) * l_6_s1_18[g, t];
      target += log(chi_7_s1_18[g, t]) * l_7_s1_18[g, t];
      target += log(chi_7_s3_18[g, t]) * l_7_s3_18[g, t];
      target += log(chi_7_s4_18[g, t]) * l_7_s4_18[g, t];
      target += log(chi_8_s1_18[g, t]) * l_8_s1_18[g, t];
      target += log(chi_8_s5_18[g, t]) * l_8_s5_18[g, t];
      target += log(chi_9_s1_18[g, t]) * l_9_s1_18[g, t];
      target += log(chi_9_s5_18[g, t]) * l_9_s5_18[g, t];
      target += log(chi_10_s1_18[g, t]) * l_10_s1_18[g, t];
    }
  }

  for (t in 1:max_T){
    for (u in t:max_T){
      {//2014
        target += (log(lambda_t1s1_t4s1_14[t, u]) + log(p_s1_t4_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t1s1_t4s1_14[t, u];
        target += (log(lambda_t1s1_t5s1_14[t, u]) + log(p_s1_t5_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t1s1_t5s1_14[t, u];

        target += (log(lambda_t4s1_t5s1_14[t, u]) + log(p_s1_t5_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t4s1_t5s1_14[t, u];

        target += (log(lambda_t5s1_t6s1_14[t, u]) + log(p_s1_t6_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t5s1_t6s1_14[t, u];
        target += (log(lambda_t5s1_t7s1_14[t, u]) + log(p_s1_t7_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t5s1_t7s1_14[t, u];
        target += (log(lambda_t5s1_t7s3_14[t, u]) + log(p_s3_t7_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t5s1_t7s3_14[t, u];
        target += (log(lambda_t5s1_t7s4_14[t, u]) + log(p_s4_t7_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t5s1_t7s4_14[t, u];
        target += (log(lambda_t5s1_t9s1_14[t, u]) + log(p_s1_t9_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t5s1_t9s1_14[t, u];

        target += (log(lambda_t6s1_t7s1_14[t, u]) + log(p_s1_t7_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t6s1_t7s1_14[t, u];
        target += (log(lambda_t6s1_t7s3_14[t, u]) + log(p_s3_t7_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t6s1_t7s3_14[t, u];
        target += (log(lambda_t6s1_t7s4_14[t, u]) + log(p_s4_t7_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t6s1_t7s4_14[t, u];
        target += (log(lambda_t6s1_t9s1_14[t, u]) + log(p_s1_t9_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t6s1_t9s1_14[t, u];

        target += (log(lambda_t7s1_t8s1_14[t, u]) + log(p_s1_t8_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t7s1_t8s1_14[t, u];
        target += (log(lambda_t7s1_t8s5_14[t, u]) + log(p_s5_t8_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t7s1_t8s5_14[t, u];
        target += (log(lambda_t7s1_t9s1_14[t, u]) + log(p_s1_t9_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t7s1_t9s1_14[t, u];

        target += (log(lambda_t7s3_t9s1_14[t, u]) + log(p_s1_t9_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t7s3_t9s1_14[t, u];
        target += (log(lambda_t7s4_t9s1_14[t, u]) + log(p_s1_t9_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t7s4_t9s1_14[t, u];

        target += (log(lambda_t8s1_t9s1_14[t, u]) + log(p_s1_t9_14[u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t8s1_t9s1_14[t, u];
        target += (log(lambda_t8s1_t11s1_14[t, u])  + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t8s1_t11s1_14[t, u];

        target += (log(lambda_t9s1_t11s1_14[t, u]) + log(tag_surv_prob_14[u]) - log(tag_surv_prob_14[t])) * m_t9s1_t11s1_14[t, u];
      }

      {//2015
        target += (log(lambda_t1s1_t4s1_15[t, u]) + log(p_s1_t4_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t1s1_t4s1_15[t, u];
        target += (log(lambda_t1s1_t5s1_15[t, u]) + log(p_s1_t5_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t1s1_t5s1_15[t, u];
        target += (log(lambda_t1s1_t7s1_15[t, u]) + log(p_s1_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t1s1_t7s1_15[t, u];
        target += (log(lambda_t1s1_t7s3_15[t, u]) + log(p_s3_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t1s1_t7s3_15[t, u];
        target += (log(lambda_t1s1_t7s4_15[t, u]) + log(p_s4_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t1s1_t7s4_15[t, u];
        target += (log(lambda_t1s1_t8s1_15[t, u]) + log(p_s1_t8_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t1s1_t8s1_15[t, u];

        target += (log(lambda_t4s1_t5s1_15[t, u]) + log(p_s1_t5_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t4s1_t5s1_15[t, u];
        if (u >= 21){
          target += (log(lambda_t4s1_t6s1_15[t, u]) + log(p_s1_t6_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t4s1_t6s1_15[t, u];
        }
        target += (log(lambda_t4s1_t7s1_15[t, u]) + log(p_s1_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t4s1_t7s1_15[t, u];
        target += (log(lambda_t4s1_t7s3_15[t, u]) + log(p_s3_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t4s1_t7s3_15[t, u];
        target += (log(lambda_t4s1_t7s4_15[t, u]) + log(p_s4_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t4s1_t7s4_15[t, u];
        target += (log(lambda_t4s1_t8s1_15[t, u]) + log(p_s1_t8_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t4s1_t8s1_15[t, u];
        target += (log(lambda_t4s1_t8s5_15[t, u]) + log(p_s5_t8_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t4s1_t8s5_15[t, u];
        target += (log(lambda_t4s1_t9s1_15[t, u]) + log(p_s1_t9_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t4s1_t9s1_15[t, u];
        if (u >= 21){
          target += (log(lambda_t5s1_t6s1_15[t, u]) + log(p_s1_t6_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t5s1_t6s1_15[t, u];
        }
        target += (log(lambda_t5s1_t7s1_15[t, u]) + log(p_s1_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t5s1_t7s1_15[t, u];
        target += (log(lambda_t5s1_t7s3_15[t, u]) + log(p_s3_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t5s1_t7s3_15[t, u];
        target += (log(lambda_t5s1_t7s4_15[t, u]) + log(p_s4_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t5s1_t7s4_15[t, u];
        target += (log(lambda_t5s1_t8s1_15[t, u]) + log(p_s1_t8_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t5s1_t8s1_15[t, u];
        target += (log(lambda_t5s1_t8s5_15[t, u]) + log(p_s5_t8_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t5s1_t8s5_15[t, u];
        target += (log(lambda_t5s1_t9s1_15[t, u]) + log(p_s1_t9_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t5s1_t9s1_15[t, u];

        target += (log(lambda_t6s1_t7s1_15[t, u]) + log(p_s1_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t6s1_t7s1_15[t, u];
        target += (log(lambda_t6s1_t7s3_15[t, u]) + log(p_s3_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t6s1_t7s3_15[t, u];
        target += (log(lambda_t6s1_t7s4_15[t, u]) + log(p_s4_t7_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t6s1_t7s4_15[t, u];

        target += (log(lambda_t7s1_t8s1_15[t, u]) + log(p_s1_t8_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t7s1_t8s1_15[t, u];
        target += (log(lambda_t7s1_t8s5_15[t, u]) + log(p_s5_t8_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t7s1_t8s5_15[t, u];
        target += (log(lambda_t7s1_t9s1_15[t, u]) + log(p_s1_t9_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t7s1_t9s1_15[t, u];

        target += (log(lambda_t7s3_t9s1_15[t, u]) + log(p_s1_t9_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t7s3_t9s1_15[t, u];
        target += (log(lambda_t7s3_t11s1_15[t, u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t7s3_t11s1_15[t, u];

        target += (log(lambda_t7s4_t9s1_15[t, u]) + log(p_s1_t9_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t7s4_t9s1_15[t, u];
        target += (log(lambda_t7s4_t11s1_15[t, u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t7s4_t11s1_15[t, u];

        target += (log(lambda_t8s1_t9s1_15[t, u]) + log(p_s1_t9_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t8s1_t9s1_15[t, u];
        target += (log(lambda_t8s1_t11s1_15[t, u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t8s1_t11s1_15[t, u];

        target += (log(lambda_t8s5_t9s5_15[t, u]) + log(p_s5_t9_15[u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t8s5_t9s5_15[t, u];
        target += (log(lambda_t8s5_t11s1_15[t, u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t8s5_t11s1_15[t, u];

        target += (log(lambda_t9s1_t11s1_15[t, u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t9s1_t11s1_15[t, u];

        target += (log(lambda_t9s5_t11s1_15[t, u]) + log(tag_surv_prob_15[u]) - log(tag_surv_prob_15[t])) * m_t9s5_t11s1_15[t, u];
      }

      {//2016
        if(m_t1s1_t4s2_16[t, u] > 0)
          target += (log(lambda_t1s1_t4s2_16[t, u]) + log(p_s2_t4_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t1s1_t4s2_16[t, u];
        target += (log(lambda_t1s1_t4s1_16[t, u]) + log(p_s1_t4_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t1s1_t4s1_16[t, u];
        target += (log(lambda_t1s1_t5s1_16[t, u]) + log(p_s1_t5_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t1s1_t5s1_16[t, u];
        target += (log(lambda_t1s1_t7s4_16[t, u]) + log(p_s4_t7_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t1s1_t7s4_16[t, u];

        target += (log(lambda_t4s2_t9s1_16[t, u]) + log(p_s1_t9_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t4s2_t9s1_16[t, u];

        target += (log(lambda_t4s1_t5s1_16[t, u]) + log(p_s1_t5_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t4s1_t5s1_16[t, u];
        target += (log(lambda_t4s1_t6s1_16[t, u]) + log(p_s1_t6_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t4s1_t6s1_16[t, u];
        target += (log(lambda_t4s1_t7s1_16[t, u]) + log(p_s1_t7_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t4s1_t7s1_16[t, u];

        target += (log(lambda_t5s1_t6s1_16[t, u]) + log(p_s1_t6_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t5s1_t6s1_16[t, u];
        target += (log(lambda_t5s1_t7s4_16[t, u]) + log(p_s4_t7_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t5s1_t7s4_16[t, u];
        target += (log(lambda_t5s1_t9s1_16[t, u]) + log(p_s1_t9_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t5s1_t9s1_16[t, u];

        target += (log(lambda_t6s1_t7s1_16[t, u]) + log(p_s1_t7_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t6s1_t7s1_16[t, u];
        target += (log(lambda_t6s1_t7s3_16[t, u]) + log(p_s3_t7_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t6s1_t7s3_16[t, u];
        target += (log(lambda_t6s1_t7s4_16[t, u]) + log(p_s4_t7_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t6s1_t7s4_16[t, u];
        target += (log(lambda_t6s1_t9s1_16[t, u]) + log(p_s1_t9_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t6s1_t9s1_16[t, u];

        target += (log(lambda_t7s1_t9s1_16[t, u]) + log(p_s1_t9_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t7s1_t9s1_16[t, u];
        target += (log(lambda_t7s1_t9s5_16[t, u]) + log(p_s5_t9_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t7s1_t9s5_16[t, u];
        target += (log(lambda_t7s1_t10s1_16[t, u]) + log(p_s1_t10_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t7s1_t10s1_16[t, u];

        target += (log(lambda_t7s3_t9s1_16[t, u]) + log(p_s1_t9_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t7s3_t9s1_16[t, u];
        target += (log(lambda_t7s3_t10s1_16[t, u]) + log(p_s1_t10_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t7s3_t10s1_16[t, u];

        target += (log(lambda_t7s4_t9s1_16[t, u]) + log(p_s1_t9_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t7s4_t9s1_16[t, u];
        target += (log(lambda_t7s4_t10s1_16[t, u]) + log(p_s1_t10_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t7s4_t10s1_16[t, u];

        target += (log(lambda_t9s1_t10s1_16[t, u]) + log(p_s1_t10_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t9s1_t10s1_16[t, u];
        target += (log(lambda_t9s1_t11s1_16[t, u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t9s1_t11s1_16[t, u];

        target += (log(lambda_t9s5_t10s1_16[t, u]) + log(p_s1_t10_16[u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t9s5_t10s1_16[t, u];
        target += (log(lambda_t9s5_t11s1_16[t, u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t9s5_t11s1_16[t, u];

        target += (log(chipps2benecia) + log(alpha_s1_t10_16[t, u]) + log(tag_surv_prob_16[u]) - log(tag_surv_prob_16[t])) * m_t10s1_t11s1_16[t, u];
      }

          
      {//2017
        target += (log(lambda_t1s1_t2s1_17[t, u]) + log(p_s1_t2_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t1s1_t2s1_17[t, u];
        target += (log(lambda_t1s1_t3s1_17[t, u]) + log(p_s1_t3_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t1s1_t3s1_17[t, u];
        if(m_t1s1_t4s2_17[t, u] > 0)
          target += (log(lambda_t1s1_t4s2_17[t, u]) + log(p_s2_t4_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t1s1_t4s2_17[t, u];
        target += (log(lambda_t1s1_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t1s1_t9s1_17[t, u];

        target += (log(lambda_t2s1_t3s1_17[t, u]) + log(p_s1_t3_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t2s1_t3s1_17[t, u];

        if(m_t2s1_t4s2_17[t, u] > 0)
          target += (log(lambda_t2s1_t4s2_17[t, u]) + log(p_s2_t4_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t2s1_t4s2_17[t, u];
          
        if (u >= 54){
          target += (log(lambda_t2s1_t4s1_17[t, u]) + log(p_s1_t4_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t2s1_t4s1_17[t, u];
        }
        target += (log(lambda_t2s1_t5s1_17[t, u]) + log(p_s1_t5_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t2s1_t5s1_17[t, u];
        target += (log(lambda_t2s1_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t2s1_t9s1_17[t, u];
        target += (log(lambda_t2s1_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t2s1_t10s1_17[t, u];
        
        if (u >= 54){
          target += (log(lambda_t3s1_t4s1_17[t, u]) + log(p_s1_t4_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t4s1_17[t, u];
        }
        target += (log(lambda_t3s1_t5s1_17[t, u]) + log(p_s1_t5_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t5s1_17[t, u];
        target += (log(lambda_t3s1_t6s1_17[t, u]) + log(p_s1_t6_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t6s1_17[t, u];
        target += (log(lambda_t3s1_t7s1_17[t, u]) + log(p_s1_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t7s1_17[t, u];
        if (u >= 35){
          target += (log(lambda_t3s1_t7s3_17[t, u]) + log(p_s3_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t7s3_17[t, u];
        }
        target += (log(lambda_t3s1_t7s4_17[t, u]) + log(p_s4_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t7s4_17[t, u];
        if (u <= 40){
          target += (log(lambda_t3s1_t8s1_17[t, u]) + log(p_s1_t8_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t8s1_17[t, u];
        }
        target += (log(lambda_t3s1_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t9s1_17[t, u];
        target += (log(lambda_t3s1_t9s5_17[t, u]) + log(p_s5_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t9s5_17[t, u];
        target += (log(lambda_t3s1_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t3s1_t10s1_17[t, u];

        target += (log(lambda_t4s2_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t4s2_t9s1_17[t, u];
        target += (log(lambda_t4s2_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t4s2_t10s1_17[t, u];
        target += (log(lambda_t4s2_t11s1_17[t, u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t4s2_t11s1_17[t, u];

        target += (log(lambda_t4s1_t5s1_17[t, u]) + log(p_s1_t5_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t4s1_t5s1_17[t, u];
        target += (log(lambda_t4s1_t6s1_17[t, u]) + log(p_s1_t6_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t4s1_t6s1_17[t, u];
        target += (log(lambda_t4s1_t7s1_17[t, u]) + log(p_s1_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t4s1_t7s1_17[t, u];

        target += (log(lambda_t5s1_t6s1_17[t, u]) + log(p_s1_t6_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t5s1_t6s1_17[t, u];
        target += (log(lambda_t5s1_t7s1_17[t, u]) + log(p_s1_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t5s1_t7s1_17[t, u];
        if (u >= 35){
          target += (log(lambda_t5s1_t7s3_17[t, u]) + log(p_s3_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t5s1_t7s3_17[t, u];
        }
        target += (log(lambda_t5s1_t7s4_17[t, u]) + log(p_s4_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t5s1_t7s4_17[t, u];
        if (u <= 40){
          target += (log(lambda_t5s1_t8s1_17[t, u]) + log(p_s1_t8_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t5s1_t8s1_17[t, u];
        }
        target += (log(lambda_t5s1_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t5s1_t9s1_17[t, u];
        target += (log(lambda_t5s1_t9s5_17[t, u]) + log(p_s5_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t5s1_t9s5_17[t, u];
        target += (log(lambda_t5s1_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t5s1_t10s1_17[t, u];

        target += (log(lambda_t6s1_t7s1_17[t, u]) + log(p_s1_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t6s1_t7s1_17[t, u];
        if (u >= 35){
          target += (log(lambda_t6s1_t7s3_17[t, u]) + log(p_s3_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t6s1_t7s3_17[t, u];
        }
        target += (log(lambda_t6s1_t7s4_17[t, u]) + log(p_s4_t7_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t6s1_t7s4_17[t, u];
        if (u <= 40){
          target += (log(lambda_t6s1_t8s1_17[t, u]) + log(p_s1_t8_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t6s1_t8s1_17[t, u];
        }
        target += (log(lambda_t6s1_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t6s1_t9s1_17[t, u];
        target += (log(lambda_t6s1_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t6s1_t10s1_17[t, u];

        if (u <= 40){
          target += (log(lambda_t7s1_t8s1_17[t, u]) + log(p_s1_t8_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s1_t8s1_17[t, u];
        }
        target += (log(lambda_t7s1_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s1_t9s1_17[t, u];
        target += (log(lambda_t7s1_t9s5_17[t, u]) + log(p_s5_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s1_t9s5_17[t, u];
        target += (log(lambda_t7s1_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s1_t10s1_17[t, u];

        target += (log(lambda_t7s3_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s3_t9s1_17[t, u];
        target += (log(lambda_t7s3_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s3_t10s1_17[t, u];
        target += (log(lambda_t7s3_t11s1_17[t, u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s3_t11s1_17[t, u];

        target += (log(lambda_t7s4_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s4_t9s1_17[t, u];
        target += (log(lambda_t7s4_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s4_t10s1_17[t, u];
        target += (log(lambda_t7s4_t11s1_17[t, u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t7s4_t11s1_17[t, u];

        target += (log(lambda_t8s1_t9s1_17[t, u]) + log(p_s1_t9_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t8s1_t9s1_17[t, u];
        target += (log(lambda_t8s1_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t8s1_t10s1_17[t, u];

        target += (log(lambda_t9s1_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t9s1_t10s1_17[t, u];
        target += (log(lambda_t9s1_t11s1_17[t, u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t9s1_t11s1_17[t, u];

        target += (log(lambda_t9s5_t10s1_17[t, u]) + log(p_s1_t10_17[u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t9s5_t10s1_17[t, u];

        target += (log(chipps2benecia) + log(alpha_s1_t10_17[t, u]) + log(tag_surv_prob_17[u]) - log(tag_surv_prob_17[t])) * m_t10s1_t11s1_17[t, u];
      }
          
      for (g in 1:2){//2018
        if (m_t1s1_t2s1_18[g, t, u] > 0)
          target += (log(lambda_t1s1_t2s1_18[t, u]) + log(p_s1_t2_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t1s1_t2s1_18[g, t, u];
        if (m_t1s1_t3s1_18[g, t, u] > 0)
          target += (log(lambda_t1s1_t3s1_18[t, u]) + log(p_s1_t3_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t1s1_t3s1_18[g, t, u];
        target += (log(lambda_t1s1_t4s1_18[t, u]) + log(p_s1_t4_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t1s1_t4s1_18[g, t, u];
        target += (log(lambda_t1s1_t5s1_18[t, u]) + log(p_s1_t5_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t1s1_t5s1_18[g, t, u];
        target += (log(lambda_t1s1_t6s1_18[t, u]) + log(p_s1_t6_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t1s1_t6s1_18[g, t, u];
        target += (log(lambda_t1s1_t7s1_18[t, u]) + log(p_s1_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t1s1_t7s1_18[g, t, u];
        if (m_t2s1_t3s1_18[g, t, u] > 0)
          target += (log(lambda_t2s1_t3s1_18[t, u]) + log(p_s1_t3_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t2s1_t3s1_18[g, t, u];

        target += (log(lambda_t2s1_t4s1_18[t, u]) + log(p_s1_t4_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t2s1_t4s1_18[g, t, u];
        target += (log(lambda_t2s1_t5s1_18[t, u]) + log(p_s1_t5_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t2s1_t5s1_18[g, t, u];
        target += (log(lambda_t2s1_t6s1_18[t, u]) + log(p_s1_t6_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t2s1_t6s1_18[g, t, u];
        target += (log(lambda_t2s1_t7s1_18[t, u]) + log(p_s1_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t2s1_t7s1_18[g, t, u];
        target += (log(lambda_t2s1_t7s3_18[t, u]) + log(p_s3_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t2s1_t7s3_18[g, t, u];
        target += (log(lambda_t2s1_t7s4_18[t, u]) + log(p_s4_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t2s1_t7s4_18[g, t, u];
        target += (log(lambda_t2s1_t9s1_18[t, u]) + log(p_s1_t9_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t2s1_t9s1_18[g, t, u];

        target += (log(lambda_t3s1_t4s1_18[t, u]) + log(p_s1_t4_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t3s1_t4s1_18[g, t, u];
        target += (log(lambda_t3s1_t5s1_18[t, u]) + log(p_s1_t5_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t3s1_t5s1_18[g, t, u];
        target += (log(lambda_t3s1_t6s1_18[t, u]) + log(p_s1_t6_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t3s1_t6s1_18[g, t, u];
        target += (log(lambda_t3s1_t7s1_18[t, u]) + log(p_s1_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t3s1_t7s1_18[g, t, u];
        target += (log(lambda_t3s1_t9s1_18[t, u]) + log(p_s1_t9_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t3s1_t9s1_18[g, t, u];

        target += (log(lambda_t4s1_t5s1_18[t, u]) + log(p_s1_t5_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t4s1_t5s1_18[g, t, u];
        target += (log(lambda_t4s1_t6s1_18[t, u]) + log(p_s1_t6_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t4s1_t6s1_18[g, t, u];
        target += (log(lambda_t4s1_t7s1_18[t, u]) + log(p_s1_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t4s1_t7s1_18[g, t, u];
        target += (log(lambda_t4s1_t7s3_18[t, u]) + log(p_s3_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t4s1_t7s3_18[g, t, u];
        target += (log(lambda_t4s1_t7s4_18[t, u]) + log(p_s4_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t4s1_t7s4_18[g, t, u];

        target += (log(lambda_t5s1_t6s1_18[t, u]) + log(p_s1_t6_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t5s1_t6s1_18[g, t, u];
        target += (log(lambda_t5s1_t7s1_18[t, u]) + log(p_s1_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t5s1_t7s1_18[g, t, u];
        target += (log(lambda_t5s1_t7s3_18[t, u]) + log(p_s3_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t5s1_t7s3_18[g, t, u];
        target += (log(lambda_t5s1_t7s4_18[t, u]) + log(p_s4_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t5s1_t7s4_18[g, t, u];

        target += (log(lambda_t6s1_t7s1_18[t, u]) + log(p_s1_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t6s1_t7s1_18[g, t, u];
        target += (log(lambda_t6s1_t7s3_18[t, u]) + log(p_s3_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t6s1_t7s3_18[g, t, u];
        target += (log(lambda_t6s1_t7s4_18[t, u]) + log(p_s4_t7_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t6s1_t7s4_18[g, t, u];
        target += (log(lambda_t6s1_t8s5_18[t, u]) + log(p_s5_t8_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t6s1_t8s5_18[g, t, u];
        
        if (m_t7s1_t8s1_18[g, t, u] > 0)
          target += (log(lambda_t7s1_t8s1_18[t, u]) + log(p_s1_t8_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t7s1_t8s1_18[g, t, u];

        target += (log(lambda_t7s1_t8s5_18[t, u]) + log(p_s5_t8_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t7s1_t8s5_18[g, t, u];
        target += (log(lambda_t7s1_t9s1_18[t, u]) + log(p_s1_t9_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t7s1_t9s1_18[g, t, u];
        target += (log(lambda_t7s1_t10s1_18[t, u]) + log(p_s1_t10_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t7s1_t10s1_18[g, t, u];
        
        target += (log(lambda_t7s3_t9s1_18[t, u]) + log(p_s1_t9_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t7s3_t9s1_18[g, t, u];
        target += (log(lambda_t7s3_t10s1_18[t, u]) + log(p_s1_t10_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t7s3_t10s1_18[g, t, u];
        
        target += (log(lambda_t7s4_t9s1_18[t, u]) + log(p_s1_t9_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t7s4_t9s1_18[g, t, u];
        target += (log(lambda_t7s4_t10s1_18[t, u]) + log(p_s1_t10_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t7s4_t10s1_18[g, t, u];

        target += (log(lambda_t8s1_t9s1_18[t, u]) + log(p_s1_t9_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t8s1_t9s1_18[g, t, u];
        target += (log(lambda_t8s1_t10s1_18[t, u]) + log(p_s1_t10_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t8s1_t10s1_18[g, t, u];

        target += (log(lambda_t8s5_t9s5_18[t, u]) + log(p_s5_t9_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t8s5_t9s5_18[g, t, u];

        target += (log(lambda_t9s1_t10s1_18[t, u]) + log(p_s1_t10_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t9s1_t10s1_18[g, t, u];
        target += (log(lambda_t9s1_t11s1_18[t, u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t9s1_t11s1_18[g, t, u];

        target += (log(lambda_t9s5_t10s1_18[t, u]) + log(p_s1_t10_18[u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t9s5_t10s1_18[g, t, u];

        target += (log(chipps2benecia) + log(alpha_s1_t10_18[t, u]) + log(tag_surv_prob_18[g, u]) - log(tag_surv_prob_18[g, t])) * m_t10s1_t11s1_18[g, t, u];
      }

    }
    

  }
}

generated quantities{
  //These relate to the actual releases
  real overall_survival_14;
  real overall_survival_15;
  real overall_survival_16;
  real overall_survival_17;
  real overall_survival_18;

  real prop_sacra_14; //proportion of survivors that remained in sacramento
  real prop_sutter_14; //proportion of survivors that took sutter
  real prop_steam_14;  //proportion of survivors that took steam
  real prop_georg_14; //proportion of survivors that took georgiana

  real prop_sacra_15;
  real prop_sutter_15;
  real prop_steam_15;
  real prop_georg_15;

  real prop_sacra_16;
  real prop_yolo_16;
  real prop_sutter_16;
  real prop_steam_16;
  real prop_georg_16;

  real prop_sacra_17;
  real prop_yolo_17;
  real prop_sutter_17;
  real prop_steam_17;
  real prop_georg_17;

  real prop_sacra_18;
  real prop_yolo_18;
  real prop_sutter_18;
  real prop_steam_18;
  real prop_georg_18;

  row_vector[T] chipps_arrival_dist_14;
  row_vector[T] chipps_arrival_dist_15;
  row_vector[T] chipps_arrival_dist_16;
  row_vector[T] chipps_arrival_dist_17;
  row_vector[T] chipps_arrival_dist_18;

  //These relate to fish-at-large passing KL each day
  vector[T] overall_survival_from_KL_14;
  vector[T] overall_survival_from_KL_15;
  vector[T] overall_survival_from_KL_16;
  vector[T] overall_survival_from_KL_17;
  vector[T] overall_survival_from_KL_18;

  vector[T] route_survival_from_KL_14[5];
  vector[T] route_survival_from_KL_15[5];
  vector[T] route_survival_from_KL_16[5];
  vector[T] route_survival_from_KL_17[5];
  vector[T] route_survival_from_KL_18[5];

  vector[T] route_probability_from_KL_14[5];
  vector[T] route_probability_from_KL_15[5];
  vector[T] route_probability_from_KL_16[5];
  vector[T] route_probability_from_KL_17[5];
  vector[T] route_probability_from_KL_18[5];
  //2014 Gen Quants ----\\
  {
    row_vector[T] temp_surv;
    row_vector[T] temp1;
    row_vector[T] temp3;
    row_vector[T] temp4;
    row_vector[T] temp5;

    temp1 = phi_s1_t0_14 * alpha_s1_t0_14; //t1s1

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t1_14') * alpha_s1_t1_14; //t2s1

    temp_surv = temp1;
    temp1 = temp_surv * alpha_s1_t2_14; //t3s1

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t3_14') * alpha_s1_t3_14; //t4s1

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t4_14') * alpha_s1_t4_14; //t5s1

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t5_14')  * alpha_s1_t5_14; //t6s1

    temp_surv = (temp1 .* phi_s1_t6_14');
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_14, ((1 - psi_1to3_t6_14) .* (1 - psi_1to4_t6_14))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_14, psi_1to3_t6_14); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_14, (1 - psi_1to3_t6_14) .* psi_1to4_t6_14); //t7s4

    temp_surv = temp1 .* phi_s1_t7_14';
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_14, (1 - psi_1to5_t7_14)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_14, psi_1to5_t7_14);       //t8s5

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t8_14') * alpha_s1_t8_14; //t9s1 by sacramento
    temp_surv = temp3;
    temp3 = (temp_surv .* phi_s3_t8_14') * alpha_s3_t8_14; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = (temp_surv .* phi_s4_t8_14') * alpha_s4_t8_14; //t9s1 by steam
    temp_surv = temp5;
    temp5 = (temp_surv .* phi_s5_t8_14') * alpha_s5_t8_14; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t9_14') * alpha_s1_t9_14; //t10s1 by sacramento
    temp_surv = temp3;
    temp3 = (temp_surv .* phi_s1_t9_14') * alpha_s1_t9_14; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = (temp_surv .* phi_s1_t9_14') * alpha_s1_t9_14; //t10s1 by steam
    temp_surv = temp5;
    temp5 = (temp_surv .* phi_s5_t9_14') * alpha_s5_t9_14; //t10s1 by georgiana

    prop_sacra_14       = sum(temp1);
    prop_sutter_14      = sum(temp3);
    prop_steam_14       = sum(temp4);
    prop_georg_14       = sum(temp5);
    overall_survival_14 = prop_sacra_14 + prop_sutter_14 + prop_steam_14 + prop_georg_14;
    chipps_arrival_dist_14 = (temp1 + temp3 + temp4 + temp5)/ overall_survival_14;
  }

  {
    matrix[T, T] temp_surv;
    matrix[T, T] temp1;
    matrix[T, T] temp3;
    matrix[T, T] temp4;
    matrix[T, T] temp5;

    temp1 = diag_pre_multiply(phi_s1_t1_14, alpha_s1_t1_14); //t2s1

    temp_surv = temp1;
    temp1 = temp_surv * alpha_s1_t2_14; //t3s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t3_14) * alpha_s1_t3_14; //t4s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t4_14) * alpha_s1_t4_14; //t5s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t5_14)  * alpha_s1_t5_14; //t6s1

    temp_surv = diag_post_multiply(temp1, phi_s1_t6_14) * alpha_s1_t6_14;
    temp1 = temp_surv; //t7s1
    temp3 = temp_surv; //t7s3
    temp4 = temp_surv; //t7s4

    temp_surv = diag_post_multiply(temp1, phi_s1_t7_14) * alpha_s1_t7_14;
    temp1 = temp_surv; //t8s1
    temp5 = temp_surv; //t8s5

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t8_14) * alpha_s1_t8_14; //t9s1 by sacramento
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s3_t8_14) * alpha_s3_t8_14; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s4_t8_14) * alpha_s4_t8_14; //t9s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t8_14) * alpha_s5_t8_14; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t9_14) * alpha_s1_t9_14; //t10s1 by sacramento
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s1_t9_14) * alpha_s1_t9_14; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s1_t9_14) * alpha_s1_t9_14; //t10s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t9_14) * alpha_s5_t9_14; //t10s1 by georgiana

    route_survival_from_KL_14[1] = temp1 * rep_vector(1, T);
    route_survival_from_KL_14[2] = rep_vector(0, T);
    route_survival_from_KL_14[3] = temp3 * rep_vector(1, T);
    route_survival_from_KL_14[4] = temp4 * rep_vector(1, T);
    route_survival_from_KL_14[5] = temp5 * rep_vector(1, T);

    temp_surv = alpha_s1_t1_14 * alpha_s1_t2_14 * alpha_s1_t3_14 * alpha_s1_t5_14;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_14, ((1 - psi_1to3_t6_14) .* (1 - psi_1to4_t6_14))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_14, psi_1to3_t6_14); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_14, (1 - psi_1to3_t6_14) .* psi_1to4_t6_14); //t7s4

    temp_surv = temp1;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_14, (1 - psi_1to5_t7_14)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_14, psi_1to5_t7_14);       //t8s5

    route_probability_from_KL_14[1] = temp1 * rep_vector(1, T);
    route_probability_from_KL_14[2] = rep_vector(0, T);
    route_probability_from_KL_14[3] = temp3 * rep_vector(1, T);
    route_probability_from_KL_14[4] = temp4 * rep_vector(1, T);
    route_probability_from_KL_14[5] = temp5 * rep_vector(1, T);
    
    overall_survival_from_KL_14 = route_survival_from_KL_14[1] .* route_probability_from_KL_14[1] +
                                  //route_survival_from_KL_14[2] .* route_probability_from_KL_14[2] +
                                  route_survival_from_KL_14[3] .* route_probability_from_KL_14[3] +
                                  route_survival_from_KL_14[4] .* route_probability_from_KL_14[4] + 
                                  route_survival_from_KL_14[5] .* route_probability_from_KL_14[5];
  }
  //2015 Gen Quants ----\\
  {
    matrix[2, T] temp_surv;
    matrix[2, T] temp1;
    matrix[2, T] temp3;
    matrix[2, T] temp4;
    matrix[2, T] temp5;

    temp1 = diag_pre_multiply(phi_s1_t0_15, alpha_s1_t0_15); //t1s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t1_15) * alpha_s1_t1_15; //t2s1

    temp_surv = temp1;
    temp1 = temp_surv * alpha_s1_t2_15; //t3s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t3_15) * alpha_s1_t3_15; //t4s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t4_15) * alpha_s1_t4_15; //t5s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t5_15)  * alpha_s1_t5_15; //t6s1

    temp_surv = diag_post_multiply(temp1, phi_s1_t6_15);
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_15, ((1 - psi_1to3_t6_15) .* (1 - psi_1to4_t6_15))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_15, psi_1to3_t6_15); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_15, (1 - psi_1to3_t6_15) .* psi_1to4_t6_15); //t7s4

    temp_surv = diag_post_multiply(temp1, phi_s1_t7_15);
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_15, (1 - psi_1to5_t7_15)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_15, psi_1to5_t7_15);       //t8s5

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t8_15) * alpha_s1_t8_15; //t9s1 by sacramento
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s3_t8_15) * alpha_s3_t8_15; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s4_t8_15) * alpha_s4_t8_15; //t9s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t8_15) * alpha_s5_t8_15; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t9_15) * alpha_s1_t9_15; //t10s1 by sacramento
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s1_t9_15) * alpha_s1_t9_15; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s1_t9_15) * alpha_s1_t9_15; //t10s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t9_15) * alpha_s5_t9_15; //t10s1 by georgiana

    prop_sacra_15    = (sum(temp1[1, ]) * N_rls_grp_15[1] + sum(temp1[2, ]) * N_rls_grp_15[2])/N_15;
    prop_sutter_15   = (sum(temp3[1, ]) * N_rls_grp_15[1] + sum(temp3[2, ]) * N_rls_grp_15[2])/N_15;
    prop_steam_15    = (sum(temp4[1, ]) * N_rls_grp_15[1] + sum(temp4[2, ]) * N_rls_grp_15[2])/N_15;
    prop_georg_15    = (sum(temp5[1, ]) * N_rls_grp_15[1] + sum(temp5[2, ]) * N_rls_grp_15[2])/N_15;
    overall_survival_15 = prop_sacra_15 + prop_sutter_15 + prop_steam_15 + prop_georg_15;
    chipps_arrival_dist_15 = (((temp1[1, ] + temp3[1, ] + temp4[1, ] + temp5[1, ]) * N_rls_grp_15[1] +
                              (temp1[2, ] + temp3[2, ] + temp4[2, ] + temp5[2, ]) * N_rls_grp_15[2])/ N_15) / overall_survival_15;
  }

  {
    matrix[T, T] temp_surv;
    matrix[T, T] temp1;
    matrix[T, T] temp3;
    matrix[T, T] temp4;
    matrix[T, T] temp5;

    temp1 = diag_pre_multiply(phi_s1_t1_15, alpha_s1_t1_15); //t2s1

    temp_surv = temp1;
    temp1 = temp_surv * alpha_s1_t2_15; //t3s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t3_15) * alpha_s1_t3_15; //t4s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t4_15) * alpha_s1_t4_15; //t5s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t5_15)  * alpha_s1_t5_15; //t6s1

    temp_surv = diag_post_multiply(temp1, phi_s1_t6_15) * alpha_s1_t6_15;
    temp1 = temp_surv; //t7s1
    temp3 = temp_surv; //t7s3
    temp4 = temp_surv; //t7s4

    temp_surv = diag_post_multiply(temp1, phi_s1_t7_15) * alpha_s1_t7_15;
    temp1 = temp_surv; //t8s1
    temp5 = temp_surv;       //t8s5

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t8_15) * alpha_s1_t8_15; //t9s1 by sacramento
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s3_t8_15) * alpha_s3_t8_15; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s4_t8_15) * alpha_s4_t8_15; //t9s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t8_15) * alpha_s5_t8_15; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t9_15) * alpha_s1_t9_15; //t10s1 by sacramento
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s1_t9_15) * alpha_s1_t9_15; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s1_t9_15) * alpha_s1_t9_15; //t10s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t9_15) * alpha_s5_t9_15; //t10s1 by georgiana

    route_survival_from_KL_15[1] = temp1 * rep_vector(1, T);
    route_survival_from_KL_15[2] = rep_vector(0, T);
    route_survival_from_KL_15[3] = temp3 * rep_vector(1, T);
    route_survival_from_KL_15[4] = temp4 * rep_vector(1, T);
    route_survival_from_KL_15[5] = temp5 * rep_vector(1, T);

    temp_surv = alpha_s1_t1_15 * alpha_s1_t2_15 * alpha_s1_t3_15 * alpha_s1_t5_15;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_15, ((1 - psi_1to3_t6_15) .* (1 - psi_1to4_t6_15))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_15, psi_1to3_t6_15); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_15, (1 - psi_1to3_t6_15) .* psi_1to4_t6_15); //t7s4

    temp_surv = temp1;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_15, (1 - psi_1to5_t7_15)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_15, psi_1to5_t7_15);       //t8s5

    route_probability_from_KL_15[1] = temp1 * rep_vector(1, T);
    route_probability_from_KL_15[2] = rep_vector(0, T);
    route_probability_from_KL_15[3] = temp3 * rep_vector(1, T);
    route_probability_from_KL_15[4] = temp4 * rep_vector(1, T);
    route_probability_from_KL_15[5] = temp5 * rep_vector(1, T);
    
    overall_survival_from_KL_15 = route_survival_from_KL_15[1] .* route_probability_from_KL_15[1] +
                                  //route_survival_from_KL_15[2] .* route_probability_from_KL_15[2] +
                                  route_survival_from_KL_15[3] .* route_probability_from_KL_15[3] +
                                  route_survival_from_KL_15[4] .* route_probability_from_KL_15[4] + 
                                  route_survival_from_KL_15[5] .* route_probability_from_KL_15[5];
  }
  //2016 Gen Quants ----\\
  {
    matrix[2, T] temp_surv;
    matrix[2, T] temp1;
    matrix[2, T] temp2;
    matrix[2, T] temp3;
    matrix[2, T] temp4;
    matrix[2, T] temp5;

    temp1 = diag_pre_multiply(phi_s1_t0_16, alpha_s1_t0_16); //t1s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t1_16) * alpha_s1_t1_16; //t2s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, (1 - psi_1to2_t2_16)) * alpha_s1_t2_16; //t3s1
    temp2 = diag_post_multiply(temp_surv, (psi_1to2_t2_16 .* phi_s2_t3_16)) *  alpha_s2_t3_16;

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t3_16) * alpha_s1_t3_16; //t4s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t4_16) * alpha_s1_t4_16; //t5s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t5_16)  * alpha_s1_t5_16; //t6s1

    temp_surv = diag_post_multiply(temp1, phi_s1_t6_16);
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_16, ((1 - psi_1to3_t6_16) .* (1 - psi_1to4_t6_16))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_16, psi_1to3_t6_16); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_16, (1 - psi_1to3_t6_16) .* psi_1to4_t6_16); //t7s4

    temp_surv = diag_post_multiply(temp1, phi_s1_t7_16);
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_16, (1 - psi_1to5_t7_16)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_16, psi_1to5_t7_16);       //t8s5

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t8_16) * alpha_s1_t8_16; //t9s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s2_t8_16) * alpha_s2_t8_16; //t9s1 by yolo
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s3_t8_16) * alpha_s3_t8_16; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s4_t8_16) * alpha_s4_t8_16; //t9s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t8_16) * alpha_s5_t8_16; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t9_16) * alpha_s1_t9_16; //t10s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s1_t9_16) * alpha_s1_t9_16; //t10s1 by yolo
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s1_t9_16) * alpha_s1_t9_16; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s1_t9_16) * alpha_s1_t9_16; //t10s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t9_16) * alpha_s5_t9_16; //t10s1 by georgiana

    prop_sacra_16    = (sum(temp1[1, ]) * N_rls_grp_16[1] + sum(temp1[2, ]) * N_rls_grp_16[2])/N_16;
    prop_yolo_16     = (sum(temp2[1, ]) * N_rls_grp_16[1] + sum(temp2[2, ]) * N_rls_grp_16[2])/N_16;
    prop_sutter_16   = (sum(temp3[1, ]) * N_rls_grp_16[1] + sum(temp3[2, ]) * N_rls_grp_16[2])/N_16;
    prop_steam_16    = (sum(temp4[1, ]) * N_rls_grp_16[1] + sum(temp4[2, ]) * N_rls_grp_16[2])/N_16;
    prop_georg_16    = (sum(temp5[1, ]) * N_rls_grp_16[1] + sum(temp5[2, ]) * N_rls_grp_16[2])/N_16;
    overall_survival_16 = prop_sacra_16 + prop_yolo_16 + prop_sutter_16 + prop_steam_16 + prop_georg_16;
    chipps_arrival_dist_16 = (((temp1[1, ] + temp2[1, ] + temp3[1, ] + temp4[1, ] + temp5[1, ]) * N_rls_grp_16[1] +
                              (temp1[2, ]  + temp2[1, ] + temp3[2, ] + temp4[2, ] + temp5[2, ]) * N_rls_grp_16[2])/ N_16) / overall_survival_16;
  }

  {
    matrix[T, T] temp_surv;
    matrix[T, T] temp1;
    matrix[T, T] temp2;
    matrix[T, T] temp3;
    matrix[T, T] temp4;
    matrix[T, T] temp5;

    temp1 = diag_pre_multiply(phi_s1_t1_16, alpha_s1_t1_16); //t2s1

    temp_surv = temp1;
    temp1 = temp_surv * alpha_s1_t2_16; //t3s1
    temp2 = diag_post_multiply(temp_surv, phi_s2_t3_16) *  alpha_s2_t3_16; //t4s2

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t3_16) * alpha_s1_t3_16; //t4s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t4_16) * alpha_s1_t4_16; //t5s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t5_16)  * alpha_s1_t5_16; //t6s1

    temp_surv = diag_post_multiply(temp1, phi_s1_t6_16) * alpha_s1_t6_16;
    temp1 = temp_surv; //t7s1
    temp3 = temp_surv; //t7s3
    temp4 = temp_surv; //t7s4

    temp_surv = diag_post_multiply(temp1, phi_s1_t7_16) * alpha_s1_t7_16;
    temp1 = temp_surv; //t8s1
    temp5 = temp_surv;       //t8s5

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t8_16) * alpha_s1_t8_16; //t9s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s2_t8_16) * alpha_s2_t8_16; //t9s1 by yolo
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s3_t8_16) * alpha_s3_t8_16; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s4_t8_16) * alpha_s4_t8_16; //t9s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t8_16) * alpha_s5_t8_16; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t9_16) * alpha_s1_t9_16; //t10s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s1_t9_16) * alpha_s1_t9_16; //t10s1 by yolo
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s1_t9_16) * alpha_s1_t9_16; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s1_t9_16) * alpha_s1_t9_16; //t10s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t9_16) * alpha_s5_t9_16; //t10s1 by georgiana

    route_survival_from_KL_16[1] = temp1 * rep_vector(1, T);
    route_survival_from_KL_16[2] = temp2 * rep_vector(1, T);
    route_survival_from_KL_16[3] = temp3 * rep_vector(1, T);
    route_survival_from_KL_16[4] = temp4 * rep_vector(1, T);
    route_survival_from_KL_16[5] = temp5 * rep_vector(1, T);

    temp_surv = alpha_s1_t1_16;
    temp1 = diag_post_multiply(temp_surv, (1 - psi_1to2_t2_16));
    temp2 = diag_post_multiply(temp_surv, psi_1to2_t2_16);

    temp_surv = temp1 * alpha_s1_t2_16 * alpha_s1_t3_16 * alpha_s1_t5_16;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_16, ((1 - psi_1to3_t6_16) .* (1 - psi_1to4_t6_16))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_16, psi_1to3_t6_16); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_16, (1 - psi_1to3_t6_16) .*psi_1to4_t6_16); //t7s4

    temp_surv = temp1;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_16, (1 - psi_1to5_t7_16)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_16, psi_1to5_t7_16);       //t8s5

    route_probability_from_KL_16[1] = temp1 * rep_vector(1, T);
    route_probability_from_KL_16[2] = temp2 * rep_vector(1, T);
    route_probability_from_KL_16[3] = temp3 * rep_vector(1, T);
    route_probability_from_KL_16[4] = temp4 * rep_vector(1, T);
    route_probability_from_KL_16[5] = temp5 * rep_vector(1, T);
    
    overall_survival_from_KL_16 = route_survival_from_KL_16[1] .* route_probability_from_KL_16[1] +
                                  route_survival_from_KL_16[2] .* route_probability_from_KL_16[2] +
                                  route_survival_from_KL_16[3] .* route_probability_from_KL_16[3] +
                                  route_survival_from_KL_16[4] .* route_probability_from_KL_16[4] + 
                                  route_survival_from_KL_16[5] .* route_probability_from_KL_16[5];
  }
  //2017 Gen Quants ----\\
  {
    row_vector[T] temp_surv;
    row_vector[T] temp1;
    row_vector[T] temp2;
    row_vector[T] temp3;
    row_vector[T] temp4;
    row_vector[T] temp5;

    temp1 = phi_s1_t0_17 * alpha_s1_t0_17; //t1s1

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t1_17') * alpha_s1_t1_17; //t2s1

    temp_surv = temp1;
    temp1 = (temp_surv .* (1 - psi_1to2_t2_17)') * alpha_s1_t2_17; //t3s1
    temp2 = (temp_surv .* (psi_1to2_t2_17 .* phi_s2_t3_17)') *  alpha_s2_t3_17; //t4s2

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t3_17') * alpha_s1_t3_17; //t4s1

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t4_17') * alpha_s1_t4_17; //t5s1

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t5_17')  * alpha_s1_t5_17; //t6s1

    temp_surv = (temp1 .* phi_s1_t6_17');
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_17, psi_1to3_t6_17); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_17, (1 - psi_1to3_t6_17) .* psi_1to4_t6_17); //t7s4

    temp_surv = temp1 .* phi_s1_t7_17';
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_17, psi_1to5_t7_17);       //t8s5

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t8_17') * alpha_s1_t8_17; //t9s1 by sacramento
    temp_surv = temp2;
    temp2 = (temp_surv .* phi_s2_t8_17') * alpha_s2_t8_17; //t9s1 by yolo
    temp_surv = temp3;
    temp3 = (temp_surv .* phi_s3_t8_17') * alpha_s3_t8_17; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = (temp_surv .* phi_s4_t8_17') * alpha_s4_t8_17; //t9s1 by steam
    temp_surv = temp5;
    temp5 = (temp_surv .* phi_s5_t8_17') * alpha_s5_t8_17; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = (temp_surv .* phi_s1_t9_17') * alpha_s1_t9_17; //t10s1 by sacramento
    temp_surv = temp2;
    temp2 = (temp_surv .* phi_s1_t9_17') * alpha_s1_t9_17; //t10s1 by yolo
    temp_surv = temp3;
    temp3 = (temp_surv .* phi_s1_t9_17') * alpha_s1_t9_17; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = (temp_surv .* phi_s1_t9_17') * alpha_s1_t9_17; //t10s1 by steam
    temp_surv = temp5;
    temp5 = (temp_surv .* phi_s5_t9_17') * alpha_s5_t9_17; //t10s1 by georgiana

    prop_sacra_17       = sum(temp1);
    prop_yolo_17        = sum(temp2);
    prop_sutter_17      = sum(temp3);
    prop_steam_17       = sum(temp4);
    prop_georg_17       = sum(temp5);
    overall_survival_17 = prop_sacra_17 + prop_yolo_17 + prop_sutter_17 + prop_steam_17 + prop_georg_17;
    chipps_arrival_dist_17 = (temp1 + temp2 + temp3 + temp4 + temp5)/ overall_survival_17;
  }

  {
    matrix[T, T] temp_surv;
    matrix[T, T] temp1;
    matrix[T, T] temp2;
    matrix[T, T] temp3;
    matrix[T, T] temp4;
    matrix[T, T] temp5;

    temp1 = diag_pre_multiply(phi_s1_t1_17, alpha_s1_t1_17); //t2s1

    temp_surv = temp1;
    temp1 = temp_surv * alpha_s1_t2_17; //t3s1
    temp2 = diag_post_multiply(temp_surv, phi_s2_t3_17) *  alpha_s2_t3_17; //t4s2

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t3_17) * alpha_s1_t3_17; //t4s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t4_17) * alpha_s1_t4_17; //t5s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t5_17)  * alpha_s1_t5_17; //t6s1

    temp_surv = diag_post_multiply(temp1, phi_s1_t6_17) * alpha_s1_t6_17;
    temp1 = temp_surv; //t7s1
    temp3 = temp_surv; //t7s3
    temp4 = temp_surv; //t7s4

    temp_surv = diag_post_multiply(temp1, phi_s1_t7_17) * alpha_s1_t7_17;
    temp1 = temp_surv; //t8s1
    temp5 = temp_surv; //t8s5

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t8_17) * alpha_s1_t8_17; //t9s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s2_t8_17) * alpha_s2_t8_17; //t9s1 by yolo
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s3_t8_17) * alpha_s3_t8_17; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s4_t8_17) * alpha_s4_t8_17; //t9s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t8_17) * alpha_s5_t8_17; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t9_17) * alpha_s1_t9_17; //t10s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s1_t9_17) * alpha_s1_t9_17; //t10s1 by yolo
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s1_t9_17) * alpha_s1_t9_17; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s1_t9_17) * alpha_s1_t9_17; //t10s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t9_17) * alpha_s5_t9_17; //t10s1 by georgiana

    route_survival_from_KL_17[1] = temp1 * rep_vector(1, T);
    route_survival_from_KL_17[2] = temp2 * rep_vector(1, T);
    route_survival_from_KL_17[3] = temp3 * rep_vector(1, T);
    route_survival_from_KL_17[4] = temp4 * rep_vector(1, T);
    route_survival_from_KL_17[5] = temp5 * rep_vector(1, T);


    temp_surv = alpha_s1_t1_17;
    temp1 = diag_post_multiply(temp_surv, (1 - psi_1to2_t2_17));
    temp2 = diag_post_multiply(temp_surv, psi_1to2_t2_17);

    temp_surv = temp1 * alpha_s1_t2_17 * alpha_s1_t3_17 * alpha_s1_t5_17;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_17, ((1 - psi_1to3_t6_17) .* (1 - psi_1to4_t6_17))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_17, psi_1to3_t6_17); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_17, (1 - psi_1to3_t6_17) .* psi_1to4_t6_17); //t7s4

    temp_surv = temp1;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_17, (1 - psi_1to5_t7_17)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_17, psi_1to5_t7_17);       //t8s5

    route_probability_from_KL_17[1] = temp1 * rep_vector(1, T);
    route_probability_from_KL_17[2] = temp2 * rep_vector(1, T);
    route_probability_from_KL_17[3] = temp3 * rep_vector(1, T);
    route_probability_from_KL_17[4] = temp4 * rep_vector(1, T);
    route_probability_from_KL_17[5] = temp5 * rep_vector(1, T);
    
    overall_survival_from_KL_17 = route_survival_from_KL_17[1] .* route_probability_from_KL_17[1] +
                                  route_survival_from_KL_17[2] .* route_probability_from_KL_17[2] +
                                  route_survival_from_KL_17[3] .* route_probability_from_KL_17[3] +
                                  route_survival_from_KL_17[4] .* route_probability_from_KL_17[4] + 
                                  route_survival_from_KL_17[5] .* route_probability_from_KL_17[5];
  }
  //2018 Gen Quants ----\\
  {
    matrix[2, T] temp_surv;
    matrix[2, T] temp1;
    matrix[2, T] temp2;
    matrix[2, T] temp3;
    matrix[2, T] temp4;
    matrix[2, T] temp5;

    temp1 = diag_pre_multiply(phi_s1_t0_18, alpha_s1_t0_18); //t1s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t1_18) * alpha_s1_t1_18; //t2s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, 1 - psi_1to2_t2_18) * alpha_s1_t2_18; //t3s1
    temp2 = diag_post_multiply(temp_surv, psi_1to2_t2_18 .* phi_s2_t3_18) *  alpha_s2_t3_18; //t4s2

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t3_18) * alpha_s1_t3_18; //t4s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t4_18) * alpha_s1_t4_18; //t5s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t5_18)  * alpha_s1_t5_18; //t6s1

    temp_surv = diag_post_multiply(temp1, phi_s1_t6_18);
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_18, ((1 - psi_1to3_t6_18) .* (1 - psi_1to4_t6_18))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_18, psi_1to3_t6_18); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_18, (1 - psi_1to3_t6_18) .* psi_1to4_t6_18); //t7s4

    temp_surv = diag_post_multiply(temp1, phi_s1_t7_18);
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_18, (1 - psi_1to5_t7_18)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_18, psi_1to5_t7_18);       //t8s5

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t8_18) * alpha_s1_t8_18; //t9s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s2_t8_18) * alpha_s2_t8_18; //t9s1 by yolo
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s3_t8_18) * alpha_s3_t8_18; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s4_t8_18) * alpha_s4_t8_18; //t9s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t8_18) * alpha_s5_t8_18; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t9_18) * alpha_s1_t9_18; //t10s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s1_t9_18) * alpha_s1_t9_18; //t10s1 by yolo    
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s1_t9_18) * alpha_s1_t9_18; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s1_t9_18) * alpha_s1_t9_18; //t10s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t9_18) * alpha_s5_t9_18; //t10s1 by georgiana

    prop_sacra_18    = (sum(temp1[1, ]) * N_rls_grp_18[1] + sum(temp1[2, ]) * N_rls_grp_18[2])/N_18;
    prop_yolo_18     = (sum(temp2[1, ]) * N_rls_grp_18[1] + sum(temp2[2, ]) * N_rls_grp_18[2])/N_18;
    prop_sutter_18   = (sum(temp3[1, ]) * N_rls_grp_18[1] + sum(temp3[2, ]) * N_rls_grp_18[2])/N_18;
    prop_steam_18    = (sum(temp4[1, ]) * N_rls_grp_18[1] + sum(temp4[2, ]) * N_rls_grp_18[2])/N_18;
    prop_georg_18    = (sum(temp5[1, ]) * N_rls_grp_18[1] + sum(temp5[2, ]) * N_rls_grp_18[2])/N_18;
    overall_survival_18 = prop_sacra_18 + prop_yolo_18 + prop_sutter_18 + prop_steam_18 + prop_georg_18;
    chipps_arrival_dist_18 = (((temp1[1, ] + temp2[1, ] + temp3[1, ] + temp4[1, ] + temp5[1, ]) * N_rls_grp_18[1] +
                              (temp1[2, ] + temp2[2, ] + temp3[2, ] + temp4[2, ] + temp5[2, ]) * N_rls_grp_18[2])/ N_18) / overall_survival_18;
  }

  {
    matrix[T, T] temp_surv;
    matrix[T, T] temp1;
    matrix[T, T] temp2;
    matrix[T, T] temp3;
    matrix[T, T] temp4;
    matrix[T, T] temp5;

    temp1 = diag_pre_multiply(phi_s1_t1_18, alpha_s1_t1_18); //t2s1

    temp_surv = temp1;
    temp1 = temp_surv * alpha_s1_t2_18; //t3s1
    temp2 = diag_post_multiply(temp_surv, phi_s2_t3_18) *  alpha_s2_t3_18; //t4s2

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t3_18) * alpha_s1_t3_18; //t4s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t4_18) * alpha_s1_t4_18; //t5s1

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t5_18)  * alpha_s1_t5_18; //t6s1

    temp_surv = diag_post_multiply(temp1, phi_s1_t6_18) * alpha_s1_t6_18;
    temp1 = temp_surv; //t7s1
    temp3 = temp_surv; //t7s3
    temp4 = temp_surv; //t7s4

    temp_surv = diag_post_multiply(temp1, phi_s1_t7_18) * alpha_s1_t7_18;
    temp1 = temp_surv; //t8s1
    temp5 = temp_surv;       //t8s5

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t8_18) * alpha_s1_t8_18; //t9s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s2_t8_18) * alpha_s2_t8_18; //t9s1 by yolo    
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s3_t8_18) * alpha_s3_t8_18; //t9s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s4_t8_18) * alpha_s4_t8_18; //t9s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t8_18) * alpha_s5_t8_18; //t9s5 by georgiana

    temp_surv = temp1;
    temp1 = diag_post_multiply(temp_surv, phi_s1_t9_18) * alpha_s1_t9_18; //t10s1 by sacramento
    temp_surv = temp2;
    temp2 = diag_post_multiply(temp_surv, phi_s1_t9_18) * alpha_s1_t9_18; //t10s1 by yolo    
    temp_surv = temp3;
    temp3 = diag_post_multiply(temp_surv, phi_s1_t9_18) * alpha_s1_t9_18; //t10s1 by sutter
    temp_surv = temp4;
    temp4 = diag_post_multiply(temp_surv, phi_s1_t9_18) * alpha_s1_t9_18; //t10s1 by steam
    temp_surv = temp5;
    temp5 = diag_post_multiply(temp_surv, phi_s5_t9_18) * alpha_s5_t9_18; //t10s1 by georgiana

    route_survival_from_KL_18[1] = temp1 * rep_vector(1, T);
    route_survival_from_KL_18[2] = temp2 * rep_vector(1, T);
    route_survival_from_KL_18[3] = temp3 * rep_vector(1, T);
    route_survival_from_KL_18[4] = temp4 * rep_vector(1, T);
    route_survival_from_KL_18[5] = temp5 * rep_vector(1, T);

    temp_surv = alpha_s1_t1_18;
    temp1 = diag_post_multiply(temp_surv, (1 - psi_1to2_t2_18));
    temp2 = diag_post_multiply(temp_surv, psi_1to2_t2_18);

    temp_surv = temp1 * alpha_s1_t2_18 * alpha_s1_t3_18 * alpha_s1_t5_18;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t6_18, ((1 - psi_1to3_t6_18) .* (1 - psi_1to4_t6_18))); //t7s1
    temp3 = temp_surv * diag_post_multiply(alpha_s1_t6_18, psi_1to3_t6_18); //t7s3
    temp4 = temp_surv * diag_post_multiply(alpha_s1_t6_18, (1 - psi_1to3_t6_18) .* psi_1to4_t6_18); //t7s4

    temp_surv = temp1;
    temp1 = temp_surv * diag_post_multiply(alpha_s1_t7_18, (1 - psi_1to5_t7_18)); //t8s1
    temp5 = temp_surv * diag_post_multiply(alpha_s1_t7_18, psi_1to5_t7_18);       //t8s5

    route_probability_from_KL_18[1] = temp1 * rep_vector(1, T);
    route_probability_from_KL_18[2] = temp2 * rep_vector(1, T);
    route_probability_from_KL_18[3] = temp3 * rep_vector(1, T);
    route_probability_from_KL_18[4] = temp4 * rep_vector(1, T);
    route_probability_from_KL_18[5] = temp5 * rep_vector(1, T);
    
    overall_survival_from_KL_18 = route_survival_from_KL_18[1] .* route_probability_from_KL_18[1] + 
                                  route_survival_from_KL_18[2] .* route_probability_from_KL_18[2] +
                                  route_survival_from_KL_18[3] .* route_probability_from_KL_18[3] +
                                  route_survival_from_KL_18[4] .* route_probability_from_KL_18[4] + 
                                  route_survival_from_KL_18[5] .* route_probability_from_KL_18[5];
  }

}
