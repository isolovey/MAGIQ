/******************************************************************************
preproc.cpp

Project: 4t_cv
    Conversion routine for 4T Binary data to NMR286 / Fitman text format
    Copyright (c) November 1996 - Robert Bartha

File description:
    Preprocessing functions.

Supports:
    4t_cv.cpp

Revision: May 18, 2010 by John Drozd
          added error bounds checking for pointer arithmetic in scale function
          and added debugging printf statements in scale function


******************************************************************************/

#include "prot.h"

int pre_process(int *fid, Preprocess *preprocess, Procpar_info *procpar_info,
float **out_data, float **scratch_data){
    
    int i, k;
    double temp2, scaled_point;
    
    // Check to see if scaling is required... -scale on command line must be set
    // If this is set, scale the data such that the first point of the real data
    // falls between 0 and 1    
    
    for (i=0; i<=*fid; i++){
        
        // Find the maximum value of the first 25 points of the time domain data
        
        temp2 = sqrt(((double)out_data[i][0]*(double)out_data[i][0])+
        ((double)out_data[i][1]*(double)out_data[i][1]));
        scaled_point = temp2;
        for (k=2; k<=48; k+=2){
            
            // For debugging
            // printf("point %d   real = %f   imag = %f  \n", i, out_data[i][0], out_data[i][1]);
            
            temp2 = sqrt(((double)out_data[i][k]*(double)out_data[i][k])+
            ((double)out_data[i][k+1]*(double)out_data[i][k+1]));
            if (temp2 > scaled_point) {
                scaled_point = temp2;
            }
        }
        
        if (preprocess[i].fid_scale) {
            
            /*temp2 = sqrt(((double)out_data[i][0]*(double)out_data[i][0])+
             ((double)out_data[i][1]*(double)out_data[i][1]));
             scaled_point = temp2;*/
            if (scaled_point < 1) {
                while (scaled_point < 1) {
                    scaled_point = scaled_point * 10;
                    preprocess[i].scale_factor = preprocess[i].scale_factor * 10;
                }
            } else if (scaled_point > 10){
                while (scaled_point > 10) {
                    scaled_point = scaled_point / 10;
                    preprocess[i].scale_factor = preprocess[i].scale_factor / 10;
                }
            }
            
        }
        // always runs scale, if no scaling is required, scaling constant should = 1
        
        scale(out_data[i], &procpar_info[0], &preprocess[i]);
        
        // execute baseline correction if requested
        
        if (preprocess[i].bc) {
            baseline_correct(out_data[i], &procpar_info[i]);
        }
        
        // execute normalization if requested
        
        if (preprocess[i].max_normalize) {
            normalize(out_data[i], scratch_data[i], &procpar_info[i]);
        }
        
    }
    
    // execute eddy current correction if requested
    
    if (preprocess[0].pre_ecc) {
        
        ecc_correction(out_data[0], out_data[1], &procpar_info[0], &preprocess[0]);
    }
    
    // execute QUALITY deconvolution if requested
    if (preprocess[0].pre_quality) {
        
        quality(out_data[0], out_data[1], scratch_data[1], &procpar_info[0], &preprocess[0]);
    }
    // execute combined QUALITY deconvolution and ECC correction if requested
    
    if (preprocess[0].pre_quecc) {
        
        quecc(out_data[0], out_data[1], scratch_data[1], &procpar_info[0], &preprocess[0]);
    }
    
    //  if baseline correction is requested following ECC, QUALITY, or QUECC
    //  this is done under the guise of "-tilt"
    
    for (i=0; i<=*fid; i++){
        
        if (preprocess[i].tilt) {
            
            baseline_correct(out_data[i], &procpar_info[i]);
        }
    }
    
    //  apply filter if requested
    
    
    if (preprocess[0].comp_filter || preprocess[1].comp_filter ||
    preprocess[0].pre_quecc_if == YES || preprocess[1].pre_quecc_if == YES ) {
        
        filter(out_data[0], out_data[1], &procpar_info[0], &preprocess[0]);
    }
    
    // zero fill files if requested
    
    for (i=0; i<=*fid; i++){
        if (preprocess[i].data_zero_fill){
            zero_fill(out_data[0], out_data[1], &procpar_info[0], &preprocess[0]);
        }
    }
    
    return 1;
}




int scale(float *data, Procpar_info *procpar_info, Preprocess *preprocess) {
    
    int j=0;
    float next_point, half=(float)0.5;
    
    for (j=0; j<(procpar_info->num_points); j+=2) {

      /*
        printf("j = %d\n",j);
        printf("*(data+j) = %g\n",*(data+j));
        printf("*(data+j+1) = %g\n",*(data+j+1));
        printf("preprocess->scale_factor = %g\n",preprocess->scale_factor);
      */  

        *(data+j) = *(data+j) * preprocess->scale_factor;
        *(data+j+1) = *(data+j+1) * preprocess->scale_factor;
        
        if(fabs(*(data+j)) < 0.000000001){
           
	  /*
             printf("fabs(*(data+j)) < 0.000000001 \n");
             printf("*(data+j+2) = %g\n",*(data+j+2));            
	  */
             next_point = (float)fabs(*(data+j+2))*preprocess->scale_factor;
	  /*
             printf("next_point = %g\n", next_point);

             if(j>=2){
             printf("*(data+j-2) = %g\n",*(data+j-2));
             }
	  */

            //if(!(fabs(*(data+j-2)) < 0.000000001) && !(next_point < 0.000000001) && ((j+2)<(procpar_info->num_points))){
   
            if(j>=2){ //new line                          
 
            if(!(fabs(*(data+j-2)) < 0.000000001) && !(next_point < 0.000000001) && ((j+2)<(procpar_info->num_points))){            
                
	      /*    

                printf("! branch\n");

                printf("*(data+j-2) = %g\n", *(data+j-2));
 
	      */
 
                *(data+j) = (float)sqrt(0.5*((*(data+j-2))*(*(data+j-2))+
                (next_point *next_point)));
            }else{

	      /*  

                printf("else of ! branch\n");

                printf("*(data+j-2) = %g\n", *(data+j-2));

	      */ 

                *(data+j) = *(data+j-2) * (float)0.5;
            }} //added brace
        }
       
        //printf("*(data+j+1) = %g\n",*(data+j+1));
 
        if(fabs(*(data+j+1)) < 0.000000001){
  
	  /*

            printf("fabs(*(data+j+1)) < 0.000000001 \n");

            printf("*(data+j+1+2) = %g\n",*(data+j+1+2));
  
	  */

            next_point = (float)fabs(*(data+j+1+2))*preprocess->scale_factor;
            
            /*

            printf("next_point = %g\n", next_point);
 
            if(j>=1){
            printf("*(data+j+1-2) = %g\n",*(data+j+1-2));
            }
 
            */

            if(j>=1){ //new line
            if(!(fabs(*(data+j+1-2)) < 0.000000001) && !(next_point < 0.000000001) && ((j+2)<(procpar_info->num_points))){
              
	      /*
  
                printf("second ! branch\n");

                printf("*(data+j+1) = %g\n", *(data+j+1));
                printf("*(data+j+1-2) = %g\n", *(data+j+1-2));

	      */

                *(data+j+1) = (float)sqrt((0.5)*((*(data+j+1-2))*(*(data+j+1-2))+
                (next_point *next_point)));
            }else{

	      /*

                printf("else of second ! branch\n");

                printf("*(data+j+1) = %g\n", *(data+j+1));
                printf("*(data+j+1-2) = %g\n", *(data+j+1-2));

	      */

                *(data+j+1) = *(data+j+1-2) * (float)0.5;
            }} //added brace
        }
    }
    
    return 1;
}



int baseline_correct(float *data, Procpar_info *procpar_info){
    
    int last_eighth=0, j=0;
    float average_offset_real=float(0.0);
    float average_offset_imag=float(0.0);
    
    last_eighth = procpar_info->num_points/8;
    
    for(j=(procpar_info->num_points-last_eighth); j<procpar_info->num_points;j+=2){
        
        average_offset_real = average_offset_real + *(data+j);
        average_offset_imag = average_offset_imag + *(data+j+1);
    }
    
    average_offset_real = average_offset_real/(float)(last_eighth/2.0);
    average_offset_imag = average_offset_imag/(float)(last_eighth/2.0);
    
    for(j=0; j<(procpar_info->num_points);j+=2){
        
        *(data+j) = *(data+j) - average_offset_real;
        *(data+j+1) = *(data+j+1) - average_offset_imag;
    }
    
    return 1;
}



int ecc_correction(float *sup_data, float *unsup_data, Procpar_info *procpar_info,
Preprocess *preprocess){
    
    int j=0;
    double mag_sup=0, mag_unsup=0;
    double phase_sup=0, phase_unsup=0, result_phase_cor_sup=0;
    double result_phase_cor_unsup=0;
    
    
    for(j=0; j<(procpar_info->num_points); j+=2){
        
        mag_sup = sqrt((double)*(sup_data+j)*(double)*(sup_data+j)
        +(double)*(sup_data+j+1)*(double)*(sup_data+j+1));
        
        phase_sup = atan2((double)*(sup_data+j+1), (double)*(sup_data+j));
        
        mag_unsup = sqrt((double)*(unsup_data+j)*(double)*(unsup_data+j)
        +(double)*(unsup_data+j+1)*(double)*(unsup_data+j+1));
        
        phase_unsup = atan2((double)*(unsup_data+j+1), (double)*(unsup_data+j));
        
        result_phase_cor_sup = phase_sup - phase_unsup;
        result_phase_cor_unsup = phase_unsup - phase_unsup;
        
        *(sup_data+j) = (float)mag_sup * (float)cos(result_phase_cor_sup);
        *(sup_data+j+1) = (float)mag_sup * (float)sin(result_phase_cor_sup);
        
        *(unsup_data+j) = (float)mag_unsup * (float)cos(result_phase_cor_unsup);
        *(unsup_data+j+1) = (float)mag_unsup * (float)sin(result_phase_cor_unsup);
    }
    
    return 1;
    
}

int zero_fill(float *sup_data, float *unsup_data, Procpar_info *procpar_info,
Preprocess *preprocess){
    
    int j=0;
    
    // if zero fill is less than number of points send a warning
    
    if (preprocess->data_zero_fill < procpar_info->num_points){
        printf("\n Zero fill less than number of points... zero fill ignored\n");
        return 0;
    }
    
    for(j=procpar_info->num_points; j<preprocess->data_zero_fill; j+=2){
        
        *(sup_data+j) = (float)0.0;
        *(sup_data+j+1) = (float)0.0;
    }
    
    return 1;
    
}

int normalize(float *data, float *scratch, Procpar_info *procpar_info){
    
    int i=0;
    double magnitude = 0.0;
    double max_magnitude = 0.0;
    
    // determine magnitude of largest point...this may be shifted slighly if the
    // data is collected as an echo...therefore the first 100 points are searched
    // to find the largest
    
    for(i=0; i<100; i+=2){
        
        magnitude = sqrt((double)*(data+i)*(double)*(data+i)
        +(double)*(data+i+1)*(double)*(data+i+1));
        
        if(magnitude > max_magnitude) max_magnitude = magnitude;
        
    }
    
    for(i=0; i<(procpar_info->num_points); i+=2){
        
        *(scratch+i) = *(data+i) / (float)max_magnitude;
        *(scratch+i+1) = *(data+i+1) / (float)max_magnitude;
        
    }
    
    return 1;
    
}

int quality(float *sup_data, float *unsup_data, float *scratch, Procpar_info *procpar_info,
Preprocess *preprocess){
    
    int i=0;
    double mag_sup=0, mag_unsup=0, div_sup_mag=0, div_unsup_mag=0;
    double phase_sup=0, phase_unsup=0, div_sup_phase=0, div_unsup_phase=0;
    double phase_scratch=0, mag_scratch=0;
    
    for(i=0; i<(procpar_info->num_points); i+=2){
        
        mag_sup = sqrt((double)*(sup_data+i)*(double)*(sup_data+i)
        +(double)*(sup_data+i+1)*(double)*(sup_data+i+1));
        
        phase_sup = atan2((double)*(sup_data+i+1), (double)*(sup_data+i));
        
        mag_unsup = sqrt((double)*(unsup_data+i)*(double)*(unsup_data+i)
        +(double)*(unsup_data+i+1)*(double)*(unsup_data+i+1));
        
        phase_unsup = atan2((double)*(unsup_data+i+1), (double)*(unsup_data+i));
        
        mag_scratch = sqrt((double)*(scratch+i)*(double)*(scratch+i)
        +(double)*(scratch+i+1)*(double)*(scratch+i+1));
        
        phase_scratch = atan2((double)*(scratch+i+1), (double)*(scratch+i));
        
        // Apply QUALITY to suppressed file using scratch (normalized water)
        
        div_sup_mag = mag_sup/mag_scratch;
        div_sup_phase = phase_sup - phase_scratch;
        
        *(sup_data+i) = (float)div_sup_mag * (float)cos(div_sup_phase);
        *(sup_data+i+1) = (float)div_sup_mag * (float)sin(div_sup_phase);
        
        // Apply QUALITY to unsuppressed file using scratch (normalized water)
        
        div_unsup_mag = mag_unsup/mag_scratch;
        div_unsup_phase = phase_unsup - phase_scratch;
        
        *(unsup_data+i) = (float)div_unsup_mag * (float)cos(div_unsup_phase);
        *(unsup_data+i+1) = (float)div_unsup_mag * (float)sin(div_unsup_phase);
    }
    
    return 1;
    
}

int filter(float *sup_data, float *unsup_data, Procpar_info *procpar_info,
Preprocess *preprocess){
    
    int i=0;
    float dwell, delay;
    float last_point_quality_unsup=(float)0, first_point_ecc_unsup=(float)0;
    
    dwell = procpar_info->acquision_time/(procpar_info->num_points/2);
    delay = preprocess[0].pre_delay_time/float(1e6);
    
    // If a combined QUALITY/ECC correction is being done, make sure that the
    // estimated filter will eliminate the hard edge at the junction between
    // QUALITIED points and ECCed points
    
    if (preprocess->pre_quecc && preprocess[0].pre_quecc_if == YES){
        last_point_quality_unsup = *(unsup_data + (preprocess->pre_quecc_points - 2));
        first_point_ecc_unsup = *(unsup_data + (preprocess->pre_quecc_points));
        
        preprocess[0].comp_filter = (float)((-1.0 * (float)log(first_point_ecc_unsup / last_point_quality_unsup))/
        (fabs(((preprocess->pre_quecc_points-2)/2)*dwell+delay)*PI));
    }
    
    if (preprocess->pre_quecc && preprocess[1].pre_quecc_if == YES){
        last_point_quality_unsup = *(unsup_data + (preprocess->pre_quecc_points - 2));
        first_point_ecc_unsup = *(unsup_data + (preprocess->pre_quecc_points));
        
        preprocess[1].comp_filter = (float)((-1.0 * (float)log(first_point_ecc_unsup / last_point_quality_unsup))/
        (fabs(((preprocess->pre_quecc_points-2)/2)*dwell+delay)*PI));
    }
    
    
    
    if( preprocess[0].comp_filter != 0){
        
        for(i=0; i<(procpar_info->num_points); i+=2){
            if (preprocess->pre_quecc && i<preprocess->pre_quecc_points){
                *(sup_data+i) = *(sup_data+i) * (float)exp(-1*PI*fabs((double)(i/2) * dwell + delay) * preprocess[0].comp_filter);
                *(sup_data+i+1) = *(sup_data+i+1) * (float)exp(-1*PI*fabs((double)(i/2) * dwell + delay) * preprocess[0].comp_filter);
            } else if (!preprocess->pre_quecc){
                *(sup_data+i) = *(sup_data+i) * (float)exp(-1*PI*fabs((double)(i/2) * dwell + delay) * preprocess[0].comp_filter);
                *(sup_data+i+1) = *(sup_data+i+1) * (float)exp(-1*PI*fabs((double)(i/2) * dwell + delay) * preprocess[0].comp_filter);
            }
        }
        
    }
    
    if( preprocess[1].comp_filter != 0){
        
        for(i=0; i<(procpar_info->num_points); i+=2){
            if (preprocess->pre_quecc && i<preprocess->pre_quecc_points){
                *(unsup_data+i) = *(unsup_data+i) * (float)exp(-1*PI*fabs((double)(i/2) * dwell + delay) * preprocess[1].comp_filter);
                *(unsup_data+i+1) = *(unsup_data+i+1) * (float)exp(-1*PI*fabs((double)(i/2) * dwell + delay) * preprocess[1].comp_filter);
            } else if (!preprocess->pre_quecc){
                *(unsup_data+i) = *(unsup_data+i) * (float)exp(-1*PI*fabs((double)(i/2) * dwell + delay) * preprocess[1].comp_filter);
                *(unsup_data+i+1) = *(unsup_data+i+1) * (float)exp(-1*PI*fabs((double)(i/2) * dwell + delay) * preprocess[1].comp_filter);
                
            }
        }
    }
    
    return 1;
    
}

int quecc(float *sup_data, float *unsup_data, float *scratch, Procpar_info *procpar_info,
Preprocess *preprocess){
    
    int i=0;
    double mag_sup=0, mag_unsup=0, div_sup_mag=0, div_unsup_mag=0;
    double phase_sup=0, phase_unsup=0, div_sup_phase=0, div_unsup_phase=0;
    double phase_scratch=0, mag_scratch=0;
    double result_phase_cor_sup=0, result_phase_cor_unsup=0;
    
    for(i=0; i<(procpar_info->num_points); i+=2){
        
        mag_sup = sqrt((double)*(sup_data+i)*(double)*(sup_data+i)
        +(double)*(sup_data+i+1)*(double)*(sup_data+i+1));
        
        phase_sup = atan2((double)*(sup_data+i+1), (double)*(sup_data+i));
        
        mag_unsup = sqrt((double)*(unsup_data+i)*(double)*(unsup_data+i)
        +(double)*(unsup_data+i+1)*(double)*(unsup_data+i+1));
        
        phase_unsup = atan2((double)*(unsup_data+i+1), (double)*(unsup_data+i));
        
        if (i < preprocess->pre_quecc_points){
            
            mag_scratch = sqrt((double)*(scratch+i)*(double)*(scratch+i)
            +(double)*(scratch+i+1)*(double)*(scratch+i+1));
            
            phase_scratch = atan2((double)*(scratch+i+1), (double)*(scratch+i));
            
            
            // Apply QUALITY to suppressed file using scratch (normalized water)
            
            // Check to make sure that no values blow up...INF
            if (mag_scratch < 0.000001){
                div_sup_mag = 0.0;
            } else {
                div_sup_mag = mag_sup/mag_scratch;
            }
            
            div_sup_phase = phase_sup - phase_scratch;
            
            *(sup_data+i) = (float)div_sup_mag * (float)cos(div_sup_phase);
            *(sup_data+i+1) = (float)div_sup_mag * (float)sin(div_sup_phase);
            
            // Apply QUALITY to unsuppressed file using scratch (normalized water)
            
            //Check to make sure that no values blow up...INF
            if (mag_scratch < 0.000001){
                div_unsup_mag = 0.0;
            } else {
                div_unsup_mag = mag_unsup/mag_scratch;
            }
            
            div_unsup_phase = phase_unsup - phase_scratch;
            
            *(unsup_data+i) = (float)div_unsup_mag * (float)cos(div_unsup_phase);
            *(unsup_data+i+1) = (float)div_unsup_mag * (float)sin(div_unsup_phase);
        } else {
            
            result_phase_cor_sup = phase_sup - phase_unsup;
            result_phase_cor_unsup = phase_unsup - phase_unsup;
            
            *(sup_data+i) = (float)mag_sup * (float)cos(result_phase_cor_sup);
            *(sup_data+i+1) = (float)mag_sup * (float)sin(result_phase_cor_sup);
            
            *(unsup_data+i) = (float)mag_unsup * (float)cos(result_phase_cor_unsup);
            *(unsup_data+i+1) = (float)mag_unsup * (float)sin(result_phase_cor_unsup);
        }
    }
    
    return 1;
}
