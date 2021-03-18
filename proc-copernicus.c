#define VERSION_NUMBER 1.6

// proc-copernicus.c
// -------------------------------------------------------------------
// Converts Chilbolton Raw NetCDF into CF-compliant Processed format.
// Applies Z calibration if specified.
//
// Version 1.6
//
// 22-10-2004 : EGP : Adapted from proc-acrobat
// 28-11-2006 : JCN : Include merging of coded and uncoded modes (no LDR)
// 27-06-2007 : JCN : Noise levels adapted for data post-10/05/07
//
// 20080517 : OTD : introduction of LDR and ZDR parameters
// 20080606 : JCN : adapted UNCODED_FIRST_GOOD from 32 to 37 (greater leakage from 13/03/08 inclusive)
//                  adapted CODED_FIRST_GOOD from 46 to 47    
// 20081023 : JCN : added ability to change range offset
// 20101004 : JCN : reduce UNCODED_FIRST_GOOD from 37 to 1 (investigate leakage)
// 20110727 : JCN : include ZDR processing

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>

#define MISSING_VALUE -999
#define CODED_FIRST_GOOD 48
#define UNCODED_FIRST_GOOD 28 
#define END_BAD 19

int handle_error(int status)
{
  printf("NetCDF Error %i\n",status);
  return(status);
}

int add_text_att(int ncid, const char *varname, const char *attname, const char *tp)
{
   int status;
   int varid;
   int length;
   
   status = nc_inq_varid(ncid, varname, &varid);
   length = strlen(tp);
   status = nc_put_att_text(ncid, varid, attname, length, tp);
	
   return (status);
}


int add_text_global(int ncid, const char *attname, const char *tp)
{
   int status;
   int length;
   
   length=strlen(tp);
   status = nc_put_att_text(ncid, NC_GLOBAL, attname, length, tp);

   return(status);
}

void despeckle(float *data, int nx, int ny)
{
  int neighbours, dx, dy, i, j;

  // Despeckle data
  printf("Despeckling...\n");

  for(i=1; i<nx-1; i++){
    for(j=1; j<ny-1; j++){
      if (data[j+i*ny] != MISSING_VALUE) {
      neighbours=0;
      for(dy=-1; dy<=1; dy++)
        for(dx=-1; dx<=1; dx++){
          if(data[(j+dy)+(i+dx)*ny] != MISSING_VALUE) neighbours++;
        }
      if(neighbours<3) data[j+i*ny] = MISSING_VALUE;
      }
    }
  }
}


int process_data(int ncid_in, int ncid_out, size_t range_length, size_t time_length, const char *varname, float offset, int first_good, int end_bad, char *snr_var, float snr_thresh, float *corr)
{
   float *data, *SNRdata;
   int i,j;
   int temp_id,snr_id;
   int status;
   
   printf("Processing %s...\n",varname);

   // Allocate data memory
   data = (float *)malloc(range_length*time_length*sizeof(float));
   SNRdata = (float *)malloc(range_length*time_length*sizeof(float));
   
   // Read input values
   status = nc_inq_varid(ncid_in, varname, &temp_id);
   if(status != NC_NOERR) handle_error(status);
   status = nc_get_var_float(ncid_in, temp_id, data);
   if(status != NC_NOERR) handle_error(status);
   status = nc_inq_varid(ncid_in, snr_var, &snr_id);
   if(status != NC_NOERR) handle_error(status);
   status = nc_get_var_float(ncid_in, snr_id, SNRdata);
   if(status != NC_NOERR) handle_error(status);
   status = nc_inq_varid(ncid_out, varname, &temp_id);
   if(status != NC_NOERR) handle_error(status);

   // Loop through rays
   for(i=0; i<time_length; i++)
     {
	// Add cal offset
	for(j=0; j<range_length; j++){
	  data[j+i*range_length] += (offset + corr[j]);
          if(SNRdata[j+i*range_length] < snr_thresh)
            data[j+i*range_length] = MISSING_VALUE;
        }
	
	// Mask off missing data
	for(j=0; j<first_good; j++)
	  data[j+i*range_length] = MISSING_VALUE;
	if(end_bad>0)
	  {
	     for(j=range_length-end_bad; j<range_length; j++)
	       data[j+i*range_length] = MISSING_VALUE;
	  }
	
     }
   
   despeckle(data, (int)time_length, (int)range_length);

   status = nc_put_var_float(ncid_out, temp_id, data);
   if(status != NC_NOERR) handle_error(status);
   
   // Free data memory
   free(data);
   free(SNRdata);

   return(0);
}

void copy_var(int ncid_in, int ncid_out, char *name)
{
  int temp_id,status;
  int ndims, natts;
  nc_type xtype;
  int dimids[NC_MAX_VAR_DIMS];
  int newid;
  int i;
  char buffer[255];

  status = nc_inq_varid(ncid_in, name, &temp_id);
  if(status != NC_NOERR) handle_error(status);
  status = nc_inq_var(ncid_in, temp_id, 0, &xtype, &ndims, dimids, &natts);
  if(status != NC_NOERR) handle_error(status);
  status = nc_def_var(ncid_out, name, xtype, ndims, dimids, &newid);
  if(status != NC_NOERR) handle_error(status);

  // Now copy the attributes
  for(i=0; i<natts; i++){
    status = nc_inq_attname(ncid_in, temp_id, i, buffer);
    status = nc_copy_att(ncid_in, temp_id, buffer, ncid_out, newid);
  }

  nc_enddef(ncid_out);
  // Now copy the data
  status = nc_copy_var(ncid_in, temp_id, ncid_out);
  if(status != NC_NOERR) handle_error(status);
  nc_redef(ncid_out);
}

void copy_var_def(int ncid_in, int ncid_out, char *name)
{
  int temp_id,status;
  int ndims, natts;
  nc_type xtype;
  int dimids[NC_MAX_VAR_DIMS];
  int newid;
  int i;
  char buffer[255];

  printf("copying variable definition for: %s\n", name);

  status = nc_inq_varid(ncid_in, name, &temp_id);
  if(status != NC_NOERR) handle_error(status);
  status = nc_inq_var(ncid_in, temp_id, 0, &xtype, &ndims, dimids, &natts);
  if(status != NC_NOERR) handle_error(status);
  status = nc_def_var(ncid_out, name, xtype, ndims, dimids, &newid);
  if(status != NC_NOERR) handle_error(status);

  // Now copy the attributes
  for(i=0; i<natts; i++){
    status = nc_inq_attname(ncid_in, temp_id, i, buffer);
    status = nc_copy_att(ncid_in, temp_id, buffer, ncid_out, newid);
  }
  printf("completed copying\n");
}

void copy_var_def_rename(int ncid_in, int ncid_out, char *name, char *new_name, int copy_atts)
{
  int temp_id,status;
  int ndims, natts;
  nc_type xtype;
  int dimids[NC_MAX_VAR_DIMS];
  int newid;
  int i;
  char buffer[255];

  status = nc_inq_varid(ncid_in, name, &temp_id);
  if(status != NC_NOERR) handle_error(status);
  status = nc_inq_var(ncid_in, temp_id, 0, &xtype, &ndims, dimids, &natts);
  if(status != NC_NOERR) handle_error(status);
  status = nc_def_var(ncid_out, new_name, xtype, ndims, dimids, &newid);
  if(status != NC_NOERR) handle_error(status);
  if (copy_atts==1) {
    // Now copy the attributes
    for(i=0; i<natts; i++){
      status = nc_inq_attname(ncid_in, temp_id, i, buffer);
      status = nc_copy_att(ncid_in, temp_id, buffer, ncid_out, newid);
    }
  }
}

void copy_var_data(int ncid_in, int ncid_out, char *name)
{
  int status,temp_id,temp_id_out;
  char string[100];
  nc_type xtype;
  int ndims, natts;
  int dimids[NC_MAX_VAR_DIMS];
  float *buffer;
  size_t nvalues,length;
  int i;

  status = nc_inq_varid(ncid_in, name, &temp_id);
  if(status != NC_NOERR) handle_error(status);

  status = nc_inq_var(ncid_in, temp_id, string, &xtype, &ndims, dimids, &natts);
  if(status != NC_NOERR) handle_error(status);

  if(xtype != NC_FLOAT){
    printf("** ERROR! copy_var_data only works for float variables.");
    nc_close(ncid_in); nc_close(ncid_out);
    exit(1);
  }

  nvalues=1;
  for(i=0; i<ndims; i++){
    status=nc_inq_dim(ncid_in,dimids[i],string,&length);
    if(status != NC_NOERR) handle_error(status);
    nvalues=nvalues*length;
  }
  buffer = (float *)malloc(nvalues*sizeof(float));
  status = nc_get_var_float(ncid_in, temp_id, buffer);
  if(status != NC_NOERR) handle_error(status);
  status = nc_inq_varid(ncid_out, name, &temp_id_out);
  if(status != NC_NOERR) handle_error(status);
  status = nc_enddef(ncid_out);
  status = nc_put_var_float(ncid_out, temp_id_out, buffer);
  if(status != NC_NOERR) handle_error(status);

  free(buffer);
}

int create_mask(int ncid_in, int ncid_out, size_t range_length, size_t time_length, float snr_thresh, int first_good)
{
  float *mask, *SNR_data, *SNR_data_coded;
  float maxval;
  int i,j,k,stind,endind,tmp_int;
  int temp_id;
  int status;
  float SNR_noise_diff = 10.0;
  float diff_thresh = 12.0;
  float fill_value;
  int lobe_range = 19;

  static char longname[] = "Mask for merging modes (1=Coded, 2=Uncoded)";
  static char shortname[] = "Mask (1=Coded, 2=Uncoded)";
  static char units[] = " ";

  fill_value= (float)MISSING_VALUE;

// Mask=0 (Null-range sidelobes), =1 (Coded), =2 (Uncoded)
  printf("Creating mask for merging modes...\n");

  status= nc_redef(ncid_out);
  if(status != NC_NOERR) handle_error(status);
  status = nc_inq_varid(ncid_out,"MASK",&temp_id);
  if(status != NC_NOERR) handle_error(status);
  status = nc_put_att_text(ncid_out,temp_id,"long_name",strlen(longname),longname);
  if(status != NC_NOERR) handle_error(status);
  status = nc_put_att_text(ncid_out,temp_id,"short_name",strlen(shortname),shortname);
  if(status != NC_NOERR) handle_error(status);
  status = nc_put_att_text(ncid_out,temp_id,"units",strlen(units),units);
  if(status != NC_NOERR) handle_error(status);
  status = nc_put_att_float(ncid_out,temp_id,"missing_value", NC_FLOAT, 1, &fill_value );
  if(status != NC_NOERR) handle_error(status);
  status = nc_put_att_float(ncid_out,temp_id,"_FillValue", NC_FLOAT, 1, &fill_value );
  if(status != NC_NOERR) handle_error(status);

  status = nc_enddef(ncid_out);
  if(status != NC_NOERR) handle_error(status);

//Allocate data memory
  mask = (float*)malloc(range_length*time_length*sizeof(float));
  SNR_data = (float*)malloc(range_length*time_length*sizeof(float));
  SNR_data_coded = (float*)malloc(range_length*time_length*sizeof(float));

//Read SNR values
  status = nc_inq_varid(ncid_in,"SNR_HC",&temp_id);
  if(status != NC_NOERR) handle_error(status);
  status = nc_get_var_float(ncid_in,temp_id,SNR_data);
  if(status != NC_NOERR) handle_error(status);
  status = nc_inq_varid(ncid_in,"SNR_HCP",&temp_id);
  if(status != NC_NOERR) handle_error(status);
  status = nc_get_var_float(ncid_in,temp_id,SNR_data_coded);
  if(status != NC_NOERR) handle_error(status);
  status = nc_inq_varid(ncid_out,"MASK",&temp_id);
  if(status != NC_NOERR) handle_error(status);
  
//Loop through rays
for (i=0;i<time_length;i++)
{
  for (j=0;j<range_length;j++) {
       mask[j+i*range_length]=MISSING_VALUE; 
       if ((SNR_data_coded[j+i*range_length] > snr_thresh) & (j >= first_good)) mask[j+i*range_length]=1.0; 
       if ((SNR_data[j+i*range_length] > 6.0) | (j < first_good)) mask[j+i*range_length]=2.0;
      
       if (mask[j+i*range_length]==1.0) {
          tmp_int=0;
          if (j-lobe_range>tmp_int) tmp_int=j-lobe_range;
          stind=tmp_int+i*range_length;
          tmp_int=range_length-1;
          if (j+lobe_range<range_length-1) tmp_int=j+lobe_range;
          endind=tmp_int+i*range_length;
          maxval=MISSING_VALUE;
          for (k=stind;k<endind;k++) {
             if ( maxval<SNR_data[k]) maxval=SNR_data[k];
          }
          if ((maxval+SNR_noise_diff-diff_thresh > SNR_data[j+i*range_length]) & (SNR_data_coded[j+i*range_length] < 6.0)) mask[j+i*range_length]=0.0;
        } 
  }
}

  status=nc_put_var_float(ncid_out,temp_id,mask);
  if(status != NC_NOERR) handle_error(status);

//Free data memory
free(mask);
free(SNR_data);
free(SNR_data_coded);

return(0);
}

int apply_mask(int ncid_out, size_t range_length, size_t time_length, char *var1_name, char *var2_name, char *var_out_name, float level)
{
  float *var1_data, *var2_data, *var_out_data, *mask;
  int i,j;
  int temp_id;
  int status;

//Allocate memory
  var1_data=(float*)malloc(range_length*time_length*sizeof(float));
  var2_data=(float*)malloc(range_length*time_length*sizeof(float));
  var_out_data=(float*)malloc(range_length*time_length*sizeof(float));
  mask=(float*)malloc(range_length*time_length*sizeof(float));

//Read data
  status=nc_inq_varid(ncid_out,var1_name,&temp_id);
  if(status != NC_NOERR)handle_error(status);
  status=nc_get_var_float(ncid_out,temp_id,var1_data);
  if(status != NC_NOERR)handle_error(status);
  status=nc_inq_varid(ncid_out,var2_name,&temp_id);
  if(status != NC_NOERR)handle_error(status);
  status=nc_get_var_float(ncid_out,temp_id,var2_data);
  if(status != NC_NOERR)handle_error(status);
  status=nc_inq_varid(ncid_out,"MASK",&temp_id);
  if(status != NC_NOERR)handle_error(status);
  status=nc_get_var_float(ncid_out,temp_id,mask);
  if(status != NC_NOERR)handle_error(status);
  status=nc_inq_varid(ncid_out,var_out_name,&temp_id);
  if(status != NC_NOERR)handle_error(status);

//Loop through rays
for (i=0;i<time_length;i++)
{   
  for (j=0;j<range_length;j++) {
      if((mask[j+i*range_length]>=level) & (var2_data[j+i*range_length]!=MISSING_VALUE)) var_out_data[j+i*range_length]=var2_data[j+i*range_length];
      else if(mask[j+i*range_length]==1.0) var_out_data[j+i*range_length]=var1_data[j+i*range_length];
      else var_out_data[j+i*range_length]=MISSING_VALUE; 
   }
}

  status=nc_put_var_float(ncid_out,temp_id,var_out_data);
  if(status != NC_NOERR)handle_error(status);

// Free memory
  free(var1_data);
  free(var2_data);
  free(var_out_data);
  free(mask);
  return(0);
}

int main(int argc, char *argv[])
{
  char ncFileName[256];
  char outFileName[256];
  int status;
  int ncid_out, ncid_in;
  int ncid = 0;
  int s_len;
  int len;
  nc_type xtype;
  int ndims, nvars, ngatts, unlimdimid;
  char buffer[2000];
  char *buffer2;
  int SNR_id;
  size_t range_length,time_length,attlen;
  int temp_id;
  float *data;
  int i,j;
  char *char_time;
  struct tm *ptr;
  time_t lt;  
  float new_ZED_HC_cal, new_ZED_HCP_cal;
  float zdr_cal_offset;
  float new_range_offset;
  struct stat file_info;

  float ZED_H_cal =-999; 
  float ZED_HC_cal=-999;
  float ZED_HCP_cal=-999;
  float range_offset=-999;
  float az_offset=0;
  float temp_float;
  float *ZED_corr, *NULL_corr;
   
  // Check command-line arguments
  if(argc<3){
    printf("Usage:\n");
    printf("proc-copernicus infile.nc outfile.nc [-Zcal n -ZPcal n -rangeoffset n]\n");
    return(1);
  }

  if(argc>3)
     {
	for(i=3; i<argc; i++)
	  {	
	     if(strcmp(argv[i],"-ZED_HC_cal")==0)
	       {
		  i++;
		  ZED_HC_cal=atof(argv[i]);
	       }
             else if(strcmp(argv[i],"-ZED_HCP_cal")==0)
               {
                  i++;
                  ZED_HCP_cal =atof(argv[i]);
               }
             else if(strcmp(argv[i],"-ZED_H_cal")==0)
               {
                  i++;
                  ZED_H_cal =atof(argv[i]);
               }
             else if(strcmp(argv[i],"-rangeoffset")==0)
               {
                  i++;
                  range_offset=atof(argv[i]);
               } 
             else if(strcmp(argv[i],"-azoffset")==0)
               {
                  i++;
                  az_offset=atof(argv[i]);
               }
	     else
	       {
		  printf("Unknown command line argument %s\n",argv[i]);
		  return(1);
	       }
	  }
     }
   
   
  strcpy(ncFileName,argv[1]);
  strcpy(outFileName,argv[2]);
   
  // Open the input file
  printf("Opening input file: %s\n",ncFileName);
  status = nc_open(ncFileName, NC_NOWRITE, &ncid_in);
  if (status != NC_NOERR) handle_error(status);

  // Read size of range and time dimensions
  status = nc_inq_dimid(ncid_in, "time", &temp_id);
  status = nc_inq_dimlen(ncid_in, temp_id, &time_length);
  status = nc_inq_dimid(ncid_in, "range", &temp_id);
  status = nc_inq_dimlen(ncid_in, temp_id, &range_length); 
  printf("%d rays x %d gates\n",time_length,range_length); 

  if(time_length<1){
    printf("No rays in file.\n");
    sprintf(buffer,"rm %s",outFileName);
    system(buffer);
    status=nc_close(ncid_in);
    exit(1);
  }

  // Open the output file
  printf("Creating output file: %s\n",outFileName);
  status = nc_create(outFileName, NC_CLOBBER, &ncid_out);  
  if (status != NC_NOERR) handle_error(status);

  // Go to define mode
  status = nc_redef(ncid_out);
    
  // Copy dimensions from input to output
  printf("Copying dimensions...\n");
  status = nc_inq(ncid_in, &ndims, &nvars, &ngatts, &unlimdimid);
  for(i=0; i<ndims; i++) {
    status = nc_inq_dim(ncid_in, i, buffer, &j);
    status = nc_def_dim(ncid_out, buffer, j, &temp_id);
  }

  // Copy global attributes
  printf("Copying global attributes...\n");
  for(i=0; i<ngatts; i++){
    status = nc_inq_attname(ncid_in, NC_GLOBAL, i, buffer);
    status = nc_copy_att(ncid_in, NC_GLOBAL, buffer, ncid_out, NC_GLOBAL);
  } 

  // Add new global attributes
  printf("Adding global attributes...\n");
  add_text_global(ncid_out, "Conventions", "CF-1.0");
  add_text_global(ncid_out, "comment", "Data restrictions: For academic research use only");
  add_text_global(ncid_out, "references", "Documentation may be found at: http://www.badc.nerc.ac.uk/data/chilbolton\nFor additional information please contact: chilbolton-netcdf@listserv.rl.ac.uk");
  add_text_global(ncid_out, "institution", "Original data recorded at Chilbolton Observatory, part of the Radio Communications Research Unit, CCLRC, UK: http://www.rcru.rl.ac.uk/chilbolton\nInstrument owned by the Department of Meteorology, University of Reading, UK: http://www.met.rdg.ac.uk/radar \nData recorded on behalf of the Universities Facility for Atmospheric Measurements (UFAM): http://www.env.leeds.ac.uk/ufam \nData held at the British Atmospheric Data Centre (BADC): http://www.badc.nerc.ac.uk");
  add_text_global(ncid_out, "title", "UFAM/Chilbolton 35-GHz Cloud Radar (Copernicus)"); 

  // Update history attribute
  printf("Updating history...\n");
  status = nc_get_att_text(ncid_in, NC_GLOBAL, "history", buffer);
  status = nc_inq_attlen(ncid_in, NC_GLOBAL, "history", &attlen);
  buffer[attlen] = '\n';
  buffer2 = &buffer[attlen+1];
  // Append new history string
  time(&lt);
  char_time=asctime(gmtime(&lt));
  strncpy(buffer2,char_time,strlen(char_time) - 1);
  buffer2[strlen(char_time) - 1] = '\0';
  strcat(buffer2," - ");
  for (i = 0; i < argc; i++ )
    {
      strcat(buffer2, argv[i]);
      strcat(buffer2, " ");
    }
  // Write history global attribute
  add_text_global(ncid_out, "history", buffer);  


  // Copy non-observable variables
  // Note: won't work for "unlimited" dimensions
  printf("Copying non-observable variables...\n");
  copy_var(ncid_in, ncid_out, "pulse_compression_code");
  copy_var(ncid_in, ncid_out, "file_state");
  copy_var(ncid_in, ncid_out, "clock");
  copy_var(ncid_in, ncid_out, "transmit_power");
  copy_var(ncid_in, ncid_out, "pulse_period");
  copy_var(ncid_in, ncid_out, "antenna_diameter");
  copy_var(ncid_in, ncid_out, "beamwidthV");
  copy_var(ncid_in, ncid_out, "beamwidthH");
  copy_var(ncid_in, ncid_out, "prf");
  copy_var(ncid_in, ncid_out, "frequency");
  copy_var(ncid_in, ncid_out, "height");
  copy_var(ncid_in, ncid_out, "longitude");
  copy_var(ncid_in, ncid_out, "latitude");
  copy_var(ncid_in, ncid_out, "range");

  // Copy time, elevation and azimuth variables
  // NB: This is different because they were unlimited in dimension
  printf("Copying time, elevation and azimuth variables...\n");
  copy_var_def(ncid_in, ncid_out, "time");
  copy_var_def(ncid_in, ncid_out, "elevation");
  copy_var_def(ncid_in, ncid_out, "azimuth");
  copy_var_data(ncid_in, ncid_out, "time");
  copy_var_data(ncid_in, ncid_out, "elevation");
  copy_var_data(ncid_in, ncid_out, "azimuth");

  status = nc_redef(ncid_out);

  copy_var_def(ncid_in, ncid_out, "SNR_HC");
  copy_var_def(ncid_in, ncid_out, "SNR_HCP");
  copy_var_data(ncid_in, ncid_out, "SNR_HC");
  copy_var_data(ncid_in, ncid_out, "SNR_HCP");
	
  // Enter define mode
  status = nc_redef(ncid_out);

  // Add CF standard names
  printf("Adding CF standard names...\n");
  add_text_att(ncid_out, "time", "standard_name", "time");
  add_text_att(ncid_out, "height", "standard_name", "height");
  add_text_att(ncid_out, "latitude", "standard_name", "latitude");
  add_text_att(ncid_out, "longitude", "standard_name", "longitude");


  printf("Fixing units...\n");
  // Make units CF-compliant
  add_text_att(ncid_out, "latitude", "units", "degree_east");
  add_text_att(ncid_out, "longitude", "units","degree_north");
  
  // Check for old-style time units
  status = nc_inq_varid(ncid_out, "time", &temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_get_att_text(ncid_out, temp_id, "long_name", buffer);
  if (status != NC_NOERR) handle_error(status);
  status = nc_inq_attlen(ncid_out, temp_id, "long_name", &attlen);
  if (status != NC_NOERR) handle_error(status);
  buffer[attlen]='\0';
  if(strcmp(buffer,"time")!=0){
    // Apply correct CF time units
    status = nc_put_att_text(ncid_out,temp_id,"units",attlen,buffer);
    if (status != NC_NOERR) handle_error(status);
    status = nc_put_att_text(ncid_out,temp_id,"long_name",4,"time");
    if (status != NC_NOERR) handle_error(status);
  }

  //Add ZDR calibration offset
  zdr_cal_offset=-3.0;
  status = nc_inq_varid(ncid_out,"ZDR_C",&temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_put_att_float(ncid_out,temp_id,"applied_calibration_offset", NC_FLOAT, 1, &zdr_cal_offset );
  if(status != NC_NOERR) handle_error(status);
  status = nc_inq_varid(ncid_out,"ZDR_CP",&temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_put_att_float(ncid_out,temp_id,"applied_calibration_offset", NC_FLOAT, 1, &zdr_cal_offset );
  if(status != NC_NOERR) handle_error(status);
  

  // Copy observable definitions
  copy_var_def(ncid_in, ncid_out, "ZED_HC");
  copy_var_def(ncid_in, ncid_out, "VEL_HC");
  copy_var_def(ncid_in, ncid_out, "SPW_HC");
  copy_var_def(ncid_in, ncid_out, "VEL_HCD");
  copy_var_def(ncid_in, ncid_out, "ZED_HCP");
  copy_var_def(ncid_in, ncid_out, "VEL_HCP");
  copy_var_def(ncid_in, ncid_out, "SPW_HCP");
  copy_var_def(ncid_in, ncid_out, "VEL_HCDP");
  copy_var_def(ncid_in, ncid_out, "LDR_HC");
  copy_var_def(ncid_in, ncid_out, "LDR_HCP");
  copy_var_def(ncid_in, ncid_out, "ZDR_C");
  copy_var_def(ncid_in, ncid_out, "ZDR_CP");

  // Create merged observable definitions
  copy_var_def_rename(ncid_in, ncid_out, "ZED_HC", "ZED_HCM",1);
  copy_var_def_rename(ncid_in, ncid_out, "VEL_HC", "VEL_HCM",1);
  copy_var_def_rename(ncid_in, ncid_out, "SPW_HC", "SPW_HCM",1);
  copy_var_def_rename(ncid_in, ncid_out, "VEL_HCD", "VEL_HCDM",1);
  copy_var_def_rename(ncid_in, ncid_out, "LDR_HC", "LDR_HCM",1);
  copy_var_def_rename(ncid_in, ncid_out, "ZDR_C", "ZDR_CM",1);
  copy_var_def_rename(ncid_in, ncid_out, "ZED_HC", "MASK",0);

  // Add attributes to aid in plotting
  add_text_att(ncid_out, "ZED_HC", "short_name", "Radar reflectivity factor");
  add_text_att(ncid_out, "ZED_HC", "units_html", "dBZ");
  add_text_att(ncid_out, "VEL_HC", "short_name", "Doppler velocity");
  add_text_att(ncid_out, "VEL_HC", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "VEL_HCD", "short_name", "Standard deviation of mean velocity");
  add_text_att(ncid_out, "VEL_HCD", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "SPW_HC", "short_name", "Spectral width");
  add_text_att(ncid_out, "SPW_HC", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "ZED_HCP", "short_name", "Radar reflectivity factor (coded)");
  add_text_att(ncid_out, "ZED_HCP", "units_html", "dBZ");
  add_text_att(ncid_out, "VEL_HCP", "short_name", "Doppler velocity (coded)");
  add_text_att(ncid_out, "VEL_HCP", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "VEL_HCDP", "short_name", "Standard deviation of mean velocity (coded)");
  add_text_att(ncid_out, "VEL_HCDP", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "SPW_HCP", "short_name", "Spectral width (coded)");
  add_text_att(ncid_out, "SPW_HCP", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "ZED_HCM", "short_name", "Radar reflectivity factor (merged)");
  add_text_att(ncid_out, "ZED_HCM", "units_html", "dBZ");
  add_text_att(ncid_out, "VEL_HCM", "short_name", "Doppler velocity (merged)");
  add_text_att(ncid_out, "VEL_HCM", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "VEL_HCDM", "short_name", "Standard deviation of mean velocity (merged)");
  add_text_att(ncid_out, "VEL_HCDM", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "SPW_HCM", "short_name", "Spectral width (merged)");
  add_text_att(ncid_out, "SPW_HCM", "units_html", "m s<sup>-1</sup>");
  add_text_att(ncid_out, "SNR_HC", "short_name", "SNR");
  add_text_att(ncid_out, "SNR_HC", "units_html", "dB");
  add_text_att(ncid_out, "SNR_HCP", "short_name", "SNR (coded)");
  add_text_att(ncid_out, "SNR_HCP", "units_html", "dB");
  add_text_att(ncid_out, "LDR_HC", "short_name", "LDR");
  add_text_att(ncid_out, "LDR_HC", "units_html", "dB");
  add_text_att(ncid_out, "LDR_HCP", "short_name", "LDR (coded)");
  add_text_att(ncid_out, "LDR_HCP", "units_html", "dB");
  add_text_att(ncid_out, "LDR_HCM", "short_name", "LDR (merged)");
  add_text_att(ncid_out, "LDR_HCM", "units_html", "dB");
  add_text_att(ncid_out, "ZDR_C", "short_name", "ZDR");
  add_text_att(ncid_out, "ZDR_C", "units_html", "dB");
  add_text_att(ncid_out, "ZDR_CP", "short_name", "ZDR (coded)");
  add_text_att(ncid_out, "ZDR_CP", "units_html", "dB");
  add_text_att(ncid_out, "ZDR_CM", "short_name", "ZDR (merged)");
  add_text_att(ncid_out, "ZDR_CM", "units_html", "dB");

  // Copy these observables if they exist
  if(nc_inq_varid(ncid_in,"ZED_HCD",&temp_id)==NC_NOERR){
    copy_var_def(ncid_in, ncid_out, "ZED_HCD");
    add_text_att(ncid_out, "ZED_HCD", "short_name", "Standard deviation of radar reflectivity");
    add_text_att(ncid_out, "ZED_HCD", "units_html", "dB");
    }
  if(nc_inq_varid(ncid_in,"ZED_HCDP",&temp_id)==NC_NOERR){
    copy_var_def(ncid_in, ncid_out, "ZED_HCDP");
    add_text_att(ncid_out, "ZED_HCDP", "short_name", "Standard deviation of radar reflectivity (coded)");
    add_text_att(ncid_out, "ZED_HCDP", "units_html", "dB");
    copy_var_def_rename(ncid_in, ncid_out, "ZED_HCD", "ZED_HCDM",1);
    add_text_att(ncid_out, "ZED_HCDM", "short_name", "Standard deviation of radar reflectivity (merged)");
    add_text_att(ncid_out, "ZED_HCDM", "units_html", "dB");
    }

  // Remove unwanted attributes
  status = nc_inq_varid(ncid_out,"ZED_HCM",&temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_del_att(ncid_out,temp_id,"applied_calibration_offset");
  if (status != NC_NOERR) handle_error(status);
  
  if(nc_inq_varid(ncid_out,"ZED_HCDM",&temp_id)==NC_NOERR){
     status = nc_inq_varid(ncid_out,"ZED_HCDM",&temp_id);
     if (status != NC_NOERR) handle_error(status);
     status = nc_del_att(ncid_out,temp_id,"applied_calibration_offset");
     if (status != NC_NOERR) handle_error(status);
  }

  // Come out of define mode
  status = nc_enddef(ncid_out);

  status = nc_inq_varid(ncid_in, "ZED_HC", &temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_get_att_float(ncid_in, temp_id, "applied_calibration_offset", &temp_float);
  if (status != NC_NOERR) handle_error(status);
  if(ZED_HC_cal==-999) ZED_HC_cal=temp_float;
  new_ZED_HC_cal=ZED_HC_cal-temp_float;

  status = nc_inq_varid(ncid_in, "ZED_HCP", &temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_get_att_float(ncid_in, temp_id, "applied_calibration_offset", &temp_float);
  if (status != NC_NOERR) handle_error(status);
  if(ZED_HCP_cal==-999) ZED_HCP_cal=temp_float;
  new_ZED_HCP_cal=ZED_HCP_cal-temp_float;
	printf("calibrations : %f %f\n", new_ZED_HC_cal,  new_ZED_HCP_cal);

  status = nc_inq_varid(ncid_in, "range", &temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_get_att_float(ncid_in, temp_id, "applied_range_offset", &temp_float);
  if (status != NC_NOERR) handle_error(status);
  if(range_offset==-999) range_offset=temp_float;
  new_range_offset=range_offset-temp_float;
  NULL_corr = (float *)malloc(range_length*sizeof(float));
  ZED_corr = (float *)malloc(range_length*sizeof(float));
  for (j=0;j<range_length; j++) {
     NULL_corr[j]=0.0;
     ZED_corr[j]=0.0;
   }

  if (new_range_offset != 0.0) {
     data = (float *)malloc(range_length*sizeof(float));
     status = nc_inq_varid(ncid_in,"range",&temp_id); 
     if (status != NC_NOERR) handle_error(status);
     status = nc_get_var_float(ncid_in,temp_id,data);
     if (status != NC_NOERR) handle_error(status);
     status = nc_inq_varid(ncid_out,"range",&temp_id); 
     if (status != NC_NOERR) handle_error(status);
     for (j=0;j<range_length;j++) {
              ZED_corr[j] = 20*log10(1+new_range_offset/data[j]); 
              data[j] += new_range_offset;
     }
     status = nc_put_var_float(ncid_out, temp_id, data);  
     if (status != NC_NOERR) handle_error(status);
     status = nc_put_att_float(ncid_out,temp_id,"applied_range_offset",NC_FLOAT,1,&range_offset);
     if (status != NC_NOERR) handle_error(status);
     free(data);
  }


  // Mark gates at beginning and end as missing
  // and add cal offset if needed
  process_data(ncid_in, ncid_out, range_length, time_length, "ZED_HC", new_ZED_HC_cal, UNCODED_FIRST_GOOD, END_BAD, "SNR_HC", 1.25, ZED_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "ZED_HCP", new_ZED_HCP_cal, CODED_FIRST_GOOD, END_BAD, "SNR_HCP", 0.85, ZED_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "VEL_HC", 0, UNCODED_FIRST_GOOD, END_BAD, "SNR_HC", 1.25, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "VEL_HCD", 0, UNCODED_FIRST_GOOD, END_BAD, "SNR_HC", 1.25, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "VEL_HCP", 0, CODED_FIRST_GOOD, END_BAD, "SNR_HCP", 0.85, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "VEL_HCDP", 0, CODED_FIRST_GOOD, END_BAD, "SNR_HCP", 0.85, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "SPW_HC", 0, UNCODED_FIRST_GOOD, END_BAD, "SNR_HC", 1.25, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "SPW_HCP", 0, CODED_FIRST_GOOD, END_BAD, "SNR_HCP", 0.85, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "LDR_HC", 0, UNCODED_FIRST_GOOD, END_BAD, "SNR_HC", 1.25, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "LDR_HCP", 0, CODED_FIRST_GOOD, END_BAD, "SNR_HCP", 0.85, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "ZDR_C", zdr_cal_offset, UNCODED_FIRST_GOOD, END_BAD, "SNR_HC", 1.25, NULL_corr);
  process_data(ncid_in, ncid_out, range_length, time_length, "ZDR_CP", zdr_cal_offset, CODED_FIRST_GOOD, END_BAD, "SNR_HCP", 0.85, NULL_corr);
  if(nc_inq_varid(ncid_in,"ZED_HCD",&temp_id)==NC_NOERR){
    process_data(ncid_in, ncid_out, range_length, time_length, "ZED_HCD", 0, UNCODED_FIRST_GOOD, END_BAD, "SNR_HC", 1.25, NULL_corr);
  }
  if(nc_inq_varid(ncid_in,"ZED_HCDP",&temp_id)==NC_NOERR){
    process_data(ncid_in, ncid_out, range_length, time_length, "ZED_HCDP", 0, CODED_FIRST_GOOD, END_BAD, "SNR_HCP", 0.85, NULL_corr);
  }
 
  free(ZED_corr);
  free(NULL_corr);  

  printf("Changing calibration offset attribute...\n");
  // Rectify calibration offset attribute
  status = nc_inq_varid(ncid_out,"ZED_HC",&temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_put_att_float(ncid_out, temp_id, "applied_calibration_offset", NC_FLOAT, 1, &ZED_HC_cal);
  if (status != NC_NOERR) handle_error(status);
  status = nc_inq_varid(ncid_out,"ZED_HCP",&temp_id);
  if (status != NC_NOERR) handle_error(status);
  status = nc_put_att_float(ncid_out, temp_id, "applied_calibration_offset", NC_FLOAT, 1, &ZED_HCP_cal);
  if (status != NC_NOERR) handle_error(status);

  // Create mask for merged data
  status =  create_mask(ncid_in, ncid_out, range_length, time_length, 0.85, CODED_FIRST_GOOD);
  if (status != NC_NOERR) handle_error(status);

  // Apply the mask
  status = apply_mask(ncid_out, range_length, time_length, "ZED_HCP", "ZED_HC", "ZED_HCM", 2.0);
  if (status != NC_NOERR) handle_error(status);
  status = apply_mask(ncid_out, range_length, time_length, "VEL_HCP", "VEL_HC", "VEL_HCM", 2.0);
  if (status != NC_NOERR) handle_error(status);
  status = apply_mask(ncid_out, range_length, time_length, "SPW_HCP", "SPW_HC", "SPW_HCM", 2.0);
  if (status != NC_NOERR) handle_error(status);
  status = apply_mask(ncid_out, range_length, time_length, "VEL_HCDP", "VEL_HCD", "VEL_HCDM", 2.0);
  if (status != NC_NOERR) handle_error(status);
  status = apply_mask(ncid_out, range_length, time_length, "ZED_HCDP", "ZED_HCD", "ZED_HCDM", 2.0);
  if (status != NC_NOERR) handle_error(status);
  status = apply_mask(ncid_out, range_length, time_length, "LDR_HCP", "LDR_HC", "LDR_HCM", 2.0);
  if (status != NC_NOERR) handle_error(status);
  status = apply_mask(ncid_out, range_length, time_length, "ZDR_CP", "ZDR_C", "ZDR_CM", 2.0);
  if (status != NC_NOERR) handle_error(status);



  // Close the output file
  printf("Closing Files.\n");
  status = nc_close(ncid_in);
  status = nc_close(ncid_out);
  printf("All done.\n");

  /* change the group and permissions of the file */
  /* first we have to found who actually owns it */
  status = stat( outFileName, &file_info);
  if (status != 0) {
        printf("an error while finding out information on the file has occurred\n");
        exit(1);
  }

  /* 701 is the gid of radar_dt */
  status = chown ( outFileName, file_info.st_uid, 701);
  if (status != 0) {
        printf("there has been an error while setting the gid\n");
        exit(1);
  }
  /* change permissions on group */
  status = chmod ( outFileName, S_IRUSR | S_IWUSR | S_IWGRP | S_IRGRP | S_IROTH  );
  if (status != 0) {
        printf("there has been an error while setting the group permission to write\n");
  }


  return(0);
}

