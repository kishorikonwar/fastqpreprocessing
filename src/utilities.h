#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__
#include "fastqprocess.h"
#include <math.h>


typedef struct _input_options {
  vector<string> I1s, R1s, R2s;
  string white_list_file;
  unsigned int barcode_length, umi_length; 
  double bam_size;
} INPUT_OPTIONS;


int process_inputs(const INPUT_OPTIONS &,  const WHITE_LIST_DATA *) ;
long get_num_blocks(const INPUT_OPTIONS &options);

WHITE_LIST_DATA *read_white_list(const string &white_list_file);

void read_options(int, char **, INPUT_OPTIONS &); 
long filesize(const char *filename);
long getFileSize(const std::string &fileName);




#endif
