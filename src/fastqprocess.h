#ifndef __FASTQ_PROCESS_H__
#define __FASTQ_PROCESS_H__

#include <FastQFile.h>
#include "FastQStatus.h"
#include "BaseAsciiMap.h"
#include "SamFile.h"
#include "SamValidation.h"

#include <semaphore.h>
#include <thread>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <vector>
#include <functional>
#include <mutex>         


using namespace std;
typedef std::pair<string, bool>  STRING_BOOL_PAIR;

typedef vector<string>  STRING_VECTOR;

typedef std::unordered_map <std::string, int> STRING_BOOL_MAP;

typedef struct _white_list_data {
    STRING_BOOL_MAP mutations;
    STRING_VECTOR barcodes;
} WHITE_LIST_DATA;
 
typedef struct SamRecordBins {
    SamRecord **samrecords;
    int *num_records;
    
    vector<int> **file_index;

    int block_size;
    short int num_files;

    bool stop;
    short int num_threads;
    int active_thread_no;
} SAM_RECORD_BINS;

typedef struct _input_options {
  vector<string> I1s, R1s, R2s;
  string white_list_file;

} INPUT_OPTIONS;

void process_file(int tindex, String filename, String filename1, String filename2, \
                      const WHITE_LIST_DATA *, SAM_RECORD_BINS *) ;

int process_inputs(const INPUT_OPTIONS &,  const WHITE_LIST_DATA *) ;
void write_to_bam(int , SAM_RECORD_BINS *) ;



#endif
