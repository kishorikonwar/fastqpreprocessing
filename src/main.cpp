#include "FastQFile.h"
#include "FastQStatus.h"
#include "BaseAsciiMap.h"
#include "SamFile.h"
#include "SamValidation.h"
#include <thread>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <vector>

using namespace std;
typedef std::unordered_map<std::string, bool> STRING_BOOL_MAP;
typedef std::pair<string, bool>  STRING_BOOL_PAIR;
typedef vector<string>  STRING_VEC;


void process_file(String filename, String filename1, String filename2, const STRING_BOOL_MAP *) ;
int process_inputs(const STRING_VEC &, const STRING_VEC &, const STRING_VEC &, const  STRING_BOOL_MAP *) ;
STRING_BOOL_MAP *read_write_list(const string &white_list_file);

/* Flag set by ‘--verbose’. */
static int verbose_flag;

int main (int argc, char **argv)
{
  int c;
  int i;
  vector<string> I1s, R1s, R2s;
  string white_list_file;

  while (1)
    {
      static struct option long_options[] =
        {
          /* These options set a flag. */
          {"verbose", no_argument,  &verbose_flag, 1},
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"I1",  required_argument,  0, 'I'},
          {"R1",  required_argument,  0, 'R'},
          {"R2",    required_argument, 0, 'r'},
          {"white-list",    required_argument, 0, 'w'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "I:R:r:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;
        case 'I':
          I1s.push_back(string(optarg)); 
          break;
        case 'R':
          R1s.push_back(string(optarg)); 
          break;
        case 'w':
          white_list_file = string(optarg); 
          break;
        case 'r':
          R2s.push_back(string(optarg)); 
          break;

        case '?':
        case 'h':
          i = 0;
          printf("Usage %s  [options] \n", argv[0]);
          while(long_options[i].name != 0) {
             printf("\t--%s       %s\n", long_options[i].name, \
                 long_options[i].has_arg==no_argument? "no argument" : "required_argument" );
             i = i + 1;
          }
          /* getopt_long already printed an error message. */
          return 0;
          break;
        default:
          abort ();
        }
    }

  printf("hello\n");
  if( R1s.size() != R2s.size() || R1s.size() ==0  )  {
     printf("R1 and R2 files mismatch i input\n");
     exit(0);
  }

  if (verbose_flag) {
       std::cout << "I1/R1/R2 files" << std::endl;
       for(int i= 0; i < I1s.size(); i++) {
           std::cout << "\t" << I1s[i] << " " << R1s[i] << " " << R2s[i] <<  std::endl;
       }
  }

  STRING_BOOL_MAP *white_list = read_write_list(white_list_file);
  process_inputs(I1s, R1s, R2s, white_list);
}


STRING_BOOL_MAP *read_write_list(const string &white_list_file) {
   char ATCG[] = {'A', 'T', 'C', 'G'};

   fstream newfile;
   STRING_BOOL_MAP *whitelist = new  STRING_BOOL_MAP;

   newfile.open(white_list_file, ios::in); //open a file to perform read operation using file object
   int k = 0;
   if (newfile.is_open()){   //checking whether the file is open
      string tp;
      while(getline(newfile, tp)){ //read data from file object and put it into string.
       //  cout << tp << "\n"; //print the data of the string
         for(int i=0; i < tp.size(); i++) {
             for(int j=0; j < 4; j++) {
                 char c = tp[i];
                 tp[i] = ATCG[j];
                 whitelist->insert( {tp, true} );
                 tp[i] = c;
             }
         }
         if( k%100000 == 0 ) printf("%d\n", k);
         k++;
      }
      newfile.close(); //close the file object.
   }

   std::cout << "size of whitelist " << whitelist->size() << std::endl;
   return whitelist;

} 

int process_inputs(const STRING_VEC &I1s, const STRING_VEC &R1s, const  STRING_VEC &R2s, const STRING_BOOL_MAP *whitelist) {
   std::thread  *workers;

   workers = new std::thread[R1s.size()];

   for (int i = 0; i < R1s.size(); i++)
   {
      workers[i] = std::thread(process_file, I1s[i].c_str(), R1s[i].c_str(), R2s[i].c_str(), whitelist);
   }


   // everyone comes home.
   for (int i = 0; i < R1s.size(); i++) 
   {
      workers[i].join();
   }

  delete [] workers;


}

void process_file(String filename1, String filename2, String filename3, const STRING_BOOL_MAP * whitelist) {
   FastQFile fastQFile1(4, 4);
   FastQFile fastQFile2(4, 4);
   FastQFile fastQFile3(4, 4);

   if(fastQFile1.openFile(filename1, BaseAsciiMap::UNKNOWN) != FastQStatus::FASTQ_SUCCESS)
   {
      std::cerr << "Failed to open file: " <<  filename1.c_str();
      return;
   }

   if(fastQFile2.openFile(filename2, BaseAsciiMap::UNKNOWN) != FastQStatus::FASTQ_SUCCESS)
   {
      std::cerr << "Failed to open file: " <<  filename2.c_str();
      return;
   }

   if(fastQFile3.openFile(filename3, BaseAsciiMap::UNKNOWN) != FastQStatus::FASTQ_SUCCESS)
   {
      std::cerr << "Failed to open file: " <<  filename3.c_str();
      return;
   }

   SamFile samOut;
   SamFileHeader samHeader;

   std::thread::id  thread_id = std::this_thread::get_id();
     
   std::vector<SamRecord *> records;

   string outputfile;
   char buf[200];
   sprintf(buf, "checking-%ld.bam", thread_id);
   outputfile = buf;
   
   samOut.OpenForWrite(outputfile.c_str());
   // Write the sam header.

   samHeader.setHDTag("VN:1.6",	"SO:unsorted");
   samOut.WriteHeader(samHeader);

   // Keep reading the file until there are no more fastq sequences to process.
   SamRecord samRecord;
   int i = 0;
   int n_barcode_errors = 0;
   while (fastQFile1.keepReadingFile())
   {
     // samRecord = new SamRecord;
      // Read one sequence. This call will read all the lines for 
      // one sequence.
      /////////////////////////////////////////////////////////////////
      // NOTE: It is up to you if you want to process only for success:
      //    if(readFastQSequence() == FASTQ_SUCCESS)
      // or for FASTQ_SUCCESS and FASTQ_INVALID: 
      //    if(readFastQSequence() != FASTQ_FAILURE)
      // Do NOT try to process on a FASTQ_FAILURE
      /////////////////////////////////////////////////////////////////
      if(fastQFile1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS && 
         fastQFile2.readFastQSequence() == FastQStatus::FASTQ_SUCCESS && 
         fastQFile3.readFastQSequence() == FastQStatus::FASTQ_SUCCESS 
      )
      {
         i = i + 1;
         std::string a = std::string(fastQFile2.myRawSequence.c_str());
         std::string b = std::string(fastQFile2.myQualityString.c_str());

         string barcode = a.substr(0,16);
         string UMI  = a.substr(16,26);

         string barcodeQString = a.substr(0,16);
         string UMIQString  = a.substr(16,26);

         if( whitelist->find(barcode) != whitelist->end() ) {
          //   printf("ERROR: %d %s\n", i,  barcode.c_str());
            samRecord.resetRecord();
            samRecord.setReadName(fastQFile3.mySequenceIdentifier.c_str());
            samRecord.setSequence(fastQFile3.myRawSequence.c_str());
            samRecord.setQuality(fastQFile3.myQualityString.c_str());

            samRecord.addTag("CR", 'Z', barcode.c_str());
            samRecord.addTag("CY", 'Z', barcodeQString.c_str());

            samRecord.addTag("UB", 'Z', UMI.c_str());
            samRecord.addTag("UY", 'Z', UMIQString.c_str());

            std::string indexseq = std::string(fastQFile1.myRawSequence.c_str());
            std::string indexSeqQual = std::string(fastQFile1.myQualityString.c_str());
            samRecord.addTag("SR", 'Z', indexseq.c_str());
            samRecord.addTag("SY", 'Z', indexSeqQual.c_str());

         //   records.push_back(samRecord);
            samRecord.setFlag(4);
            samOut.WriteRecord(samHeader, samRecord);
         }
         else {
             n_barcode_errors += 1;
         }

         if(i % 1000000 == 0) {  
             printf("%d\n", i);
             std::string a = std::string(fastQFile2.myRawSequence.c_str());

             printf("%s\n", fastQFile1.mySequenceIdLine.c_str());
             printf("%s\n", fastQFile2.mySequenceIdLine.c_str());
             printf("%s\n", fastQFile3.mySequenceIdLine.c_str());
 //            printf("%s\n", fastQFile1.myRawSequence.c_str());
         }
         if( i== 10000000 ) break;
      }
   }
   // Finished processing all of the sequences in the file.
   // Close the input file.
   fastQFile1.closeFile();
   printf("Total reads %d, errors %d %f\n", i, n_barcode_errors, n_barcode_errors/(double)i *100);
}


/*
#include "SamFile.h"
#include "SamValidation.h"

int sam(int argc, char ** argv)
{
   // Open the input file for reading.
   SamFile samIn;
   samIn.OpenForRead(argv[1]);

   // Open the output file for writing.
   SamFile samOut;
   samOut.OpenForWrite(argv[2]);

   // Read the sam header.
   SamFileHeader samHeader;
   samIn.ReadHeader(samHeader);

   // Write the sam header.
   samOut.WriteHeader(samHeader);

   SamRecord samRecord;

    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;

    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file, so write it.
        samOut.WriteRecord(samHeader, samRecord);
    }

    std::cout << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;
    std::cout << "Number of records written = " << 
        samOut.GetCurrentRecordCount() << std::endl;

    // Return success since a failure would have thrown
    // an exception.
    return(returnStatus);
 }
*/
