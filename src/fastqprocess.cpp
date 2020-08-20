#include "fastqprocess.h"
#include "utilities.h"

std::mutex mtx; 
sem_t *semaphores = 0;
sem_t *semaphores_workers = 0;

using namespace std;


void process_file(int tindex, String filename, String filename1, String filename2, \
                      const WHITE_LIST_DATA *, SAM_RECORD_BINS *) ;

int process_inputs(const STRING_VECTOR &, const STRING_VECTOR &, const STRING_VECTOR &, \
                   const WHITE_LIST_DATA *) ;
void write_to_bam(int , SAM_RECORD_BINS *) ;


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

  if( R1s.size() != R2s.size() || R1s.size() ==0  )  {
     printf("R1 and R2 files mismatch i input\n");
     exit(0);
  }


  if (verbose_flag) {
       std::cout << "I1/R1/R2 files" << std::endl;
       for(int i= 0; i < I1s.size(); i++) {
           std::cout << "\t" << I1s[i] << " " << filesize(I1s[i].c_str()) << " " << R1s[i] << " " << R2s[i] <<  std::endl;
       }
  }

  WHITE_LIST_DATA *white_list_data = read_write_list(white_list_file);

  process_inputs(I1s, R1s, R2s, white_list_data);
}

SAM_RECORD_BINS * create_samrecord_holders(short int nthreads, int block_size, \
                                           short int num_files) {

   SAM_RECORD_BINS *samrecord_data = new SAM_RECORD_BINS;

   if( (samrecord_data->samrecords = new SamRecord *[nthreads]) == 0 ) {
      std::cerr << "Failed to allocate memory for the samRecords pointer arrays" << std::endl ;
      return 0;
   }

   for(int i = 0; i < nthreads; i++) {
      if( (samrecord_data->samrecords[i] = new SamRecord[block_size]) == 0 ) {
         std::cerr << "Failed to allocate memory for the samRecords" << std::endl ;
         return 0;
        }
   }
   if( (samrecord_data->num_records = new int[nthreads]) == 0 ) {
      std::cerr << "Failed to allocate memory for the num records array" << std::endl ;
      return 0;
   }

   if( (samrecord_data->file_index = new vector<int> *[nthreads]) == 0 ) {
      std::cerr << "Failed to allocate memory for the pointer for  array of vectors" << std::endl ;
      return 0;
   }
   for(int i = 0; i < nthreads; i++) {
     if( (samrecord_data->file_index[i] = new vector<int>[num_files]) == 0 ) {
        std::cerr << "Failed to allocate memory for the vectors for index of file" << std::endl ;
        return 0;
     }
   }

   samrecord_data->block_size = block_size;
   samrecord_data->num_files = num_files;
   samrecord_data->stop = false;

   return samrecord_data;
}

int process_inputs(const STRING_VECTOR &I1s, const STRING_VECTOR &R1s, \
                   const  STRING_VECTOR &R2s, const WHITE_LIST_DATA *white_list_data) {

   int num_files = 200;
   int block_size = 100000;


   // create the data for the threads 
   SAM_RECORD_BINS *samrecord_data = create_samrecord_holders(R1s.size(), block_size, num_files);

     
   semaphores_workers = new sem_t[num_files];
   for(int i = 0; i < num_files; i++) {
      sem_init((semaphores_workers + i), 0, 0);
   }

   // create the bam file writers semaphores
   semaphores = new sem_t[num_files];
   for(int i = 0; i < num_files; i++) {
      sem_init((semaphores + i), 0, 0);
   }

   // execute the bam file writers threads
   std::thread *writers = new std::thread[num_files];
   for (int i = 0; i < num_files; i++)
   {
      writers[i] = std::thread(write_to_bam, i, samrecord_data);
   }

    
   // execute the fastq readers threads
   std::thread  *readers = new std::thread[R1s.size()];
   for (int i = 0; i < R1s.size(); i++) {
      readers[i] = std::thread(process_file, i, I1s[i].c_str(), R1s[i].c_str(), R2s[i].c_str(), \
                     white_list_data, samrecord_data);
   }

   // everyone comes home.
   for (int i = 0; i < R1s.size(); i++) {
      readers[i].join();
   }

   //set the stop flag for the writers
   samrecord_data->stop = true;

   //ask the writers to make one more loop in the whilte loop
   for(int j = 0; j < samrecord_data->num_files; j++) {
       sem_post(&semaphores[j]);
   }

   // wait for the writers to stop after they have seen the stop flag
   for (int i = 0; i < samrecord_data->num_files; i++) {
      writers[i].join();
   }

  //destroy the semaphores 
  for(int i = 0; i < samrecord_data->num_files; i++) {
         sem_destroy(&semaphores[i]);
  }

  // destroy the semaphores for semaphores_workers
  for(int i = 0; i < samrecord_data->num_files; i++) {
         sem_destroy(&semaphores_workers[i]);
  }


  delete [] readers;

  delete [] writers;
}

void write_to_bam(int windex, SAM_RECORD_BINS *samrecord_data) {
   std::thread::id  thread_id = std::this_thread::get_id();
     

   SamFile samOut;
   string outputfile;

   char buf[200];
   sprintf(buf, "checking-%ld.bam", thread_id);
   outputfile = buf;
   
   samOut.OpenForWrite(outputfile.c_str());
   // Write the sam header.

   SamFileHeader samHeader;
   samHeader.setHDTag("VN", "1.6");
   samHeader.setHDTag("SO", "unsorted");

   //samHeader.setRGTag("ID", "A");
  // samHeader.setRGTag("SM", "unknown");

   samOut.WriteHeader(samHeader);


   while(true) {
     sem_wait(&semaphores[windex]);
     
     SamRecord *samRecord  = samrecord_data->samrecords[samrecord_data->active_thread_no];

     //printf("\tWriting out : %d  :  active thread : %d\n", windex, samrecord_data->active_thread_no);
     for(auto index: samrecord_data->file_index[samrecord_data->active_thread_no][windex]) {
          samOut.WriteRecord(samHeader, samRecord[index]);
     }

     sem_post(&semaphores_workers[windex]);

     if(samrecord_data->stop) break;  
   }

   printf("closing %d\n", windex);
   samOut.Close();
}


void process_file(int tindex, String filename1, String filename2, String filename3, \
                   const WHITE_LIST_DATA* white_list_data, SAM_RECORD_BINS *samrecord_data) {
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

   std::thread::id  thread_id = std::this_thread::get_id();
     
   SamRecord *samRecord  = samrecord_data->samrecords[tindex];
   int block_size = samrecord_data->block_size;

   // Keep reading the file until there are no more fastq sequences to process.
   int i = 0;
   int n_barcode_errors = 0;
   int r =0;
   printf(" opening the thread in %d\n", tindex);
   while (fastQFile1.keepReadingFile())
   {
      if(fastQFile1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS && 
         fastQFile2.readFastQSequence() == FastQStatus::FASTQ_SUCCESS && 
         fastQFile3.readFastQSequence() == FastQStatus::FASTQ_SUCCESS 
      )
      {
         i = i + 1;
         //check the sequence names matching
         std::string a = std::string(fastQFile2.myRawSequence.c_str());
         std::string b = std::string(fastQFile2.myQualityString.c_str());

         string barcode = a.substr(0,16);
         string UMI  = a.substr(16,26);

         string barcodeQString = b.substr(0,16);
         string UMIQString  = b.substr(16,26);

         if( white_list_data->mutations.find(barcode) != white_list_data->mutations.end()) {
            //printf("record %d\n", r);
            samRecord[r].resetRecord();
            samRecord[r].setReadName(fastQFile3.mySequenceIdentifier.c_str());
            samRecord[r].setSequence(fastQFile3.myRawSequence.c_str());
            samRecord[r].setQuality(fastQFile3.myQualityString.c_str());

            samRecord[r].addTag("CR", 'Z', barcode.c_str());
            samRecord[r].addTag("CY", 'Z', barcodeQString.c_str());

            samRecord[r].addTag("UB", 'Z', UMI.c_str());
            samRecord[r].addTag("UY", 'Z', UMIQString.c_str());

            std::string indexseq = std::string(fastQFile1.myRawSequence.c_str());
            std::string indexSeqQual = std::string(fastQFile1.myQualityString.c_str());
            samRecord[r].addTag("SR", 'Z', indexseq.c_str());
            samRecord[r].addTag("SY", 'Z', indexSeqQual.c_str());
            samRecord[r].setFlag(4);
            
            samrecord_data->num_records[tindex]++;
            int bucket = std::hash<std::string>{}(barcode.c_str()) % samrecord_data->num_files;
            samrecord_data->file_index[tindex][bucket].push_back(r);

            r = r + 1;
            if(r == block_size || !fastQFile1.keepReadingFile()) {
               mtx.lock();
               samrecord_data->active_thread_no = tindex;

         //      printf("Ready to write %d\n", tindex);
               for(int j = 0; j < samrecord_data->num_files; j++) {
                  sem_post(&semaphores[j]);
               }

               // there is where I wait white the writers are writing
               for(int j = 0; j < samrecord_data->num_files; j++) {
                  sem_wait(&semaphores_workers[j]);
               }

               // they are done writing 
               
               for(int j = 0; j < samrecord_data->num_files; j++) {
                     samrecord_data->file_index[tindex][j].clear();
               }
               r = 0;
               samrecord_data->num_records[tindex] = 0;
               mtx.unlock();
            }
           }
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

         }
         if( i== 5000000 ) {  
             for(int j = 0; j < samrecord_data->num_files; j++) {
                 sem_destroy(&semaphores[j]);
             }
            break;
         }
   }
   // Finished processing all of the sequences in the file.
   // Close the input file.
   fastQFile1.closeFile();
   fastQFile2.closeFile();
   fastQFile3.closeFile();
   printf("Total reads %d, errors %d %f\n", i, n_barcode_errors, n_barcode_errors/(double)i *100);

}

