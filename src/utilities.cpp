#include "utilities.h"

using namespace std;

long filesize(const char *filename)
{
  FILE *f = fopen(filename,"rb");  /* open the file in read only */

  long size = 0;
  if (fseek(f,0,SEEK_END)==0) /* seek was successful */
      size = ftell(f);
  fclose(f);
  return size;
}

long getFileSize(const std::string &fileName)
{
    ifstream file(fileName, ifstream::in | ifstream::binary);

    if(!file.is_open())
    {
        return -1;
    }

    file.seekg(0, ios::end);
    long fileSize = file.tellg();
    file.close();

    return fileSize;
}

WHITE_LIST_DATA *read_white_list(const string &white_list_file) {
   // be careful of caps 
   char ATCG[] = {'A', 'C', 'G', 'T', 'N'};

   fstream newfile;
   WHITE_LIST_DATA *white_list_data = new  WHITE_LIST_DATA;

   newfile.open(white_list_file, ios::in); //open a file to perform read operation using file object
   int k = 0;
   if (newfile.is_open()){   //checking whether the file is open
      string tp;
      while(getline(newfile, tp)){ //read data from file object and put it into string.
       //  cout << tp << "\n"; //print the data of the string
         white_list_data->barcodes.push_back(tp);

         for(int i=0; i < tp.size(); i++) {
           for(int j=0; j < 5; j++) {
             char c = tp[i];
             tp[i] = ATCG[j];
             if(white_list_data->mutations.find(tp) == white_list_data->mutations.end()) {
               white_list_data->mutations.insert( {tp, k} );
             } else {
               white_list_data->mutations[tp] = k;
             }
             tp[i] = c;
           }
         }
         //printf("%s\n", tp.c_str());
         white_list_data->mutations.at(tp) = -1;

         if( k%100000 == 0 && k != 0 ) printf("%d\n", k);
         k++;

      }
      newfile.close(); //close the file object.
   }

   std::cout << "Size of whitelist " << white_list_data->mutations.size() << std::endl;
   return white_list_data;
} 

void read_options(int argc, char **argv, INPUT_OPTIONS &options) {
  int c;
  int i;

  int verbose_flag;

  while (true) {
      static struct option long_options[] =
        {
          /* These options set a flag. */
          {"verbose", no_argument,  &verbose_flag, 1},
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"barcode-length",    required_argument, 0, 'b'},
          {"umi-length",    required_argument, 0, 'u'},
          {"bam-size",    required_argument, 0, 'B'},
          {"sample-id",    required_argument, 0, 's'},
          {"I1",  required_argument,  0, 'I'},
          {"R1",  required_argument,  0, 'R'},
          {"R2",    required_argument, 0, 'r'},
          {"white-list",    required_argument, 0, 'w'},
          {0, 0, 0, 0}
        };

      const char *help_messages[] = {
           "verbose messages  ", 
           "barcode length [required]", 
           "UMI length [required]",
           "output BAM file in GB [optional: default 1 GB]", 
           "sample id [required]", 
           "I1", 
           "R1 [required]",
           "R2 [required]", 
           "white list (from cellranger) of correct barcodes [required]", 
      };


      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "b:u:B:s:I:R:r:w:", long_options, &option_index);

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
        case 'b':
            options.barcode_length = atoi(optarg); 
            break;
        case 'u':
            options.umi_length = atoi(optarg); 
            break;
        case 'B':
            options.bam_size = atof(optarg); 
            break;
        case 's':
            options.sample_id = string(optarg); 
            break;
        case 'I':
            options.I1s.push_back(string(optarg)); 
            break;
        case 'R':
            options.R1s.push_back(string(optarg)); 
            break;
        case 'r':
            options.R2s.push_back(string(optarg)); 
            break;
        case 'w':
            options.white_list_file = string(optarg); 
            break;

        case '?':
        case 'h':
          i = 0;
          printf("Usage: %s [options] \n", argv[0]);
          while(long_options[i].name != 0) {
             printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name, \
                     long_options[i].has_arg==no_argument? "no argument" : "required_argument", \
                     help_messages[i]
                   );
             i = i + 1;
          }
          /* getopt_long already printed an error message. */
          return;
          break;
        default:
          abort ();
        }
    }

  if( options.R1s.size() != options.R2s.size() || options.R1s.size() ==0  )  {
     printf("R1 and R2 files mismatch i input\n");
     exit(0);
  }

  if( options.bam_size <= 0 )  {
     printf("Size of a bam file (in GB) cannot be negative\n");
     exit(0);
  }

  if( options.sample_id.size() == 0 )  {
     printf("Must provide a sample id or sampe name\n");
     exit(0);
  }

  if( options.barcode_length <= 0 )  {
     printf("Barcode length must be a positive integer\n");
     exit(0);
  }

  if( options.umi_length <= 0 )  {
     printf("UMI length must be a positive integer\n");
     exit(0);
  }

  if (verbose_flag) {
       std::cout << "I1/R1/R2 files" << std::endl;
       for(int i= 0; i < options.I1s.size(); i++) {
           std::cout << "\t" << options.I1s[i] << " " << filesize(options.I1s[i].c_str()) \
                     << " "  << options.R1s[i] << " " << filesize(options.R1s[i].c_str()) \
                     << " "  << options.R2s[i] << " " << filesize(options.R2s[i].c_str()) \
                     <<  std::endl;
       }
  }

}

long get_num_blocks(const INPUT_OPTIONS &options) {

    double tot_size = 0;
    for(int i= 0; i < options.R1s.size(); i++) {
        tot_size +=  filesize(options.I1s[i].c_str());
        tot_size +=  filesize(options.R1s[i].c_str());
        tot_size +=  filesize(options.R2s[i].c_str());
    }

    //printf("ceil %f\n",  tot_size/(1024*1024*1024)/(double)options.bam_size);
    return ceil( ( tot_size/(1024*1024*1024))/(double)options.bam_size);
}
 
