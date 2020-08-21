#include "fastqprocess.h"
#include "utilities.h"

using namespace std;

/* Flag set by ‘--verbose’. */
static int verbose_flag;

int main (int argc, char **argv)
{
  int c;
  int i;

  INPUT_OPTIONS options;

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

  if( options.R1s.size() != options.R2s.size() || options.R1s.size() ==0  )  {
     printf("R1 and R2 files mismatch i input\n");
     exit(0);
  }


  if (verbose_flag) {
       std::cout << "I1/R1/R2 files" << std::endl;
       for(int i= 0; i < options.I1s.size(); i++) {
           std::cout << "\t" << options.I1s[i] << " " << filesize(options.I1s[i].c_str()) \
                     << " "  << options.R1s[i] << " " << filesize(options.R1s[i].c_str()) \
                     << " "  << options.R2s[i] << " " << filesize(options.R2s[i].c_str()) <<  std::endl;
       }
  }

  WHITE_LIST_DATA *white_list_data = read_white_list(options.white_list_file);

  process_inputs(options, white_list_data);
}

