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

WHITE_LIST_DATA *read_write_list(const string &white_list_file) {
   // be careful of caps 
   char ATCG[] = {'A', 'T', 'C', 'G'};

   fstream newfile;
   WHITE_LIST_DATA *white_list_data = new  WHITE_LIST_DATA;

   newfile.open(white_list_file, ios::in); //open a file to perform read operation using file object
   int k = 0;
   if (newfile.is_open()){   //checking whether the file is open
      string tp;
      white_list_data->barcodes.push_back(tp);
      while(getline(newfile, tp)){ //read data from file object and put it into string.
       //  cout << tp << "\n"; //print the data of the string
         for(int i=0; i < tp.size(); i++) {
             for(int j=0; j < 4; j++) {
                 char c = tp[i];
                 tp[i] = ATCG[j];
                 white_list_data->mutations.insert( {tp, true} );
                 tp[i] = c;
             }
         }
         if( k%100000 == 0 && k != 0 ) printf("%d\n", k);
         k++;
      }
      newfile.close(); //close the file object.
   }

   std::cout << "size of whitelist " << white_list_data->mutations.size() << std::endl;
   return white_list_data;
} 


