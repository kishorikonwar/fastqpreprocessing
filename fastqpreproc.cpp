#include "SamFile.h"
#include "SamValidation.h"

int main(int argc, char ** argv)
{
   // Check for the appropriate number of arguments.
   if(argc != 3)
   {
      printf("./bam <inputFile> <outputFile.sam/bam>\n");
      exit(-1);
   }

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
