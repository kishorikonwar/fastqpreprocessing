#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__
#include "fastqprocess.h"

WHITE_LIST_DATA *read_white_list(const string &white_list_file);

long filesize(const char *filename);
long getFileSize(const std::string &fileName);




#endif
