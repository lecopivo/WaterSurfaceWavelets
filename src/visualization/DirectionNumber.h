#define DIR_NUM 32

#if (DIR_NUM%4)!=0 || DIR_NUM<=0

#error "Number of direction has to be a non-zero multiple of four!"

#endif
