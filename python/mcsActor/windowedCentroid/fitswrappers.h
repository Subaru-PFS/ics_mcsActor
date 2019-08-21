#ifndef fitswrappers_H_   /* Include guard */
#define fitswrappers_H_

int write_fits_simple_int(char filename[],long nx,long ny,int *image);
int write_fits_simple_long(char filename[],long nx,long ny,long *image);
int write_fits_simple_longlong(char filename[],long nx,long ny,long *image);
int *read_fits_int(char *filename, int *nx, int *ny);
long *read_fits_long(char *filename, long *nx, long *ny);

#endif // fitswrappers_H_
