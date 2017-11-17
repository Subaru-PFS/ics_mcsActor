
#ifndef fibreid_H_   /* Include guard */
#define fibreid_H_

void add_node(struct id_list **head,int val,int ind);

void add_node_cur(struct id_list **head,int ind);

void add_node_list(struct id_list **head,int val);

void print_list_match(struct match_list *head);

int print_list(struct id_list *head);

int remove_node(struct id_list **head,int val);

int remove_node_list(struct match_list **head,int val);

int remove_node_cur(struct id_list **head,int ind);

int *matchrot(struct fibreid *homepos,struct fibreid *curpos,double theta,double xc,double yc,double xt,double yt,int nhome,int ncur);

void xy_from_thetaphi(double l1,double l2,double theta,double phi,double xc,double yc,double *xf,double *yf);

long *fibremap(struct fibreid *homepos,int nx,int ny,double maxarm, int nhome);

int is_hidden(int *arc_image,long nx,long ny,int xc,int yc,int inner_r, int outer_r,int thresh,float maxd);

int *wherefrom(int *arc_image,long n_x,long n_y,int xc,int yc,int *spiral,int np);

int *create_spiral(int minr,int maxr,int *nval);

int write_fits_simple_long(char filename[],long nx,long ny,long *image);


struct fibreid *read_homes(char *filename,int npos);

int print_currentpos(struct fibreid *currentpos, int ncur);

int adjacent_assigned(struct id_list *un_fibre,long index,int *adjacent_fibre_list);

int fibre_identify_fine(struct fibreid *homepos,struct fibreid *currentpos,long *patrolregion, int nx, int ny, int ncur, int nhome);

int fibre_identify_coarse(struct fibreid *homepos,struct fibreid *currentpos,long *patrolregion,int *arc_image,int *spiral,int *adjacent_fibre_list,int npval,int nx,int ny,int ncur,int nhome);


long *read_fits_long(char *filename, long *nx, long *ny);

int *create_adjacent_fibre_list(struct fibreid *homepos,int nhome, double maxarm);

void  nearest_point(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int nx,int ny);

void  initial_position(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,long *patrolregion,int nx,int ny);

void  newly_singled(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int nx,int ny);

void update_list(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,int idnum,int index,int ind,int indm,fibreid *homepos,fibreid *currentpos);

void  zero_matches(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int *arc_image,int nx,int ny);

void  newly_singled_coarse(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int *arc_image,int nx,int ny,int *adjacent_fibre_list);

void  confused_fibres(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int *arc_image,int nx,int ny,int *spiral,int npval);

double *get_homes(int *image,int n_x, int n_y, int *np,int hmin, double fwhm, int boxsize);

fibreid *centroid_coarse(int *image, int *arc_image,double *homes,int n_x,int n_y,int hmin, double fwhm, int boxsize,int nhome);

fibreid *centroid_fine(int *image, double *homes,double *xp,double *yp,int n_x,int n_y,int hmin, double fwhm, int boxsize,int nhome);

#endif // fibreid_H_
