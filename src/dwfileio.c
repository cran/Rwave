
/*******************************************************************
*              (c) Copyright  1997                                 *
*                         by                                       *
* Author: Rene Carmona, Bruno Torresani, Wen L. Hwang and A. Wang  *
*                 Princeton University                             *
*                 All right reserved                               *
*******************************************************************/



#include "Swave.h"


/****************************************************************
*  Function: input_signal:
*  -----------------------
*  Read a string of floating number.
*
*    fname: filename
*    Pic: place to store
*    size: size to read 
*
****************************************************************/
void input_signal(fname,Pic,size)
     char *fname;
     int size;
     float **Pic;
{
  FILE *fp;
  int k;
  float tmp;
  
  if ((fp = fopen(fname,"r")) == NULL) {
    printf("Can't open file %s\n", fname);
  }

  if(!(*Pic = (float *)malloc(sizeof(float) * size)))
    error("Memory allocation failed for *Pic in input.c \n");

  for(k = 0; k < size; k++){
    if (!fscanf(fp,"%f\n",&tmp)) error("error in reading \n");
    (*Pic)[k] = tmp;
  }
  
  fclose(fp);
}


/****************************************************************
*  Function: init_filename:
*  ------------------------
*  Initialize a string for filename.
*
*    filename: filename
*
****************************************************************/

init_filename(filename)
     char filename[];
{
  filename[0] = '\0';
}

/****************************************************************
*  Function: filename_given:
*  -------------------------
*  Given a filename
*
****************************************************************/


filename_given(filename,name)
     char filename[];
     char *name;
{
  int lg;
  int i;

  lg = strlen(name);
  for(i = 0 ; i <= lg; i++)
    filename[i] = name[i];
}

/****************************************************************
*  Function: filename_given:
*  -------------------------
*  Concatenate t to s, resultant as buff
*
****************************************************************/

strconcate(s , t, buff)
     char s[], t[], buff[];
{
  int i, j;
  for (i = 0; s[i] != '\0'; i++)
    buff[i] = s[i];
  for (j= 0; t[j] != '\0'; j++)
    buff[i++] = t[j];
  buff[i] = '\0';
}


/****************************************************************
*  Function: filename_inc:
*  -------------------------
*  "Increment" a filename by appending inc to its end
*
****************************************************************/

filename_inc(filename, inc)
     char filename[];
     int inc;
{
  char t[4];
  char buff[STRING_SIZE];
  int tmp;

  t[0] = '.';
  tmp = inc/10;
  if (tmp == 0) {
    t[1] = inc + '0';
    t[2] = '\0';
  }
  else { 
    t[1] = tmp + '0';
    t[2] = inc % 10 + '0';
    t[3] = '\0';
  }
  strconcate(filename, t, filename);
}

      
/****************************************************************
*  Function: output_signal:
*  ------------------------
*  Write signal to disk
*  
*  s: signal
*  np: size of signal
*  fname: filename
*
****************************************************************/

void output_signal( s, np, fname)
     char *fname;
     float *s;
     int np;
{
  FILE *fp;
  int i;
  
  if ((fp = fopen(fname,"w")) == NULL) {
    printf("Can't open file %s\n", fname);
  }
  
  for ( i = 0 ; i < np; i++ )
    fprintf( fp, "%f\n", s[i] );
  fclose( fp );
}

/****************************************************************
*  Function: output_array:
*  -----------------------
*  Output array of floating numbers
*  
*  array: floating numbers
*  nrow: number of row in array
*  ncol: number of column in array
*  file_name: filename for output
*
****************************************************************/


void output_array( array, nrow, ncol, file_name )
     float **array;
     int nrow, ncol;
     char *file_name;
{
  FILE *fp;
  int x, y;

  if ((fp = fopen(file_name,"w")) == NULL) {
    printf("Can't open file %s\n", file_name);
  }

  for ( x = 0; x < nrow; x++ )  {
    for ( y = 0; y < ncol; y++ )
      fprintf( fp, "%f ", array[x][y] );
  }
  fclose( fp );
}




