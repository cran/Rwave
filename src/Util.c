
/****************************************************************
*               (c) Copyright  1997                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  Princeton University                         *
*                  All right reserved                           *
****************************************************************/






int find2power(n)
int n;
{
   long m, m2;

   m = 0;
   m2 = 1<<m; /* 2 to the power of m */
   while (m2-n < 0) {
        m++;
        m2 <<= 1; /* m2 = m2*2 */
   }
   return(m);
}

void error(char error_text[])
{
  printf("%s\n",error_text);
  exit(0);
  return;
}
