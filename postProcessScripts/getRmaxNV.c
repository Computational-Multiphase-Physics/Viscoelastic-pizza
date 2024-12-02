/* Title: Calculating the velocity of the maximum radius of the pizza blob
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"
#include "tag.h"

char filename[80], nameTrack[80];
scalar f[];
vector u;

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  sprintf(nameTrack, "%s", arguments[2]);

  restore (file = filename);
  boundary((scalar *){f, u.x, u.y});

  // tag all liquid parts starts
  scalar d[];
  double threshold = 1e-4;
  foreach(){
    d[] = (f[] > threshold);
  }
  int n = tag (d), size[n];
  for (int i = 0; i < n; i++){
    size[i] = 0;
  }
  foreach_leaf(){
    if (d[] > 0){
      size[((int) d[]) - 1]++;
    }
  }
  int MaxSize = 0;
  int MainPhase = 0;
  for (int i = 0; i < n; i++){
    // fprintf(ferr, "%d %d\n",i, size[i]);
    if (size[i] > MaxSize){
      MaxSize = size[i];
      MainPhase = i+1;
    }
  }
  // tag all liquid parts ends

  face vector s[];
  s.x.i = -1;
  double yMax = -HUGE;
  double vTip, xTP;
  foreach(reduction(max:yMax)){
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d[] == MainPhase) {
      coord n1 = facet_normal (point, f, s);
      double alpha1 = plane_alpha (f[], n1);
      coord segment1[2];
      if (facets (n1, alpha1, segment1) == 2){
        double x1 = x + (segment1[0].x+segment1[1].x)*Delta/2.;
        double y1 = y + (segment1[0].y+segment1[1].y)*Delta/2.;
        if (y1 > yMax){
          yMax = y1;
          xTP = x1;
          vTip = u.y[];
        }
      }
    }
  }
  FILE * fp = ferr;
  fprintf(ferr, "%f %7.6e %7.6e %7.6e\n", t, xTP, yMax, vTip);
  fflush (fp);
  fclose (fp);

  FILE *fp2;
  fp2 = fopen (nameTrack, "a");
  fprintf(fp2, "%f %7.6e %7.6e %7.6e\n", t, xTP, yMax, vTip);
  fclose(fp2);

}
