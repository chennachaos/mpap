
#ifndef incl_MathGeom_h
#define incl_MathGeom_h


#include <cmath>

#include "MathBasic.h"



double inline triangleArea2D(double *x1, double *x2, double *x3)
{
  return 0.5 * (x1[0]*(x2[1]-x3[1])+x2[0]*(x3[1]-x1[1])+x3[0]*(x1[1]-x2[1]));
}




double inline triangleArea3D(double *x1, double *x2, double *x3, double *n = NULL)
{
  double cross[3] = { -(x2[2]-x1[2])*(x3[1]-x1[1])+(x2[1]-x1[1])*(x3[2]-x1[2]),
                      +(x2[2]-x1[2])*(x3[0]-x1[0])-(x2[0]-x1[0])*(x3[2]-x1[2]),
                      -(x2[1]-x1[1])*(x3[0]-x1[0])+(x2[0]-x1[0])*(x3[1]-x1[1]) },

         crossNorm = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);

  if (n != NULL) for (int i=0; i<3; i++) n[i] = cross[i] / crossNorm; 

  return 0.5 * crossNorm;
}




double inline triangleAspectRatio(double *x1, double *x2, double *x3)
{
  double s2[3] = {(x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]),
                  (x2[0]-x3[0])*(x2[0]-x3[0]) + (x2[1]-x3[1])*(x2[1]-x3[1]),
                  (x3[0]-x1[0])*(x3[0]-x1[0]) + (x3[1]-x1[1])*(x3[1]-x1[1])},
         s[3]  = { sqrt(s2[0]), sqrt(s2[1]), sqrt(s2[2]) },
         h1, h2, h3, h4, ri2, ro2;

  h1 = s2[0]*s[0] + s2[1]*s[1] + s2[2]*s[2];
  h2 = 2. * s[0] * s[1] * s[2];
  h3 =      s2[0]*s[1] + s2[1]*s[0] + s2[1]*s[2] + s2[2]*s[1] + s2[2]*s[0] + s2[0]*s[2];
  h4 = 4. * (s[0] + s[1] + s[2]);

  ri2 = (h3 - h1 - h2) / h4;

  h1 = s2[0] * s2[1] * s2[2];
  h2 = s2[0] * s2[1] + s2[1] * s2[2] + s2[2] * s2[0];
  h3 = s2[0] * s2[0] + s2[1] * s2[1] + s2[2] * s2[2];

  ro2 = h1 / (h2 + h2 - h3);

  return sqrt(ri2 / ro2);
}




template<typename Type> bool inline pointInTriangle2D(Type *x1, Type *x2, Type *x3, Type *xp)
{
  if ((x1[1]-x2[1])*(xp[0]-x1[0]) + (x2[0]-x1[0])*(xp[1]-x1[1]) < 0.) return false;
  if ((x2[1]-x3[1])*(xp[0]-x2[0]) + (x3[0]-x2[0])*(xp[1]-x2[1]) < 0.) return false;
  if ((x3[1]-x1[1])*(xp[0]-x3[0]) + (x1[0]-x3[0])*(xp[1]-x3[1]) < 0.) return false;

  return true;
}




bool inline pointInQuadrilateral2D(double *x1, double *x2, double *x3, double *x4, double *xp)
{
  if (pointInTriangle2D(x1,x2,x3,xp) || pointInTriangle2D(x1,x3,x4,xp)) return true;

  return false;
}




double inline dist2(double *x1, double *x2, int ndm)
{
  double d = 0;

  for (int i=0; i<ndm; i++) d += (x1[i] - x2[i]) * (x1[i] - x2[i]);

  return d;
}




void inline pointAndLine2D(double *x0, double *ray, double *xp, 
                           double *alpha, double *d = NULL, double *l = NULL)
{
  double rn = sqrt(ray[0]*ray[0]+ray[1]*ray[1]), irn = 1. / rn;

  *alpha = ((xp[0]-x0[0])*ray[0] + (xp[1]-x0[1])*ray[1])*irn*irn;
  
  if (d != NULL)
  {
    *d = (ray[0]*(xp[1]-x0[1]) - ray[1]*(xp[0]-x0[0])) * irn;

    if (l != NULL) *l = rn;
  }
  return;
}




void inline pointAndLine3D(double *x0, double *ray, double *xp, 
                           double *alpha, double *d = NULL, double *l = NULL)
{
  double rn = sqrt(ray[0]*ray[0]+ray[1]*ray[1]+ray[2]*ray[2]), irn = 1. / rn;

  *alpha = ((xp[0]-x0[0])*ray[0] + (xp[1]-x0[1])*ray[1] + (xp[2]-x0[2])*ray[2])*irn*irn;
  
  if (d == NULL) return;

  double n[3] = { xp[0]-x0[0]-*alpha*ray[0],
                  xp[1]-x0[1]-*alpha*ray[1],
                  xp[2]-x0[2]-*alpha*ray[2] };

  *d = sqrt(dot3(n,n));

  if (l != NULL) *l = rn;

  return;
}




void inline analysePointAndEdge2D(double *x1, double *x2, double *xp, 
                                  double *a, double *d, double *l)
{
  double ray[2] = { x2[0]-x1[0], x2[1]-x1[1] };

  pointAndLine2D(x1,ray,xp,a,d,l);

  return;
}




bool inline pointOnEdge2D(double *x1, double *x2, double *xp, 
                          double *a = NULL, double *d = NULL, double *l = NULL)
{
  double ray[2] = { x2[0]-x1[0], x2[1]-x1[1] }, alpha, dist, len, tol = 1.e-8, tol4 = .4;

  pointAndLine2D(x1,ray,xp,&alpha,&dist,&len);

  if (alpha < -tol || alpha > 1.+tol) return false;

  if ((alpha+tol)*(1.-alpha+tol)*tol4 < abs(dist/len)) return false;

  if (a != NULL) { *a = alpha; if (d != NULL) { *d = dist; if (l != NULL) *l = len; } }

  return true;
}




bool inline pointOnEdge3D(double *x1, double *x2, double *xp, 
                          double *a = NULL, double *d = NULL, double *l = NULL)
{
  double ray[3] = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] }, 
         alpha, dist, len, tol4 = .4;

  pointAndLine3D(x1,ray,xp,&alpha,&dist,&len);

  if (alpha < .0 || alpha > 1.) return false;

  if (alpha*(1.-alpha)*tol4 > abs(dist/len)) return false;

  if (a != NULL) { *a = alpha; if (d != NULL) { *d = dist; if (l != NULL) *l = len; } }

  return true;
}




bool inline pointInTriangle2D(double *x1, double *x2, double *x3, double *xp,
                              double *p, double *a = NULL)
{
  double ac[3];

  ac[0] = triangleArea2D(x1,x2,xp); if (ac[0] < .0) return false;
  ac[1] = triangleArea2D(x2,x3,xp); if (ac[1] < .0) return false;
  ac[2] = triangleArea2D(x3,x1,xp); if (ac[2] < .0) return false;

  double area = ac[0] + ac[1] + ac[2];

  for (int i=0; i<3; i++) p[i] = ac[i] / area;

  if (a != NULL) *a = area;

  return true;
}




bool inline pointInTriangle3D(double *x1, double *x2, double *x3, double *xp, 
                               double *p = NULL,  double *d = NULL, double *a = NULL)
{
  double n[3], h[3], d1[3], d2[3], m[3], ac[3], area, dist;

  area = triangleArea3D(x1,x2,x3,n);

  h[0] = xp[0] - x1[0];
  h[1] = xp[1] - x1[1];
  h[2] = xp[2] - x1[2];

  dist = dot3(n,h); // distance of xp to plane of triangle

  h[0] = xp[0] - dist * n[0]; // projection of xp into plane of triangle
  h[1] = xp[1] - dist * n[1];
  h[2] = xp[2] - dist * n[2];

  //std::cout << "\n" 
  //          << x1[0] << "," << x1[1] << "," << x1[2] << " x1\n"
  //          << x2[0] << "," << x2[1] << "," << x2[2] << " x2\n"
  //          << x3[0] << "," << x3[1] << "," << x3[2] << " x3\n"
  //          << xp[0] << "," << xp[1] << "," << xp[2] << " xp\n"
  //          <<  n[0] << "," <<  n[1] << "," <<  n[2] << " n\n"
  //          <<  h[0] << "," <<  h[1] << "," <<  h[2] << " h\n";

  d1[0] = x2[0] - h[0];
  d1[1] = x2[1] - h[1];
  d1[2] = x2[2] - h[2];
                  
  d2[0] = x3[0] - h[0];
  d2[1] = x3[1] - h[1];
  d2[2] = x3[2] - h[2];

  cross3(m,d1,d2);

  //std::cout << dot3(m,n) << " 23\n";

  if (dot3(m,n) < .0) return false;

  ac[0] = .5 * sqrt(dot3(m,m)) / area;

  d1[0] = x1[0] - h[0];
  d1[1] = x1[1] - h[1];
  d1[2] = x1[2] - h[2];

  cross3(m,d2,d1);

  //std::cout << dot3(m,n) << " 31\n";

  if (dot3(m,n) < .0) return false;

  ac[1] = .5 * sqrt(dot3(m,m)) / area;

  d2[0] = x2[0] - h[0];
  d2[1] = x2[1] - h[1];
  d2[2] = x2[2] - h[2];

  cross3(m,d1,d2);

  //std::cout << dot3(m,n) << " 12\n";

  if (dot3(m,n) < .0) return false; 

  ac[2] = .5 * sqrt(dot3(m,m)) / area;

  double l[3] = { sqrt(dist2(x2,x3,3)), sqrt(dist2(x3,x1,3)), sqrt(dist2(x1,x2,3)) };

  if ((ac[1]*ac[2]*l[0]+ac[2]*ac[0]*l[1]+ac[0]*ac[1]*l[2])*.4 < abs(dist)) return false;

  if (p == NULL) return true; else for (int i=0; i<3; i++) p[i] = ac[i];

  if (d == NULL) return true; else *d = dist;

  if (a == NULL) return true; else *a = area;

  return true;
}







bool inline edgesIntersect2D(double *x1, double *x2, double *X1, double *X2)
{
  if (x1 == X1 || x1 == X2) return false;
  if (x2 == X1 || x2 == X2) return false;

  double dx  = x2[0] - x1[0],
         dy  = x2[1] - x1[1],
         DX  = X2[0] - X1[0],
         DY  = X2[1] - X1[1],
         xX  = X1[0] - x1[0],
         yY  = X1[1] - x1[1],
         tol = 1.e-8,
         fct = dy * DX - dx * DY,
         alpha, beta;

  if (fabs(fct) < 1.e-10) return false;

  fct = 1./fct;

  alpha = (yY * DX - xX * DY) * fct;

  if (alpha > 1.+tol || alpha < -tol) return false;

  beta  = (dx * yY - dy * xX) * fct;

  if (beta > 1.+tol || beta < -tol) return false;

  if (beta < tol && (alpha > 1.-tol || alpha < tol)) return false;

  if (beta > 1.-tol && (alpha < tol || alpha > 1.-tol)) return false;

  return true;
}







double inline cosAngleBetweenVectors(double *v1, double *v2, int ndm)
{
  double a = v1[0] * v2[0],
         b = v1[0] * v1[0],
         c = v2[0] * v2[0];

  for (int i=1; i<ndm; i++) 
  { 
    a += v1[i] * v2[i];
    b += v1[i] * v1[i];
    c += v2[i] * v2[i];
  }

  return a / sqrt(b*c);
}







double inline boundaryAngleBetweenTwoEdges2D(double *v1, double *v2, int ndm)
{
  double angle = 3.1415927 - myAcos(cosAngleBetweenVectors(v1,v2,2));

  if (-v1[0]*v2[1] + v1[1]*v2[0] >= 0.) angle = 6.2831853 - angle;

  return angle;
}






void inline circleThroughTwoPoints
              (double* dat, double rad, double x1, double y1, double x2, double y2)
{
  dat[0] = .5*(x1+x2);
  dat[1] = .5*(y1+y2);

  double x0[2] = { .5*(x1+x2), .5*(y1+y2) },
          n[2] = { y1 - y2, x2 - x1 },
          k = rad / sqrt(n[0]*n[0]+n[1]*n[1]),
         a, r2 = rad * rad, r;

  n[0] *= k;
  n[1] *= k;

  a = .7;

  r = (dat[0]+a*n[0]-x1)*(dat[0]+a*n[0]-x1)
     +(dat[1]+a*n[1]-y1)*(dat[1]+a*n[1]-y1) - r2;

  //cout << r << "\n";

  while (abs(r) > 1.e-12)
  {
    k = (dat[0]+a*n[0]-x1)*(n[0]+n[0])
       +(dat[1]+a*n[1]-y1)*(n[1]+n[1]);

    a -= r / k;

    r = (dat[0]+a*n[0]-x1)*(dat[0]+a*n[0]-x1)
       +(dat[1]+a*n[1]-y1)*(dat[1]+a*n[1]-y1) - r2;

    //cout << r << "\n";
  }
  //exit(0);

  dat[0] += a * n[0];
  dat[1] += a * n[1];

  dat[2] = atan2(y1-dat[1],x1-dat[0]);
  dat[3] = atan2(y2-dat[1],x2-dat[0]);

  if (dat[3] < dat[2]) dat[3] += 3.1415926535897932384626433 + 3.1415926535897932384626433;

  //cout << dat[0] << "," << dat[1] << "\n";
  //cout << dat[2] << "," << dat[3] << "\n";

  return;
}




#endif

