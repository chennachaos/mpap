

  typedef double Scalar;

Scalar PointPolylineDistanceSquared2D(myPoint& p, vector<myLine>& vertices)
{
  //using Scalar = double;
  
  Scalar dSq, INFINITY;
  Scalar xMinusA, yMinusB;
  Scalar xNextMinusA, yNextMinusB;
  
  Scalar xMinusASq, yMinusBSq;
  Scalar xNextMinusASq, yNextMinusBSq;

  xMinusA = vertices[0].x - p.x;
  yMinusB = vertices[0].y - p.y;

  xMinusASq = xMinusA * xMinusA;
  yMinusBSq = yMinusB * yMinusB;

  xNextMinusA = vertices[1].x - p.x;
  yNextMinusB = vertices[1].y - p.y;

  xNextMinusASq = xNextMinusA * xNextMinusA;
  yNextMinusBSq = yNextMinusB * yNextMinusB;

  // Compute distance to first segment
  Line l = { vertices[i], vertices[i+1] - vertices[i] };

  Scalar t;

  dSq = PointLineDistanceSquared2D(p, l, FALSE, t);

  // If closest point not on segment, check appropriate end point

  if(t < 0)
  {
    dSq = min(dSq, xMinusASq + yMinusBSq);
  }
  else if (t > 1)
  {
   dSq = min(dSq, xNextMinusASq + yNextMinusBSq);
  }

  // Go through each successive segment, rejecting if possible,
  // and computing the distance squared if not rejected.

  for(i=1; i<(nSegments-1); i++)
  {
    // Rejection test
    if (((abs(xMinusASq) > dSq) && (abs(xNextMinusASq) <= dSq)
    && (xMinusA * xNextMinusA > 0)) ||
    ((abs(yMinusBSq) > dSq) && (abs(yNextMinusBSq) <= dSq)
    && (yMinusB * yNextMinusB > 0)) )
    {
      if (i != nSegments - 2)
      {
        xMinusA = xNextMinusA;
        yMinusB = yNextMinusB;

        xNextMinusA = vertices[i + 2].x - p.x;
        yNextMinusB = vertices[i + 2].y - p.y;
      }
      continue;
    }

    // Rejection test failed - check distance to line

    Line l = { vertices[i], vertices[i+1] - vertices[i] };

    Scalar t;

    dSq = PointLineDistanceSquared3D(p, l, FALSE, t);

    // If closest point not on segment, check appropriate end point

    if(t < 0)
    {
        dSq = MIN(dsq, xMinusASq + yMinusBSq );
    }
    else if (t > 1)
    {
        dSq = MIN(dsq, xNextMinusASq + yNextMinusBSq );
    }

    if(i != nSegments - 2)
    {
        xMinusA = xNextMinusA;
        yMinusB = yNextMinusB;

        xNextMinusA = vertices[i + 2].x - p.x;
        yNextMinusB = vertices[i + 2].y - p.y;
    }
  }

  return dSq;
}





Scalar PointPolylineDistanceSquared3D(Point p, Point vertices[], int nSegments)
{

Scalar dSq = INFINITY;
Scalar xMinusA, yMinusB, zMinusC;
Scalar xNextMinusA, yNextMinusB, zNextMinusC

Scalar xMinusASq, yMinusBSq, zMinusCSq;
Scalar xNextMinusASq, yNextMinusBSq, zNextMinusCSq;
xMinusA = vertices[0].x - p.x;
yMinusB = vertices[0].y - p.y;
zMinusC = vertices[0].z - p.z;
xMinusASq = xMinusA * xMinusA;
yMinusBSq = yMinusB * yMinusB;
zMinusCSq = zMinusC * zMinusC;
xNextMinusA = vertices[1].x - p.x;
yNextMinusB = vertices[1].y - p.y;
zNextMinusC = vertices[1].z - p.z;
xNextMinusASq = xNextMinusA * xMNextinusA;
yNextMinusBSq = yNextMinusB * yNextMinusB;
zNextMinusCSq = zNextMinusC * zNextMinusC;
// Compute distance to first segment
Line l = { vertices[i], vertices[i+1] - vertices[i] };
Scalar t;
dSq = PointLineDistanceSquared3D(p, l, FALSE, t)
// If closest point not on segment, check appropriate end point
if (t < 0) {
dSq = MIN(dsq, xMinusASq + yMinusBSq + zMinusCSq);
} else if (t > 1) {
dSq = MIN(dsq, xNextMinusASq + yNextMinusBSq + zNextMinusCSq);
}
// Go through each successive segment, rejecting if possible,
// and computing the distance squared if not rejected.
for (i = 1; i < nSegments - 1; i++) {
// Rejection test
if (((Abs(xMinusASq) > dSq) && (Abs(xNextMinusASq) <= dSq)
&& (xMinusA * xNextMinusA > 0)) ||
((Abs(yMinusBSq) > dSq) && (Abs(yNextMinusBSq) <= dSq)
&& (yMinusB * yNextMinusB > 0)) ||
((Abs(zMinusCSq) > dSq) && (Abs(zNextMinusCSq) <= dSq)
&& (zMinusC * zNextMinusC > 0))) {
if (i != nSegments - 2) {


xMinusA = xNextMinusA;
yMinusB = yNextMinusB;
zMinusC = zNextMinusC;
xNextMinusA = vertices[i + 2].x - p.x;
yNextMinusB = vertices[i + 2].y - p.y;
zNextMinusC = vertices[i + 2].z - p.z;
}
continue;
}
// Rejection test failed - check distance to line
Line l = { vertices[i], vertices[i+1] - vertices[i] };
Scalar t;
dSq = PointLineDistanceSquared3D(p, l, FALSE, t)
// If closest point not on segment, check appropriate end point
if (t < 0) {
dSq = MIN(dsq, xMinusASq + yMinusBSq + zMinusCSq);
} else if (t > 1) {
dSq = MIN(dsq, xNextMinusASq + yNextMinusBSq + zNextMinusCSq);
}
if (i != nSegments - 2) {
xMinusA = xNextMinusA;
yMinusB = yNextMinusB;
zMinusC = zNextMinusC;
xNextMinusA = vertices[i + 2].x - p.x;
yNextMinusB = vertices[i + 2].y - p.y;
zNextMinusC = vertices[i + 2].z - p.z;
}
}
return dSq;
}

