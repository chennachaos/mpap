#ifndef MY_CGAL_ROUTINES2_h
#define MY_CGAL_ROUTINES2_h


#include <algorithm>
#include <iostream>
#include <list>
#include <vector>
#include <fstream>
#include <limits>
#include <boost/foreach.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
 

/*
typedef CGAL::Simple_cartesian<double> CGAL_K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef CGAL_K::FT CGAL_FT;
typedef CGAL_K::Ray_3 CGAL_Ray;
typedef CGAL_K::Line_3 CGAL_Line;
typedef CGAL_K::Point_3 CGAL_Point;
typedef CGAL_K::Segment_3 CGAL_Segment;
typedef CGAL_K::Triangle_3 CGAL_Triangle;
//typedef std::list<Triangle>::iterator Iterator;
//typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
//typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
//typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef CGAL::Polyhedron_3<CGAL_K>  CGAL_Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<CGAL_Polyhedron>  CGAL_Primitive;
typedef CGAL::AABB_traits<CGAL_K, CGAL_Primitive>  CGAL_Traits;
typedef CGAL::AABB_tree<CGAL_Traits>  CGAL_Tree;

typedef CGAL::Side_of_triangle_mesh<CGAL_Polyhedron, CGAL_K>  CGAL_Point_inside;
*/

inline int  checkBoundedSideCGAL(CGAL::Bounded_side res)
{
    int c=0;
    if( (res == CGAL::ON_BOUNDED_SIDE) || (res == CGAL::ON_BOUNDARY) )
      c = 1;

    return c;
}




#endif

