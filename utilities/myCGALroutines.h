#ifndef MY_CGAL_ROUTINES_h
#define MY_CGAL_ROUTINES_h


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


typedef CGAL_Polyhedron::HalfedgeDS  HalfedgeDS;


// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS>
{
  public:
    std::vector<CGAL_Point>  &vertices;
    std::vector<int>    &faces;
    int  elType;

    polyhedron_builder( std::vector<CGAL_Point> &_coords, std::vector<int> &_tris, int eT=3 ) : vertices(_coords), faces(_tris), elType(eT) {}
    
    void operator()( HDS& hds)
    {
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
 
        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        
        if(elType == 3)
          B.begin_surface( vertices.size(), faces.size()/3 );
        else
          B.begin_surface( vertices.size(), faces.size()/4 );

       // add the polyhedron vertices
       for( int i=0; i<(int)vertices.size(); i++ )
       {
         B.add_vertex(vertices[i]);
       }
   
       // add the polyhedron triangles

       if(elType == 3)
       {
         for( int i=0; i<(int)faces.size(); i+=3 )
         {
           B.begin_facet();
           B.add_vertex_to_facet( faces[i+0] );
           B.add_vertex_to_facet( faces[i+1] );
           B.add_vertex_to_facet( faces[i+2] );
           B.end_facet();
         }
       }
       else
       {
         for( int i=0; i<(int)faces.size(); i+=4 )
         {
           B.begin_facet();
           B.add_vertex_to_facet( faces[i+0] );
           B.add_vertex_to_facet( faces[i+1] );
           B.add_vertex_to_facet( faces[i+3] );
           B.end_facet();

           B.begin_facet();
           B.add_vertex_to_facet( faces[i+0] );
           B.add_vertex_to_facet( faces[i+3] );
           B.add_vertex_to_facet( faces[i+2] );
           B.end_facet();
         }
       }

       // finish up the surface
       B.end_surface();
    }
};



inline CGAL_Point_inside*  createPointerToCGALPointInsideFromCGALTree(CGAL_Tree& tree)
{
    CGAL_Point_inside*  inside_tester = new CGAL_Point_inside(tree);
    
    return inside_tester;
}


#endif

