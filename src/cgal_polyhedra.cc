
// Adapted from CGAL example (Author: Pierre Alliez) by Vladimir Fuka.

#include <iostream>
#include <fstream>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>


typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Object_and_primitive_id Object_and_primitive_id;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef CGAL::Bbox_3 Bbox_3;

typedef struct {double x,y,z;} d3;

typedef struct {Polyhedron *poly; Tree * tree;} Polytree;

extern "C" {
 
 void polyhedron_from_file (Polytree  **pptree, const char *fname){
   Polyhedron *polyhedron = new Polyhedron;
   
   std::ifstream in(fname);
   in >> *polyhedron;
   
   Tree *tree = new Tree(polyhedron->facets_begin(),polyhedron->facets_end());
   
   std::cout << "facets: " << polyhedron->size_of_facets() << std::endl;
   std::cout << "halfedges: " << polyhedron->size_of_halfedges() << std::endl;
   std::cout << "vertices: " << polyhedron->size_of_vertices() << std::endl;
   
   tree->accelerate_distance_queries();
   
   *pptree = new Polytree;
   (*pptree)->poly = polyhedron;
   (*pptree)->tree = tree;
 }
 
 void polyhedron_closest (Polytree *ptree, const d3 *query, d3 *near){
   Point query_point(query->x,query->y,query->z);
   
   Point closest = ptree->tree->closest_point(query_point);
   
   near->x = closest.x();
   near->y = closest.y();
   near->z = closest.z();
 }
 
 bool polyhedron_inside(Polytree *ptree, const d3 *query, const d3 *outside_ref){
   Point query_point(query->x,query->y,query->z);
   Point ref_point(outside_ref->x,outside_ref->y,outside_ref->z);
   Segment seg = Segment(query_point,ref_point);
   
   int n = ptree->tree->number_of_intersected_primitives(seg);
   
   return n%2 == 1;
 }

 void polyhedron_bbox(Polytree *ptree, d3 *const min, d3 *const max){
   Bbox_3 bbox = ptree->tree->bbox();
   *min = {bbox.xmin(), bbox.ymin(), bbox.zmin()};
   *max   = {bbox.xmax(), bbox.ymax(), bbox.zmax()};
 }
 
 void polyhedron_finalize(Polytree **pptree){
  delete (*pptree)->tree; (*pptree)->tree = NULL;
  delete (*pptree)->poly; (*pptree)->poly = NULL;
  delete *pptree; *pptree = NULL;
 }
 
}
