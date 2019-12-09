#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <cstdlib>
#include <iterator>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;

typedef AT::Site_2                    Site_2;
typedef AT::Point_2                   Point_2;
typedef K::Segment_2                  Segment_2;
typedef K::Intersect_2                Intersect_2;
typedef VD::Locate_result             Locate_result;
typedef VD::Vertex_handle             Vertex_handle;
typedef VD::Face_handle               Face_handle;
typedef VD::Halfedge_handle           Halfedge_handle;
typedef VD::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
typedef std::vector<Face_handle>      facets;

using namespace std;


void print_endpoint(Halfedge_handle e)
{
    if((e->has_source()) && (e->has_target()))
    {
        std::cout <<"(("<< e->source()->point() << ") , " <<"("<< e->target()->point() <<"))"<< std::endl;
    }
    else if((e->has_source()) && !(e->has_target()))
    {
        std::cout <<"("<< e->source()->point() <<")"<<" and target at infinity"<<std::endl;
    }
    else if((e->has_target()) && !(e->has_source()))
    {
        std::cout <<"("<< e->target()->point() <<")"<<" and source at infinity"<< std::endl;
    }
    else if(!(e->has_target()) && !(e->has_source()))
    {
        std::cout << "points at infinity" << std::endl;
    }
}

void print_faces(Face_handle *f)
{
    std::cout << "The halfedges of the Voronoi face are"<< std::endl;
    Halfedge_handle h;
    Ccb_halfedge_circulator ec_start = (*f)->ccb();
    Ccb_halfedge_circulator ec = ec_start;
    do{
        print_endpoint(ec);
    }while ( ++ec != ec_start );
    std::cout << std::endl;
    std::cout << std::endl;
}

int main(int argc, char** argv)
{
  CGAL::Timer clock1;
  CGAL::Timer clock2;
  facets cells;
  char* input_voronoi;
  char* input_moves;
  for(int i=0;i<argc;i++)
  {
	string a = string(argv[i]);
	if(a == "-v")
	{
		input_voronoi = argv[i+1];
	}
	if(a == "-m")
	{
		input_moves = argv[i+1];
	}
  }
  clock1.start();
  std::ifstream ifs(input_voronoi);
  assert( ifs );
  VD vd;
  Site_2 t;
  while ( ifs >> t )
  {
      vd.insert(t);
  }
  ifs.close();
  assert( vd.is_valid() );
  clock1.stop();
  clock2.start();
  std::ifstream ifq(input_moves);
  assert( ifq );
  Point_2 p,q;
  ifq >> p;

  while ( ifq >> q )
  {
    Segment_2 seg(p,q);
    Locate_result lr = vd.locate(p);
    Face_handle* f;
    if ( Vertex_handle* v = boost::get<Vertex_handle>(&lr) )
    {
        *f = (*v)->halfedge()->face();
	cout<<"vertex"<<endl;
    }
    else if ( Halfedge_handle* e = boost::get<Halfedge_handle>(&lr) )
    {
        *f = (*e)->face();
	cout<<"halfedge"<<endl;
    }
    else if ( Face_handle* fp = boost::get<Face_handle>(&lr) )
    {
        f = fp;
	cout<<"face"<<endl;
    }
      Halfedge_handle h;
      while(1)
      {
          Ccb_halfedge_circulator ec_start = (*f)->ccb();
      	  Ccb_halfedge_circulator ec = ec_start;
      do{
        if(ec->has_source() && ec->has_target())
        {
            Segment_2 seg2(ec->source()->point(),ec->target()->point());
            CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type
            result = intersection(seg, seg2);
            if(result)
            {
                h = ec->opposite();
                cells.push_back(*f);
                *f = h->face();
                break;
            }
        }
       }while ( ++ec != ec_start );
        break;
      }
    std::cout << std::endl;
    p = q;
  }
  cout<<"num of cells that the car crosses: "<<cells.size()<<endl;
  cout<<endl;
  cout<<"print halfedges of cells that the car crosses (source_point,target_point)"<<endl;
  cout<<endl;
  for(int i=0;i<cells.size();i++)
  {
      print_faces(&cells[i]);
  }
  ifq.close();
  clock2.stop();
  
  cout<<"Computation time for making the Voronoi diagram is: "<<clock1.time()<<endl;
  cout<<"Computation time for locating the cells is: "<<clock2.time()<<endl;
  
  return 0;
}
