#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <cmath>

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
typedef VD::Face_iterator             face_iterator;
typedef std::vector<Face_handle>      facets;

using namespace std;


int new_rand(int a, int b, int c)
{
    while(a == b)
    {
        a = rand() % c + 0;
    }
    return a;
}

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
  srand(time(NULL));
  int count_points = 0;
  facets cells;
  char* input_voronoi;
  char* input_moves;
  int num_faces=0;
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
  CGAL::Timer clock1;
  CGAL::Timer clock2;
  clock1.start();
  std::ifstream ifs(input_voronoi);
  assert( ifs );
  VD vd;
  Site_2 t;
  while ( ifs >> t )
  {
      num_faces++;
      vd.insert(t);
  }
  ifs.close();
  assert( vd.is_valid() );
  clock1.stop();
  int num_forbidden;
  num_forbidden = floor((1.0/3.0)*num_faces);
  int forbiddenCells[num_forbidden];
  for(int i=0;i<num_forbidden;i++)
  {
      forbiddenCells[i] = -1;
  }
  int random[num_forbidden];
  for(int i=0;i<num_forbidden;i++)
  {
      random[i] = -1;
  }
  for(int i=0;i<num_forbidden;i++)
  {
      random[i] = rand() % num_faces + 0;
      for(int j=0;j<i;j++)
      {
          if(random[i] == forbiddenCells[j])
          {
              random[i] = new_rand(random[i],forbiddenCells[j],num_faces);
          }
      }
      forbiddenCells[i] = random[i];
  }

  clock2.start();
  std::ifstream ifq(input_moves);
  assert( ifq );
  Point_2 p,q;
  ifq >> p;
  Halfedge_handle h;
  while ( ifq >> q )
  {
    count_points++;
    bool forbidden;
    bool check_result = false;
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
      while(1)
      {
	int counter=0;
	forbidden = false;
        for(face_iterator fi = vd.faces_begin();fi != vd.faces_end();fi++)
        {
		if(Face_handle(fi) == *f)
		{
                	for(int k=0;k<num_forbidden;k++)
                	{
                    		if(forbiddenCells[k] == counter)
				{
					forbidden = true;
					break;
                    		}
                	}
           	 }
		if(forbidden == true)
		{
			break;
		}
		counter++;
         }
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
				check_result = true;
				if(forbidden == false)
				{
					h = ec->opposite();
					cells.push_back(*f);
					print_faces(f);
					*f = h->face();
					break;
				}
				else
				{
					if(count_points == 1)
					{
						cout<<"Exit! First point can't be in a forbidden shell!"<<endl;
						return -1;
					}
					cout<<"Forbidden cell! The halfedges around the face that the car will cross are:"<<endl;
					Halfedge_handle h1 = h;
					do{
						h1 = h1->next();
						print_endpoint(h1);
					}while(h1 != Halfedge_handle(ec));
					print_endpoint(ec);
					h = ec->opposite();
					*f = h->face();
					break;
				}
			}
		}
	}while ( ++ec != ec_start );
        break;
      }
    std::cout << std::endl;
    if(!((forbidden == true) && (check_result == false)))
    {
	p = q;
    }
  }
  ifq.close();
  clock2.stop();

  cout<<"Computation time for making the Voronoi diagram is: "<<clock1.time()<<endl;
  cout<<"Computation time for locating the cells is: "<<clock2.time()<<endl;

  return 0;
}
