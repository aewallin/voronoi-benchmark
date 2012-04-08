// Boost.Polygon library voronoi_benchmark.cpp file

//          Copyright Andrii Sydorchuk 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#include <iomanip>
#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/timer.hpp>

typedef boost::int32_t int32;

// Include for Boost.Polygon Voronoi library.
#include "boost/polygon/voronoi.hpp"
typedef boost::polygon::voronoi_diagram<double> VD_BOOST;

// Includes for CGAL library.
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
typedef CGAL::Quotient<CGAL::MP_Float> ENT;
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double> CK;
typedef CGAL::Simple_cartesian<ENT> EK;
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<
    CK, CGAL::Field_with_sqrt_tag, EK, CGAL::Field_tag> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt> SDT_CGAL;
typedef SDT_CGAL::Point_2 Point_CGAL;
typedef SDT_CGAL::Site_2 Site_CGAL;

// Include Boost.Polygon library.
#include "boost/polygon/polygon.hpp"
typedef boost::polygon::point_data<int32> POINT_POLYGON;
typedef boost::polygon::directed_line_segment_data<int32> SEGMENT_POLYGON;
typedef boost::polygon::directed_line_segment_set_data<int32> SSD_POLYGON;

const int RANDOM_SEED = 27;
const int NUM_TESTS = 3;
const int NUM_SEGMENTS[] = {10, 100, 1000 }; //, 10000, 100000, 1000000};
const int NUM_RUNS[] = {100, 10, 10, 10, 10, 1};
std::ofstream bf("benchmark_segments.txt", std::ios_base::out | std::ios_base::app);
boost::timer timer;

void format_line(int num_points, int num_tests, double time_per_test) {
  bf << "| " << std::setw(16) << num_points << " ";
  bf << "| " << std::setw(15) << num_tests << " ";
  bf << "| " << std::setw(17) << time_per_test << " ";
  bf << "|" << std::endl;
}

std::vector<double> get_intersection_runtime() {
  std::vector<double> running_times;
  boost::mt19937 gen(RANDOM_SEED);
  for (int i = 0; i < NUM_TESTS; ++i) {
    timer.restart();
    for (int j = 0; j < NUM_RUNS[i]; ++j) {
      SSD_POLYGON ssd;
      for (int k = 0; k < NUM_SEGMENTS[i]; ++k) {
        int32 x1 = gen();
        int32 y1 = gen();
        int32 dx = (gen() & 1023) + 1;
        int32 dy = (gen() & 1023) + 1;
        ssd.insert(SEGMENT_POLYGON(POINT_POLYGON(x1, y1),
                                   POINT_POLYGON(x1 + dx, y1 + dy)));
      }
      ssd.clean();
    }
    double time_per_test = timer.elapsed() / NUM_RUNS[i];
    running_times.push_back(timer.elapsed());
  }
  return running_times;
}

void run_voronoi_test(const std::vector<double> &running_times) {
  boost::mt19937 gen(RANDOM_SEED);
  bf << "Boost.Polygon Voronoi of Segments:\n";
  for (int i = 0; i < NUM_TESTS; ++i) {
    timer.restart();
    for (int j = 0; j < NUM_RUNS[i]; ++j) {
      SSD_POLYGON ssd;
      VD_BOOST vd;
      for (int k = 0; k < NUM_SEGMENTS[i]; ++k) {
        int32 x1 = gen();
        int32 y1 = gen();
        int32 dx = (gen() & 1023) + 1;
        int32 dy = (gen() & 1023) + 1;
        ssd.insert(SEGMENT_POLYGON(POINT_POLYGON(x1, y1), POINT_POLYGON(x1 + dx, y1 + dy)) );
      }
      ssd.clean();
      boost::polygon::construct_voronoi_segments(ssd, &vd);
    }
    double time_per_test = (timer.elapsed() - running_times[i]) / NUM_RUNS[i];
    format_line(NUM_SEGMENTS[i], NUM_RUNS[i], time_per_test);
  }
  bf << "\n";
}

void run_cgal_test(const std::vector<double> &running_times) {
  boost::mt19937 gen(RANDOM_SEED);
  bf << "CGAL Triangulation of Segments:\n";
  for (int i = 0; i < NUM_TESTS; ++i) {
    timer.restart();
    for (int j = 0; j < NUM_RUNS[i]; ++j) {
      SSD_POLYGON ssd;
      for (int k = 0; k < NUM_SEGMENTS[i]; ++k) {
        int32 x1 = gen();
        int32 y1 = gen();
        int32 dx = (gen() & 1023) + 1;
        int32 dy = (gen() & 1023) + 1;
        ssd.insert(SEGMENT_POLYGON(POINT_POLYGON(x1, y1), POINT_POLYGON(x1 + dx, y1 + dy)));
      }
      ssd.clean();
      SDT_CGAL dt;
      for (SSD_POLYGON::iterator_type it = ssd.begin(); it != ssd.end(); ++it) {
        dt.insert(Site_CGAL::construct_site_2( Point_CGAL(it->low().x(), it->low().y()), Point_CGAL(it->high().x(), it->high().y())));
      }
    }
    double time_per_test = (timer.elapsed() - running_times[i]) / NUM_RUNS[i];
    format_line(NUM_SEGMENTS[i], NUM_RUNS[i], time_per_test);
  }
  bf << "\n";
}

#include <openvoronoi/voronoidiagram.hpp>
void run_ovd_test(const std::vector<double> &running_times) {
    boost::mt19937 gen(RANDOM_SEED);
    double minimum_coordinate = std::numeric_limits<uint32_t>::min();
    double maximum_coordinate = std::numeric_limits<uint32_t>::max();
    bf << "OpenVoronoi:\n";
    for (int i = 0; i < NUM_TESTS; ++i) {
        timer.restart();
        for (int j = 0; j < NUM_RUNS[i]; ++j) {
            SSD_POLYGON ssd;
            for (int k = 0; k < NUM_SEGMENTS[i]; ++k) {
                int32 x1 = gen();
                int32 y1 = gen();
                int32 dx = (gen() & 1023) + 1;
                int32 dy = (gen() & 1023) + 1;
                ssd.insert(SEGMENT_POLYGON(POINT_POLYGON(x1, y1), POINT_POLYGON(x1 + dx, y1 + dy)));
            }
            ssd.clean();
            int bins = static_cast<int>( sqrt(NUM_SEGMENTS[i] ) );
            ovd::VoronoiDiagram vd(1, bins );
            vd.set_silent(true); // suppress warning messages/output
            typedef std::pair<int,int> segment_ids;
            typedef std::vector<segment_ids> polygon_ids;
            polygon_ids ids;
            for (SSD_POLYGON::iterator_type it = ssd.begin(); it != ssd.end(); ++it) {
                double x1 = static_cast<double>( it->low().x() );
                double y1 = static_cast<double>( it->low().y() );
                double x2 = static_cast<double>( it->high().x() );
                double y2 = static_cast<double>( it->high().y() );

                // For OpenVoronoi we scale coordinates so they fit within a unit circle.
                // e.g. a box centered at (0,0) with side-length 0.8
                x1 -= (maximum_coordinate-minimum_coordinate); // center around x=0
                y1 -= (maximum_coordinate-minimum_coordinate); // center around y=0
                x1 *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
                y1 *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
                x2 -= (maximum_coordinate-minimum_coordinate); // center around x=0
                y2 -= (maximum_coordinate-minimum_coordinate); // center around y=0
                x2 *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
                y2 *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
                //std::cout << " point site : " << x1 << " , " << y1 << "\n" << std::flush;
                segment_ids seg;
                seg.first = vd.insert_point_site( ovd::Point( x1, y1) );
                seg.second = vd.insert_point_site( ovd::Point( x2, y2) );
                ids.push_back(seg);           
            }
            // now insert line-segments
            for (unsigned int n=0;n<ids.size();++n)
                vd.insert_line_site( ids[n].first, ids[n].second );
        }
        double time_per_test = (timer.elapsed() - running_times[i]) / NUM_RUNS[i];
        format_line(NUM_SEGMENTS[i], NUM_RUNS[i], time_per_test);
    }
    bf << "\n";
}


int main() {
  bf << std::setiosflags(std::ios::right | std::ios::fixed) << std::setprecision(6);
  std::vector<double> running_times = get_intersection_runtime();
  run_ovd_test(running_times);
  run_voronoi_test(running_times);
  run_cgal_test(running_times);
  bf.close();
  return 0;
}
