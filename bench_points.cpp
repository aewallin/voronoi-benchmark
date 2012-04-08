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
typedef boost::polygon::default_voronoi_builder VB_BOOST;
typedef boost::polygon::voronoi_diagram<double> VD_BOOST;

// Includes for CGAL library.
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
typedef CGAL::Quotient<CGAL::MP_Float> ENT;
typedef CGAL::Simple_cartesian<double> CK;
typedef CGAL::Simple_cartesian<ENT> EK;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<CK, CGAL::Field_with_sqrt_tag, EK, CGAL::Field_tag> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt> SDT_CGAL;
typedef SDT_CGAL::Point_2 Point_CGAL;

std::ofstream bf("benchmark_points.csv", std::ios_base::out | std::ios_base::app);
boost::timer timer;

// print a row of results to the file
void format_line(std::string lib_name, int num_points, int num_tests, double total_time, double time_per_test) {
  bf <<  lib_name << ",";
  bf << num_points << ",";
  bf << num_tests << ",";
  bf << time_per_test << ",";
  bf <<  std::endl;
    std::cout << "| " << std::setw(16) << lib_name << " ";
    std::cout << "| " << std::setw(16) << num_points << " ";
    std::cout << "| " << std::setw(15) << num_tests << " ";
    std::cout << "| " << std::setw(17) << time_per_test << " ";
    std::cout << "| " << std::setw(17) << total_time << " ";
    std::cout << "| " << std::setw(17) << 1e6*time_per_test/(num_points*log(num_points)/log(2)) << " ";
    std::cout << "|" << std::endl << std::flush;
}

void run_boost_test(int RANDOM_SEED, int NUM_TESTS, std::vector<int> NUM_POINTS, std::vector<int> NUM_RUNS) {
  boost::mt19937 gen(RANDOM_SEED);
  //bf << "Boost.Polygon Voronoi of Points:\n";
  for (int i = 0; i < NUM_TESTS; ++i) {
    timer.restart();
    for (int j = 0; j < NUM_RUNS[i]; ++j) {
      VB_BOOST vb;
      VD_BOOST vd;
      for (int k = 0; k < NUM_POINTS[i]; ++k) {
        vb.insert_point(static_cast<int32>(gen()), static_cast<int32>(gen()));
      }
      vb.construct(&vd);
    }
    double total_time = timer.elapsed();
    double time_per_test = total_time / NUM_RUNS[i];
    format_line("Boost.Polygon", NUM_POINTS[i], NUM_RUNS[i], total_time, time_per_test);
  }
  //bf << "\n";
}


void run_cgal_test(int RANDOM_SEED, int NUM_TESTS, std::vector<int> NUM_POINTS, std::vector<int> NUM_RUNS) {
  boost::mt19937 gen(RANDOM_SEED);
  //bf << "CGAL Triangulation of Points:\n";
  for (int i = 0; i < NUM_TESTS; ++i) {
    timer.restart();
    for (int j = 0; j < NUM_RUNS[i]; ++j) {
      SDT_CGAL dt;
      for (int k = 0; k < NUM_POINTS[i]; ++k) {
        dt.insert(Point_CGAL(static_cast<int32>(gen()), static_cast<int32>(gen())));
      }
    }
    double total_time = timer.elapsed();
    double time_per_test = total_time / NUM_RUNS[i];
    format_line("CGAL", NUM_POINTS[i], NUM_RUNS[i], total_time, time_per_test);
  }
  //bf << "\n";
}

#include <openvoronoi/voronoidiagram.hpp>
#include <openvoronoi/version.hpp>

void run_ovd_test(int RANDOM_SEED, int NUM_TESTS, std::vector<int> NUM_POINTS, std::vector<int> NUM_RUNS) {
    boost::mt19937 gen(RANDOM_SEED);
    double minimum_coordinate = std::numeric_limits<uint32_t>::min();
    double maximum_coordinate = std::numeric_limits<uint32_t>::max();
    //std::cout << "OpenVoronoi " << ovd::version() << " " << ovd::build_type() << "\n";
    //bf << "OpenVoronoi:\n";
  for (int i = 0; i < NUM_TESTS; ++i) {
    timer.restart();
    for (int j = 0; j < NUM_RUNS[i]; ++j) {
        int bins = static_cast<int>( sqrt(NUM_POINTS[i] ) );
        ovd::VoronoiDiagram vd(1, bins ) ;
        for (int k = 0; k < NUM_POINTS[i]; ++k) {
            double x = static_cast<double>(gen());
            double y = static_cast<double>(gen());

            // For OpenVoronoi we scale coordinates so they fit within a unit circle.
            // e.g. a box centered at (0,0) with side-length 0.8
            x -= (maximum_coordinate-minimum_coordinate); // center around x=0
            y -= (maximum_coordinate-minimum_coordinate); // center around y=0
            x *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
            y *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
            vd.insert_point_site( ovd::Point( x, y) );
        }
    }
    double total_time = timer.elapsed();
    double time_per_test = total_time / NUM_RUNS[i];
    format_line("OpenVoronoi", NUM_POINTS[i], NUM_RUNS[i], total_time, time_per_test);
  }
  //bf << "\n";
}

int main() {
    const int RANDOM_SEED = 27;
    const int max_exponent = 40;
    std::vector<int> NUM_POINTS;
    std::vector<int> NUM_RUNS;

    for (int m=10;m<max_exponent;m++) {
        NUM_POINTS.push_back( (int)( pow(2,0.5*m) ) ); // this nicely spaced points on log-axis
        if (m<20)
            NUM_RUNS.push_back(100);
        else if (m<30)
            NUM_RUNS.push_back(10);
        else
            NUM_RUNS.push_back(2);
    }
    int NUM_TESTS = NUM_POINTS.size();

    bf << std::setiosflags(std::ios::right | std::ios::fixed) << std::setprecision(6);
    run_ovd_test(RANDOM_SEED, NUM_TESTS, NUM_POINTS, NUM_RUNS);
    run_boost_test(RANDOM_SEED, NUM_TESTS, NUM_POINTS, NUM_RUNS);
    run_cgal_test(RANDOM_SEED, NUM_TESTS, NUM_POINTS, NUM_RUNS);

    bf.close();
    return 0;
}
