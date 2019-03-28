/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <list>

#include "eckit/container/KDTree.h"
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/UnitSphere.h"

namespace fckit {

extern "C" {

double fckit__sphere_distance(const double & lonA, const double & latA, const double & lonB, const double & latB) {
  const eckit::geometry::Point2 A(lonA, latA);
  const eckit::geometry::Point2 B(lonB, latB);
  return eckit::geometry::UnitSphere::distance(A, B);
}

void fckit__sphere_lonlat2xyz(const double & lon, const double & lat, double & x, double & y, double & z) {
  const eckit::geometry::Point2 lonlat(lon, lat);
  eckit::geometry::Point3 xyz;
  eckit::geometry::UnitSphere::convertSphericalToCartesian(lonlat, xyz);
  x = xyz[0];
  y = xyz[1];
  z = xyz[2];
}

void fckit__sphere_xyz2lonlat(const double & x, const double & y, const double & z, double & lon, double & lat) {
  const eckit::geometry::Point3 xyz(x, y, z);
  eckit::geometry::Point2 lonlat;
  eckit::geometry::UnitSphere::convertCartesianToSpherical(xyz, lonlat);
  lon = lonlat[0];
  lat = lonlat[1];
}

struct TreeTrait {
    typedef eckit::geometry::Point3 Point;
    typedef double                  Payload;
};

typedef eckit::KDTreeMemory<TreeTrait> KDTree;

KDTree* fckit__kdtree_create(const int & n, const double * lon, const double * lat) {
    // Define tree
    KDTree * kd = new KDTree();

    // Define points list
    typedef KDTree::PointType Point;
    std::vector<KDTree::Value> points;
    for (unsigned int i = 0; i < n; i++) {
        eckit::geometry::Point2 lonlat(lon[i], lat[i]);
        Point xyz = Point();
        eckit::geometry::UnitSphere::convertSphericalToCartesian(lonlat, xyz);
        double index = static_cast<double>(i);
        KDTree::Value v(xyz, index);
        points.push_back(v);
    }

    // Build and return KDTree
    kd->build(points.begin(), points.end());
    return kd;
}

void fckit__kdtree_destroy (KDTree * kd) {
    delete kd;
}

void fckit__kdtree_k_nearest_neighbors(KDTree * kd, const double & lon, const double & lat, const int & nn, int * nn_index) {
    // Define central point
    eckit::geometry::Point2 lonlat(lon, lat);
    eckit::geometry::Point3 xyz;
    eckit::geometry::UnitSphere::convertSphericalToCartesian(lonlat, xyz);

    // Find nn nearest neighbors
    KDTree::NodeList list = kd->kNearestNeighbours(xyz, nn);

    // Copy nearest neighbors index
    int i = 0;
    for (KDTree::NodeList::iterator it = list.begin(); it != list.end(); ++it)
    {
       nn_index[i] = static_cast<int>(it->payload());
       i++;
    }
    ASSERT(i==nn);
}

void fckit__kdtree_find_in_sphere(KDTree * kd, const double & lon, const double & lat, const double & r, int & nn) {
    // Define central point
    eckit::geometry::Point2 lonlat(lon, lat);
    eckit::geometry::Point3 xyz;
    eckit::geometry::UnitSphere::convertSphericalToCartesian(lonlat, xyz);

    // Convert radius on unit sphere to chord
    double chord = 2.0*sin(0.5*std::min(r, M_PI));

    // Find nearest neighbors within radius r
    KDTree::NodeList list = kd->findInSphere(xyz, chord);

    // Copy number of neighbors
    nn = list.size();
}

} // extern "C"

} // namespace fckit
