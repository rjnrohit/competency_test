//
// Created by rohit on 23/03/20.
//

#include <iostream>
#include "affine_transformation.hpp"
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/systems/angle/degrees.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace boost::astronomy::coordinate;
using namespace boost::units;
using namespace boost::units::si;
using namespace boost::units::cgs;
using namespace boost::geometry;
using namespace boost::numeric::ublas;
namespace  bud = boost::units::degree;

int main(){
    affine_transformation at;
    cout << at.get_scale() << endl;
    cout << at.get_rotate() << endl;
    cout << at.get_shear() << endl;
    cout << at.get_translation_vec() << endl;

    auto point = make_cartesian_representation(1.0*meter,23.0*meter,-2.0*meter);

    at.set_translation_vector(point);
    at.set_scale(12,1,3);
    at.set_rotate(23.1,37,53);
    at.set_shear(1.3,4.5,6.7,7.8,9.1,1.2);

    cout << at.get_scale() << endl;
    cout << at.get_rotate() << endl;
    cout << at.get_shear() << endl;
    cout << at.get_translation_vec() << endl;

    at.set_translation_vector(point);
    at.set_scale(12,1,3);
    at.set_rotate(23.1* bud::degree,37.2*bud::degree ,53.1*bud::degree);
    at.set_shear(1.3,4.5,6.7,7.8,9.1,1.2);

    cout << at.get_scale() << endl;
    cout << at.get_rotate() << endl;
    cout << at.get_shear() << endl;
    cout << at.get_translation_vec() << endl;

    cout << at.Rotate(point) << endl;
    cout << at.Transform(point) << endl;

}