//
// Created by rohit on 23/03/20.
//
/*
 * Refrences: https://people.cs.clemson.edu/~dhouse/courses/401/notes/affines-matrices.pdf
 * I referenced above PDF for implementation.
 */
#ifndef ASTRONOMY_AFFINE_TRANSFORMATION_HPP
#define ASTRONOMY_AFFINE_TRANSFORMATION_HPP

#include <boost/astronomy/coordinate/representation.hpp>
#include <boost/astronomy/coordinate/arithmetic.hpp>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/plane_angle.hpp>
#include <boost/units/systems/si/prefixes.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/angle/degrees.hpp>

#include <boost/numeric/ublas/matrix.hpp>


namespace  bus = boost::units::si;
namespace  bac = boost::astronomy::coordinate;
namespace  bu = boost::units;
namespace  bnu = boost::numeric::ublas;

using  namespace std;
using namespace boost::numeric::ublas;


namespace boost{ namespace astronomy { namespace coordinate{

template
<
typename ElementType = double,
typename Representation = bac::cartesian_representation<
        ElementType,bu::quantity<bus::length>,
        bu::quantity<bus::length>,
        bu::quantity<bus::length>>
>
struct affine_transformation{
private:
    matrix<ElementType> shear = matrix<ElementType>(3, 3);
    matrix<ElementType> scale = matrix<ElementType>(3, 3);
    matrix<ElementType> rotate= matrix<ElementType>(3, 3);
    Representation translation_vec;
public:

    //! construct object with no parameters
    affine_transformation(){
        this->shear = identity_matrix<ElementType>(3,3);
        this->scale = identity_matrix<ElementType>(3,3);
        this->rotate= identity_matrix<ElementType>(3,3);
        this->translation_vec = Representation(0.0*bus::meter, 0.0*bus::meter, 0.0*bus::meter);
    }

    //! construct object with give three matrices shear, scale, rotate
    affine_transformation
    (
            matrix<ElementType> const& shear,
            matrix<ElementType> const& scale,
            matrix<ElementType> const& rotate
    )
    {
        this->shear = shear;
        this->scale = scale;
        this->rotate= rotate;
        this->translation_vec = Representation(0.0*bus::meter, 0.0*bus::meter, 0.0*bus::meter);
    }

    //!setter function to set shear_matrix with give shear matrix
    void set_shear(matrix<ElementType> shear_matrix){
        this->shear = shear_matrix;
    }

    //!setter function to set shear_matrix with give shear parameters
    void set_shear
            (
                    ElementType hxy,
                    ElementType hxz,
                    ElementType hyx,
                    ElementType hyz,
                    ElementType hzx,
                    ElementType hzy
            )
    {
        this->shear(0,1) = hxy;
        this->shear(0,2) = hxz;
        this->shear(1,0) = hyx;
        this->shear(1,2) = hyz;
        this->shear(1,0) = hzx;
        this->shear(1,0) = hzy;
    }

    //!setter function to set rotation_matrix with give rotate matrix
    void set_rotation(matrix<ElementType> const& rotation_matrix){
        this->rotate = rotation_matrix;
    }

    //!setter function to set rotation_matrix with give rotation parameters(in terms of quantities)
    template
    <typename q1, typename q2, typename q3>
    void set_rotate(bu::quantity<q1> x_angle, bu::quantity<q2> y_angle, bu::quantity<q3> z_angle){
        ElementType x_ang = static_cast<bu::quantity<bus::plane_angle>>(x_angle).value();
        ElementType y_ang = static_cast<bu::quantity<bus::plane_angle>>(y_angle).value();
        ElementType z_ang = static_cast<bu::quantity<bus::plane_angle>>(z_angle).value();

        matrix<ElementType> x_matrix;
        matrix<ElementType> y_matrix;
        matrix<ElementType> z_matrix;

        x_matrix = identity_matrix<ElementType>(3,3);
        y_matrix = identity_matrix<ElementType>(3,3);
        z_matrix = identity_matrix<ElementType>(3,3);

        x_matrix(1,1) = std::cos(x_ang);
        x_matrix(1,2) = -std::sin(x_ang);
        x_matrix(2,1) = std::sin(x_ang);
        x_matrix(2,2) = std::cos(x_ang);

        y_matrix(0,0) = std::cos(y_ang);
        y_matrix(0,2) = std::sin(y_ang);
        y_matrix(2,0) = -std::sin(y_ang);
        y_matrix(2,2) = std::cos(y_ang);

        z_matrix(0,0) = std::cos(z_ang);
        z_matrix(0,1) = -std::sin(z_ang);
        z_matrix(1,0) = std::sin(z_ang);
        z_matrix(1,1) = std::cos(z_ang);

        this->rotate = prod(x_matrix,y_matrix);
        this->rotate = prod(this->rotate,z_matrix);
    }

    //!setter function to set rotation_matrix with give rotation parameters
    void set_rotate(ElementType x_ang,ElementType y_ang,ElementType z_ang){
        matrix<ElementType> x_matrix = matrix<ElementType>(3,3);
        matrix<ElementType> y_matrix = matrix<ElementType>(3,3);
        matrix<ElementType> z_matrix = matrix<ElementType>(3,3);

        x_matrix = identity_matrix<ElementType>(3,3);
        y_matrix = identity_matrix<ElementType>(3,3);
        z_matrix = identity_matrix<ElementType>(3,3);

        x_matrix(1,1) = std::cos(x_ang);
        x_matrix(1,2) = -std::sin(x_ang);
        x_matrix(2,1) = std::sin(x_ang);
        x_matrix(2,2) = std::cos(x_ang);

        y_matrix(0,0) = std::cos(y_ang);
        y_matrix(0,2) = std::sin(y_ang);
        y_matrix(2,0) = -std::sin(y_ang);
        y_matrix(2,2) = std::cos(y_ang);

        z_matrix(0,0) = std::cos(z_ang);
        z_matrix(0,1) = -std::sin(z_ang);
        z_matrix(1,0) = std::sin(z_ang);
        z_matrix(1,1) = std::cos(z_ang);

        this->rotate = prod(x_matrix,y_matrix);
        this->rotate = prod(this->rotate,z_matrix);
    }

    //!setter function to set scale_matrix with give scale matrix
    void set_scale(matrix<ElementType> const& scale_matrix){
        this->scale = scale_matrix;
    }

    //!setter function to set scale_matrix with give scale parameters
    void set_scale(ElementType sx,ElementType sy,ElementType sz){
        this->scale(0,0) = sx;
        this->scale(1,1) = sy;
        this->scale(2,2) = sz;
    }

    //! setter function to set translation vector
    void set_translation_vector(Representation const & vec){
        this->translation_vec = vec;
    }

    //!getter function to get scale matrix
    matrix<ElementType> get_scale(){
        return this->scale;
    }

    //!getter function to get shear matrix
    matrix<ElementType> get_shear(){
        return this->shear;
    }

    //!getter function to get rotate matrix
    matrix<ElementType> get_rotate(){
        return this->rotate;
    }

    //!getter function to get translation vector
    Representation get_translation_vec(){
        return this->translation_vec;
    }

    //! translate a representation point
    template <typename representation>
    Representation Translate(representation const& point){
        return point+this->translation_vec;
    }


    //! Rotate a representation point
    template <typename representation>
    Representation Rotate(representation const& point){
        auto tempRep = boost::astronomy::coordinate::make_cartesian_representation(point);
        boost::astronomy::coordinate::cartesian_representation
        <
            ElementType,
            typename representation::quantity3,
            typename representation::quantity3,
            typename representation::quantity3
        >tempRep1
        (
                static_cast<typename representation::quantity3>(tempRep.get_x()),
                static_cast<typename representation::quantity3>(tempRep.get_y()),
                static_cast<typename representation::quantity3>(tempRep.get_z())
        );

        auto gp = tempRep1.get_point();

        bnu::vector<ElementType> tempVec(3);

        tempVec(0) = boost::geometry::get<0>(gp);
        tempVec(1) = boost::geometry::get<1>(gp);
        tempVec(2) = boost::geometry::get<2>(gp);

        auto res = prod(this->rotate,tempVec);

        boost::geometry::set<0>(gp,res(0));
        boost::geometry::set<1>(gp,res(1));
        boost::geometry::set<2>(gp,res(2));
        return Representation(gp);
    }

    //! Apply Shear to a representation point
    template <typename representation>
    Representation Shear(representation const& point){
        auto tempRep = boost::astronomy::coordinate::make_cartesian_representation(point);
        boost::astronomy::coordinate::cartesian_representation
                <
                        ElementType,
                        typename representation::quantity3,
                        typename representation::quantity3,
                        typename representation::quantity3
                >tempRep1
                (
                        static_cast<typename representation::quantity3>(tempRep.get_x()),
                        static_cast<typename representation::quantity3>(tempRep.get_y()),
                        static_cast<typename representation::quantity3>(tempRep.get_z())
                );

        auto gp = tempRep1.get_point();

        bnu::vector<ElementType> tempVec(3);

        tempVec(0) = boost::geometry::get<0>(gp);
        tempVec(1) = boost::geometry::get<1>(gp);
        tempVec(2) = boost::geometry::get<2>(gp);

        auto res = prod(this->shear,tempVec);

        boost::geometry::set<0>(gp,res(0));
        boost::geometry::set<1>(gp,res(1));
        boost::geometry::set<2>(gp,res(2));
        return Representation(gp);
    }
    //! Scale a representation point
    template <typename representation>
    Representation Scale(representation const& point){
        auto tempRep = boost::astronomy::coordinate::make_cartesian_representation(point);
        boost::astronomy::coordinate::cartesian_representation
                <
                        ElementType,
                        typename representation::quantity3,
                        typename representation::quantity3,
                        typename representation::quantity3
                >tempRep1
                (
                        static_cast<typename representation::quantity3>(tempRep.get_x()),
                        static_cast<typename representation::quantity3>(tempRep.get_y()),
                        static_cast<typename representation::quantity3>(tempRep.get_z())
                );

        auto gp = tempRep1.get_point();

        bnu::vector<ElementType> tempVec(3);

        tempVec(0) = boost::geometry::get<0>(gp);
        tempVec(1) = boost::geometry::get<1>(gp);
        tempVec(2) = boost::geometry::get<2>(gp);

        auto res = prod(this->scale,tempVec);

        boost::geometry::set<0>(gp,res(0));
        boost::geometry::set<1>(gp,res(1));
        boost::geometry::set<2>(gp,res(2));
        return Representation(gp);
    }

    //! get a transformed coordinate of a representation point
    template
    <typename representation>
    representation Transform(representation const& point){
        return this->Translate(this->Shear(this->Rotate(this->Scale(point))));
    }
};

}}}//boost::astronomy::coordinate

#endif //ASTRONOMY_AFFINE_TRANSFORMATION_HPP
