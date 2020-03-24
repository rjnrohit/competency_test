
#include <boost/units/io.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/astronomy/coordinate/cartesian_representation.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/astronomy/coordinate/io.hpp>

#define ROW_SIZE 3

/*
 * Refrences: https://people.cs.clemson.edu/~dhouse/courses/401/notes/affines-matrices.pdf
 * I referenced above PDF for implementation.
 */

namespace  bus = boost::units::si;
namespace  bac = boost::astronomy::coordinate;
namespace  bu = boost::units;
namespace  bnu = boost::numeric::ublas;
using namespace bnu;
using namespace bus;
using namespace bu;
using namespace bac;

typedef  bac::cartesian_representation<double, quantity<si::length>, quantity<si::length>, quantity<si::length>> cord_rep;

template
<
typename  ElementType = double
>
struct affine_matrix
{
public:
    /*
     * affine matrix is a class. It contains
     * three matrices scale, shear and rotation
     * which is necessary for transformation
     */
    bnu::matrix<ElementType>scale ;
    bnu::matrix<ElementType>shear ;
    bnu::matrix<ElementType>rotation ;

    //! constructor with no parameters
    //! resizing of ublas::matrix is necessary
    affine_matrix(){
        this->scale.resize(ROW_SIZE,ROW_SIZE);
        this->rotation.resize(ROW_SIZE,ROW_SIZE);
        this->shear.resize(ROW_SIZE,ROW_SIZE);
        for(int i=0;i<ROW_SIZE;++i){
            for(int j=0;j<ROW_SIZE;++j){
                this->scale(i,j) = (ElementType)(i==j);
                this->shear(i,j) = (ElementType)(i==j);
                this->rotation(i,j) = (ElementType)(i==j);
            }
        }

    }
    //! constructor with given three matrices as parameters
    affine_matrix(matrix <ElementType> const&shear_matrix,matrix <ElementType> const&scale_matrix,matrix <ElementType>const& rotation_matrix){
        this->scale.resize(ROW_SIZE,ROW_SIZE);
        this->rotation.resize(ROW_SIZE,ROW_SIZE);
        this->shear.resize(ROW_SIZE,ROW_SIZE);
        this->scale = scale_matrix;
        this->shear = shear_matrix;
        this->rotation = rotation_matrix;
    }

    void set_shear(matrix<ElementType>const& shear_matrix){
        this->shear = shear_matrix;
    }
    /*
     * setter function for shear matrix with given matrix element parameter
     */
    void set_shear(ElementType hxy,ElementType hxz,ElementType hyx, ElementType hyz,ElementType hzx, ElementType hzy){
        this->shear(0,1) = hxy;
        this->shear(0,2) = hxz;
        this->shear(1,0) = hyz;
        this->shear(1,2) = hyz;
        this->shear(1,0) = hzx;
        this->shear(1,0) = hzy;
    }
    void set_rotation(matrix<ElementType> const& rotation_matrix){
        this->rotation = rotation_matrix;
    }

    void set_scale(matrix<ElementType> const& scale_matrix){
        this->scale = scale_matrix;
    }

    /*
     * setter function for scale matrix with given matrix elements
     */
    void set_scale(ElementType sx,ElementType sy,ElementType sz){
        this->scale(0,0) = sx;
        this->scale(1,1) = sy;
        this->scale(2,2) = sz;
    }


    //! get_affine_matrix returns affine_matrix for affine transformation
    //! It multiply three matrices scale, shear and rotation
    matrix<ElementType> get_affine_matrix() const{

        auto temp =  prod(shear,scale);
        matrix<ElementType> result;
        result.resize(ROW_SIZE,ROW_SIZE);

        for(int i=0;i<3;++i){
            for(int j=0;j<3;++j){
                result(i,j) =0;
                for(int k=0;k<3;++k){
                    result(i,j) += temp(k,j)*rotation(i,k);
                }
            }
        }
        return result;
    }
};


template
<typename  ElementType = double>
struct affine_transformation{
public:
    /*
     * This affine_transformation class accepts an cartesina_representation parameter
     * and return tranformed coordinate.
     * It has affine attribute which is affine transformation matrix and
     * a translation coordinate where it contains DELTA_X, DELTA_Y,DELTA_Z
     */

    matrix<ElementType> affine;
    affine_matrix<ElementType> affine_matrix_object;
    cord_rep translation_cord;

    //! constructor without any parameter
    affine_transformation(){
        this->affine = affine_matrix_object.get_affine_matrix();
        translation_cord.set_x_y_z(0*meter,0*meter,0*meter);

    }

    //! construct object with affineMatrix parameter
    explicit  affine_transformation(affine_matrix<ElementType>const& affineMatrix){
        this->affine_matrix_object = affineMatrix;
        this->affine = affineMatrix.get_affine_matrix();
        translation_cord.set_x_y_z(0*meter,0*meter,0*meter);
    }
    //! construct object with translation_vector parameter
    explicit  affine_transformation(cord_rep const&translation_vector){
        this->translation_cord  = translation_vector;
    }
    //! construct object with shear Matrices shear, scale and rotation
    affine_transformation(matrix<ElementType>const& shear,matrix<ElementType>const& scale,matrix<ElementType>const& rotation){
        this->affine_matrix_object.set_rotation(rotation);
        this->affine_matrix_object.set_shear(shear);
        this->affine_matrix_object.set_scale(scale);
        this->affine = affine_matrix_object.get_affine_matrix();
        translation_cord.set_x_y_z(0*meter,0*meter,0*meter);
    }

    void set_affine_matrix_object(affine_matrix<ElementType> const& object){
        this->affine_matrix_object = object;
        this->affine = object.get_affine_matrix();
    }

    void set_translation_vector(cord_rep const& trans_vec){
        this->translation_cord  = trans_vec;
    }


    //! returns transformed cartesian_representation of a coordinate
    cord_rep get_transformed_coordinate(cord_rep const& test_cord){
        vector<ElementType> test_cord_params (3);
        vector<ElementType> translation_vec (3);

        test_cord_params(0) = test_cord.get_x().value();
        test_cord_params(1) = test_cord.get_y().value();
        test_cord_params(2) = test_cord.get_z().value();

        translation_vec(0) = translation_cord.get_x().value();
        translation_vec(1) = translation_cord.get_y().value();
        translation_vec(2) = translation_cord.get_z().value();

        auto result = prod(affine,test_cord_params) + translation_vec;

        return make_cartesian_representation(result(0)*meter,result(1)*meter,result(2)*meter);
    }
};



int main(){
    std::cout << "Running affine_transformation" << std::endl;
    //pick a point in cartesian coordinate
    cord_rep point = make_cartesian_representation(22*meter, 13*meter, 19*meter);

    std::cout <<"Point: " <<point << std::endl;
    //take a translation_coordinate
    cord_rep  translation_vec = make_cartesian_representation(12*meter, 13*meter,15*meter);
    std::cout <<"translation vec:" <<  translation_vec << std::endl;
    //create a affine_matrix object
    affine_matrix <double> am;
    // printing affine_matrix attributes

    std:: cout << "affine_ matrix\n" <<am.get_affine_matrix() <<std::endl;

    //create a object of affine_transformation class
    affine_transformation<double> at;

    //initialization of shear matrix
    //here shear is
    /*
        [
            1,-1, 6,
            2, 1, 5,
            3, 8, 1
        ]
    */
    matrix <double> shear(ROW_SIZE,ROW_SIZE);
    shear(0,0) = 1;
    shear(0,1) = -1;
    shear(0,2) = 6;
    shear(1,0) = 2;
    shear(1,1) = 1;
    shear(1,2) = 5;
    shear(2,0) = 3;
    shear(2,1) = 8;
    shear(2,2) = 1;

    am.set_shear(shear);

    //initialization of scale matrix
    //here scale is
    /*
        [
            2, 0, 0,
            0, 3, 0,
            0, 0, 4
        ]
    */
    matrix <double> scale(ROW_SIZE,ROW_SIZE);
    scale(0,0) = 2;
    scale(0,1) = 0;
    scale(0,2) = 0;
    scale(1,0) = 0;
    scale(1,1) = 3;
    scale(1,2) = 0;
    scale(2,0) = 0;
    scale(2,1) = 0;
    scale(2,2) = 4;

    am.set_scale(scale);
    //initialization of rotation matrix
    //here rotation is
    /*
        [
            0.2,-0.43, 0.6,
            0.24, 0.1, 0.0,
            0.03, 0.8, 0.1
        ]
    */
    matrix <double> rotation(ROW_SIZE,ROW_SIZE);
    rotation(0,0) = 0.2;
    rotation(0,1) = -0.43;
    rotation(0,2) = 0.6;
    rotation(1,0) = 0.24;
    rotation(1,1) = 0.1;
    rotation(1,2) = 0.0;
    rotation(2,0) = 0.03;
    rotation(2,1) = 0.8;
    rotation(2,2) = 0.1;

    am.set_rotation(rotation);

    std::cout << "am shear\n" << am.shear << std::endl;
    std::cout << "am scale\n" << am.scale << std::endl;
    std::cout << "am rotation\n" << am.rotation << std::endl;

    std::cout << "am affine_matrix\n" << am.get_affine_matrix() << std::endl;

    //set affine matrix object
    at.set_affine_matrix_object(am);

    std::cout << "at affine\n" << at.affine << std::endl;

    at.set_translation_vector(translation_vec);
    std::cout << "at trans_cord\n" << at.translation_cord << std::endl;

    //Printing transformed Point
    std::cout << at.get_transformed_coordinate(point) << std::endl;
    return 0;
}


