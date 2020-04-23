**Affine Transformation**

****GSOC2020 programming competency_test****

This is solution of competency test mentioned at https://github.com/boostorg/wiki/wiki/Google-Summer-of-Code%3A-2020#programming-competency-test-4.

Using this we can transform any representation by following ways:
1. Translate 
2. Shear 
3. Scale
4. Rotate
5. Transform(To apply all above three methods at once)

We can transform all representations of coordinates i.e Cartesian Coordinates, Spherical Coordinates and Spherical Equatorial Coordinates.
My affine transformation class expects One translation vector(In any representation) and a 3x3 matrix to apply Affine Transformation.
A Demo program can be seen as follows:

```
//Pick a Coordinate
auto point = make_spherical_representation(1.43*radian, 1.58*radian, 45.4*meter);

//Intialize an affine_transformation class
Affine_Transformation<decltype(point)> at;

// set translation vector
auto translate = make_spherical_representation(0.63*radian,0,89*radian, 21*meter);
at.set_translation_vector(translate);

//set shear matrix
//                  [1.0,-1.0,6.0]
//  shear matrix is [2.0,1.0, 5.0]
//                  [3.0,8.0 ,1.0]   
matrix<double> shear(3,3)= {{1.0,-1.0,6.0},{2.0,1.0,5.0},{3.0,8.0,1.0}};
at.set_shear(shear);

//set scale  matrix
//                  [0.23,0.0,0.0]
//  scale matrix is [0.0,0.87,0.0]
//                  [0.0,0.0,0.13]   
matrix<double> scale(3,3)= {{0.23,0.0,0.0},{0.0,0.87,0.0},{0.0,0.0,0.13}};
at.set_scale(scale);

// set rotate parameters
//x_angle = 23*radian, y_anlge = 24*radian, z_angle = 36*degree
at.set_rotate_parameters(x_angle, y_angle,z_angle);

// Apply Translate method to translate above coordinate
auto translate_coordinate = at.translate(point);

// Apply Shear method to above coordinate
auto shear_coordinate = at.shear(translate_coordinate);

//Apply rotation method
auto rotate_coordinate = at.rotate(shear_coordinate);

//Apply Scaling
auto scale_coordinate = at.scale(rotate_coordinate);

//Apply affine transformation all method at once

auto transformed = at.transform(point);
```

